

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include <boost/format.hpp>

#include <cstdio>

#include "constants.hpp"
#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "utils-string.hpp"
#include "utils-streamio.hpp"
#include "utils-errors.hpp"

#include "elem-iacs.hpp"
#include "specs-fit-prop-pot.hpp"
#include "compound.hpp"
#include "physconst.hpp"

#include "compoundfit.hpp"


using namespace std;
using namespace utils;
using namespace constants;
using boost::format;



CompoundStructureFit::CompoundStructureFit()
  :
  CompoundStructure()
{
  prop_use    = CompoundPropertiesUse();
  prop_readin = CompoundProperties();
  prop_pred   = CompoundProperties();

  mds_specs = MDSettings();

  use_u = CompoundPropertiesUseUncertainties();
  use_w = CompoundPropertiesUseWeights();
  prop_u = CompoundPropertiesUncertainties();
  prop_w = CompoundPropertiesWeights();
}




// Make sure we are not using the same properties in different
// disguises, since this will lead to problems in some fitting
// algorithms:


void CompoundStructureFit::check_and_fix_uses(){
  prop_use.check_and_fix();

}

int CompoundStructureFit::NData(){
  int N=0;

  if (prop_use.a) N++;
  if (prop_use.b) N++;
  if (prop_use.c) N++;
  if (prop_use.bpa) N++;
  if (prop_use.cpa) N++;
  if (prop_use.r0) N++;
  if (prop_use.angle_ab) N++;
  if (prop_use.angle_ac) N++;
  if (prop_use.angle_bc) N++;
  if (prop_use.Vatom) N++;
  if (prop_use.Ecoh) N++;
  if (prop_use.Emix) N++;
  if (prop_use.B) N++;
  if (prop_use.Bp) N++;
  for (int k=0; k<6; ++k)
    for (int p=0; p<6; ++p)
      if (prop_use.C.elem(k,p)) N++;
  if (prop_use.Fmax) N++;
  if (prop_use.Pmax) N++;
  if (prop_use.displmax) N++;
  if (prop_use.frc){
    N += 3.0 * basis_elems.size();
  }

  return N;
}





void CompoundStructureFit::read_forces(void){
  ifstream fp;
  string line;
  vector<string> args;
  istringstream strbuf;
  int ns, i;
  double fx, fy, fz, td1, td2, td3;
  

  prop_use.frc = true;

  // Allocate space for forces:
  int nb = basis_elems.size();

  prop_readin.frc.resize(nb);
  for (i=0; i<nb; ++i) prop_readin.frc[i].resize(3);

  // ***************************************************************************
  // Allocate space for predicted forces. Do this now when forces are read,
  // so we won't have to wonder later if they are allocated or not.
  // ***************************************************************************
  prop_pred.frc.resize(nb);
  for (i=0; i<nb; ++i) prop_pred.frc[i].resize(3);

  if (use_u.frc){
    prop_u.frc.resize(nb);
    for (i=0; i<nb; ++i)
      prop_u.frc[i].resize(3);
  }
  else {
    prop_w.frc.resize(nb);
    for (i=0; i<nb; ++i)
      prop_w.frc[i].resize(3);
  }


  // Prepare file:
  fp.open(filename_frc.c_str());
  if (!fp)
    aborterror("Error: Could not find/open forces file " +
	       filename + " for compound " +
	       name + ". Exiting.");
  
  // Read file for forces:
  int iat = 0, iline = 0;
  while (true){

    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");

    if (ns>=6){
      strbuf.str(args[0]); strbuf >> fx; strbuf.clear();
      strbuf.str(args[1]); strbuf >> fy; strbuf.clear();
      strbuf.str(args[2]); strbuf >> fz; strbuf.clear();
      prop_readin.frc[iat][0] = fx;
      prop_readin.frc[iat][1] = fy;
      prop_readin.frc[iat][2] = fz;
      strbuf.str(args[3]); strbuf >> td1; strbuf.clear();
      strbuf.str(args[4]); strbuf >> td2; strbuf.clear();
      strbuf.str(args[5]); strbuf >> td3; strbuf.clear();
      if (use_u.frc){
	prop_u.frc[iat][0] = td1;
	prop_u.frc[iat][1] = td2;
	prop_u.frc[iat][2] = td3;
      }
      else {
	prop_w.frc[iat][0] = td1;
	prop_w.frc[iat][1] = td2;
	prop_w.frc[iat][2] = td3;
      }
      iat++;
    }

    iline++;


    if (!fp) break;
  }
  fp.close();
  fp.clear();


  if (iat != nb)
    aborterror("ERROR: Read " + tostring(iat) + " force vectors when expecting "
	       + tostring(nb) + ". Exiting.");


}




#include "compoundfit-list.cppinc"


