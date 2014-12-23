

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






void CompoundStructureFit::check_and_fix_Cij(){
  // Warn user and exit? Or rewrite to use "standard" settings?


  csystem_sub=0;

  if (csystem=="monoclinic" || csystem=="tetragonal" || csystem=="trigonal"){
    int k=0, p=0;

    if (csystem=="monoclinic"){
      if (prop_use.C.elem(1,5)) k++;
      if (prop_use.C.elem(2,5)) k++;
      if (prop_use.C.elem(3,5)) k++;
      if (prop_use.C.elem(4,6)) k++;
      
      if (prop_use.C.elem(1,6)) p++;
      if (prop_use.C.elem(2,6)) p++;
      if (prop_use.C.elem(3,6)) p++;
      if (prop_use.C.elem(4,5)) p++;

      if (k>0) csystem_sub=0;
      if (p>0) csystem_sub=1;
      if (k>0 && p>0){
	aborterror("ERROR: Monoclinic crystal system uses either C15,C25,C35,C46 or "
		   "C16,C26,C36,C45, not Cij from both. Exiting.");
      }
    }
    else if (csystem=="tetragonal"){
      if (prop_use.C.elem(1,6)) csystem_sub=0;
      else csystem_sub=1;
    }
    else if (csystem=="trigonal"){
      if (prop_use.C.elem(2,5)) csystem_sub=0;
      else csystem_sub=1;
    }
  }



  Matrix<bool> Cuse(7,7,false);
  get_Cuse(Cuse);


  for (int k=1; k<=6; ++k){
    for (int p=1; p<=6; ++p){

      if (prop_use.C.elem(k-1,p-1)==true && Cuse.elem(k,p)==false){
	string mess = "ERROR: Crystal system " + csystem + " does not use elastic constant C"
	  + tostring(k) + tostring(p) + ". Used Cij are:";
	for (int ik=1; ik<=6; ++ik)
	  for (int ip=1; ip<=6; ++ip)
	    if (Cuse.elem(ik,ip)) mess += " C" + tostring(ik) + tostring(ip);
	aborterror(mess);
      }

    }
  }



}



void CompoundStructureFit::get_Cuse(Matrix<bool> & Cuse){
  /* "Cuse" indicates which Cij elements are used by this code.
     This method compares read-in Cij with these to find
     (1) which alternate C-tensors are used, and
     (2) if user provided correct Cij elements.
   */




  if (csystem=="cubic"){
    Cuse.elem(1,1) = true;
    Cuse.elem(1,2) = true;
    Cuse.elem(4,4) = true;
  }
  else if (csystem=="hexagonal"){
    Cuse.elem(1,1) = true;
    Cuse.elem(1,2) = true;
    Cuse.elem(1,3) = true;
    Cuse.elem(3,3) = true;
    Cuse.elem(4,4) = true;
  }
  else if (csystem=="orthorombic"){
    Cuse.elem(1,1) = true;
    Cuse.elem(1,2) = true;
    Cuse.elem(1,3) = true;
    Cuse.elem(2,2) = true;
    Cuse.elem(2,3) = true;
    Cuse.elem(3,3) = true;
    Cuse.elem(4,4) = true;
    Cuse.elem(5,5) = true;
    Cuse.elem(6,6) = true;
  }
  else if (csystem=="monoclinic"){
    Cuse.elem(1,1) = true;
    Cuse.elem(1,2) = true;
    Cuse.elem(1,3) = true;
    Cuse.elem(2,2) = true;
    Cuse.elem(2,3) = true;
    Cuse.elem(3,3) = true;
    Cuse.elem(4,4) = true;
    Cuse.elem(5,5) = true;
    Cuse.elem(6,6) = true;

    if (csystem_sub==0){
      Cuse.elem(1,5) = true;
      Cuse.elem(2,5) = true;
      Cuse.elem(3,5) = true;
      Cuse.elem(4,6) = true;
    }
    else {
      Cuse.elem(1,6) = true;
      Cuse.elem(2,6) = true;
      Cuse.elem(3,6) = true;
      Cuse.elem(4,5) = true;
    }
  }
  else if (csystem=="triclinic"){
    for (int k=1; k<=6; ++k)
      for (int p=1; p<=6; ++p)
	Cuse.elem(k,p)=true;
  }
  else if (csystem=="tetragonal"){
    Cuse.elem(1,1) = true;
    Cuse.elem(1,2) = true;
    Cuse.elem(1,3) = true;
    Cuse.elem(3,3) = true;
    Cuse.elem(4,4) = true;
    Cuse.elem(6,6) = true;

    if (csystem_sub==0){
      Cuse.elem(1,6) = true;
    }
  }
  else if (csystem=="trigonal"){
    Cuse.elem(1,1) = true;
    Cuse.elem(1,2) = true;
    Cuse.elem(1,3) = true;
    Cuse.elem(3,3) = true;
    Cuse.elem(4,4) = true;
    Cuse.elem(1,4) = true;

    if (csystem_sub==0){
      Cuse.elem(2,5) = true;
    }
  }

}





#include "compoundfit-list.cppinc"


