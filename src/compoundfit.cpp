

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
#include "utils-math.hpp"
#include "utils-vecalg.hpp"

#include "elem-iacs.hpp"
#include "specs-fit-prop-pot.hpp"
#include "compound.hpp"
#include "physconst.hpp"

#include "compoundfit.hpp"

#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;


using namespace utils;
using namespace constants;
using boost::format;









CompoundPropertiesUse::CompoundPropertiesUse()
  : C(6,6,false)
{
  a = false;
  b = false;
  c = false;
  bpa = false;
  cpa = false;
  r0 = false;
  angle_ab = false;
  angle_ac = false;
  angle_bc = false;
  Vatom = false;
  Ecoh = false;
  Ecoh_delta = false;
  Emix = false;
  Eform = false;
  B = false;
  Bp = false;
  Fmax = false;
  Pmax = false;
  displmax = false;
  frc = false;

  bondlen   = false;
  bondangle = false;
}


void CompoundPropertiesUse::check_and_fix(){
  if (bpa){
    if (a && b) b=false;
  }
  if (cpa){
    if (a && c) c=false;
  }
  if (Vatom){
    if (a && (b || bpa) && (c || cpa)) Vatom=false;
  }
  if (r0)
    a=b=c=bpa=cpa=Vatom=angle_ab=angle_ac=angle_bc=0;
}




CompoundProperties::CompoundProperties()
  : C(6,6,0)
{
  a = -1;
  b = -1;
  c = -1;
  bpa = -1;
  cpa = -1;
  r0 = -1;
  angle_ab = 0.5*PI;
  angle_ac = 0.5*PI;
  angle_bc = 0.5*PI;
  Vatom = 0;
  Ecoh = 0;
  Ecoh_delta = 0;
  Emix = 0;
  Eform = 0;
  B = 0;
  Bp = 0;
  Fmax = 0;
  Pmax = 0;
  displmax = 0;
}

CompoundPropertiesUseUncertainties::CompoundPropertiesUseUncertainties()
  : C(6,6,false)
{
  a = false;
  b = false;
  c = false;
  bpa = false;
  cpa = false;
  r0 = false;
  angle_ab = false;
  angle_ac = false;
  angle_bc = false;
  Vatom = false;
  Ecoh = false;
  Ecoh_delta = false;
  Emix = false;
  Eform = false;
  B = false;
  Bp = false;
  Fmax = false;
  Pmax = false;
  displmax = false;
  frc = false;
}

CompoundPropertiesUseWeights::CompoundPropertiesUseWeights()
  : C(6,6,true)
{
  a = true;
  b = true;
  c = true;
  bpa = true;
  cpa = true;
  r0 = true;
  angle_ab = true;
  angle_ac = true;
  angle_bc = true;
  Vatom = true;
  Ecoh = true;
  Ecoh_delta = false;
  Emix = true;
  Eform = true;
  B = true;
  Bp = true;
  Fmax = true;
  Pmax = true;
  displmax = true;
  frc = true;
}

CompoundPropertiesUncertainties::CompoundPropertiesUncertainties()
  : C(6,6,0.1)
{
  a = 0.1;
  b = 0.1;
  c = 0.1;
  bpa = 0.1;
  cpa = 0.1;
  r0 = 0.1;
  angle_ab = 0.1;
  angle_ac = 0.1;
  angle_bc = 0.1;
  Vatom = 0.1;
  Ecoh = 0.1;
  Ecoh_delta = 0.1;
  Emix = 0.1;
  Eform = 0.1;
  B = 0.1;
  Bp = 0.1;
  Fmax = 0.1;
  Pmax = 0.1;
  displmax = 0.1;
}

CompoundPropertiesWeights::CompoundPropertiesWeights()
  : C(6,6,1)
{
  a = 1;
  b = 1;
  c = 1;
  bpa = 1;
  cpa = 1;
  r0 = 1;
  angle_ab = 1;
  angle_ac = 1;
  angle_bc = 1;
  Vatom = 1;
  Ecoh = 1;
  Ecoh_delta = 1;
  Emix = 1;
  Eform = 1;
  B = 1;
  Bp = 1;
  Fmax = 1;
  Pmax = 1;
  displmax = 1;
}





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
  if (prop_use.Ecoh_delta) N++;
  if (prop_use.Emix) N++;
  if (prop_use.Eform) N++;
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
  if (prop_use.bondlen){
    N += prop_readin.bondlen.size();
  }
  if (prop_use.bondangle){
    N += prop_readin.bondangle.size();
  }

  return N;
}





void CompoundStructureFit::read_forces(void){
  std::ifstream fp;
  std::string line;
  std::vector<std::string> args;
  std::istringstream strbuf;
  int ns, i;
  double fx, fy, fz, td1, td2, td3;
  

  prop_use.frc = true;

  // Allocate space for forces:
  int nb = basis_elems.size();

  prop_readin.frc.resize(nb);
  for (i=0; i<nb; ++i) prop_readin.frc[i] = Vector3<double>(0);

  // ***************************************************************************
  // Allocate space for predicted forces. Do this now when forces are read,
  // so we won't have to wonder later if they are allocated or not.
  // ***************************************************************************
  prop_pred.frc.resize(nb);
  for (i=0; i<nb; ++i) prop_pred.frc[i] = Vector3<double>(0);

  if (use_u.frc){
    prop_u.frc.resize(nb);
    for (i=0; i<nb; ++i)
      prop_u.frc[i] = Vector3<double>(0);
  }
  else {
    prop_w.frc.resize(nb);
    for (i=0; i<nb; ++i)
      prop_w.frc[i] = Vector3<double>(0);
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
    ns = utils::get_substrings( line, args, "\t :,()[]=" );

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





// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################



void CompoundStructureFit::check_and_fix_Cij(){
  // Warn user and exit? Or rewrite to use "standard" settings?


  if ( pbc[0] && pbc[1] && pbc[2] && csystem_sub<0){
    aborterror("ERROR: Failed to find symmetry operations for compound "
	       + name + " with stated crystal system " + csystem
	       + ". Exiting.");
  }


#if 0
  // ############################################################################
  // ############################################################################


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

      if (k>0) csystem_sub=1;
      if (p>0) csystem_sub=2;
      if (k>0 && p>0){
	aborterror("ERROR: Monoclinic crystal system uses either C15,C25,C35,C46 or "
		   "C16,C26,C36,C45, not Cij from both. Exiting.");
      }
    }
    else if (csystem=="tetragonal"){
      if (prop_use.C.elem(1,6)) csystem_sub=1;
      else csystem_sub=2;
    }
    else if (csystem=="trigonal"){
      if (prop_use.C.elem(2,5)) csystem_sub=1;
      else csystem_sub=2;
    }
  }
#endif


  // ###############################################################33
  // Find standard Cij elements to use:
  // ###############################################################33
  Matrix<bool> Cuse(7,7,false);
  get_Cuse(Cuse);


  for (int k=1; k<=6; ++k){
    for (int p=1; p<=6; ++p){

      if (prop_use.C.elem(k-1,p-1)==true && Cuse.elem(k,p)==false){
	std::string mess = "Warning: Crystal system " + csystem + " does not use elastic constant C"
	  + tostring(k) + tostring(p) + ". Used Cij are:";
	for (int ik=1; ik<=6; ++ik)
	  for (int ip=1; ip<=6; ++ip)
	    if (Cuse.elem(ik,ip)) mess += " C" + tostring(ik) + tostring(ip);

	prop_use.C.elem(k-1,p-1) = false;
	std::cout << mess
		  << " Turned off wrong usage." << std::endl;
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

    if (csymaxis=="y"){
      Cuse.elem(1,5) = true;
      Cuse.elem(2,5) = true;
      Cuse.elem(3,5) = true;
      Cuse.elem(4,6) = true;
    }
    else if (csymaxis=="z"){
      Cuse.elem(1,6) = true;
      Cuse.elem(2,6) = true;
      Cuse.elem(3,6) = true;
      Cuse.elem(4,5) = true;
    }
    else {
      aborterror("Error: Crystal symmetry axis for monoclinic compound " + name
		 + " is neither y nor z. Fix and try again.");
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
    
    if (pointgroup=="4" || pointgroup=="bar(4)" || pointgroup=="4/m"){
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

    if (pointgroup=="3" || pointgroup=="bar(3)"){
      Cuse.elem(2,5) = true;
    }
  }

}






void CompoundStructureFit::get_Cresolved(Matrix<double> & Clincomb,
					 Matrix<double> & Cfull
					 ){
  double C11, C12, C13, C14, C15, C16;
  double C22, C23, C24, C25, C26;
  double C33, C34, C35, C36;
  double C44, C45, C46;
  double C55, C56;
  double C66;

  C11 = C12 = C13 = C14 = C15 = C16 = 0.0;
  C22 = C23 = C24 = C25 = C26 = 0.0;
  C33 = C34 = C35 = C36 = 0.0;
  C44 = C45 = C46 = 0.0;
  C55 = C56 = 0.0;
  C66 = 0.0;


  // ----------------------------------------------------
  // General diagonal elements:
  // ----------------------------------------------------
  C11 = Clincomb.elem(1,1) * 2;
  C22 = Clincomb.elem(2,2) * 2;
  C33 = Clincomb.elem(3,3) * 2;
  C44 = Clincomb.elem(4,4) / 2;
  C55 = Clincomb.elem(5,5) / 2;
  C66 = Clincomb.elem(6,6) / 2;


  // Make sure used diagonal elements are correct:
  if (csystem=="cubic"){
    C33 = C22 = C11;
    C66 = C55 = C44;
  }
  else if (csystem=="hexagonal" ||
	   csystem=="tetragonal" ||
	   csystem=="trigonal"){
    C22 = C11;
    C55 = C44;
  }


  // ----------------------------------------------------
  // General non-diagonal elements:
  // ----------------------------------------------------
  C12 = (Clincomb.elem(1,2) * 2 - C11 - C22 ) / 2.0;
  C13 = (Clincomb.elem(1,3) * 2 - C11 - C33 ) / 2.0;
  C14 = (Clincomb.elem(1,4) * 2 - C11 - 4 * C44 ) / 4.0;
  C15 = (Clincomb.elem(1,5) * 2 - C11 - 4 * C55 ) / 4.0;
  C16 = (Clincomb.elem(1,6) * 2 - C11 - 4 * C66 ) / 4.0;

  C23 = (Clincomb.elem(2,3) * 2 - C22 - C33 ) / 2.0;
  C24 = (Clincomb.elem(2,4) * 2 - C22 - 4 * C44 ) / 4.0;
  C25 = (Clincomb.elem(2,5) * 2 - C22 - 4 * C55 ) / 4.0;
  C26 = (Clincomb.elem(2,6) * 2 - C22 - 4 * C66 ) / 4.0;

  C34 = (Clincomb.elem(3,4) * 2 - C33 - 4 * C44 ) / 4.0;
  C35 = (Clincomb.elem(3,5) * 2 - C33 - 4 * C55 ) / 4.0;
  C36 = (Clincomb.elem(3,6) * 2 - C33 - 4 * C66 ) / 4.0;

  C45 = (Clincomb.elem(4,5) * 2 - 4 * C44 - 4 * C55 ) / 8.0;
  C46 = (Clincomb.elem(4,6) * 2 - 4 * C44 - 4 * C66 ) / 8.0;

  C56 = (Clincomb.elem(5,6) * 2 - 4 * C55 - 4 * C66 ) / 8.0;


  // Make sure used non-diagonal elements are correct:
  if (csystem=="cubic"){
    C13 = C12;
    C23 = C12;
  }
  else if (csystem=="hexagonal" ||
	   csystem=="tetragonal" ||
	   csystem=="trigonal"){
    C23 = C13;
  }



  // ----------------------------------------------------
  // Construct the upper triangular part of the elastic constants matrix:
  // ----------------------------------------------------

  Cfull.elem(0,0) = C11;
  Cfull.elem(1,1) = C22;
  Cfull.elem(2,2) = C33;
  Cfull.elem(3,3) = C44;
  Cfull.elem(4,4) = C55;
  Cfull.elem(5,5) = C66;

  Cfull.elem(0,1) = C12;
  Cfull.elem(0,2) = C13;
  Cfull.elem(0,3) = C14;
  Cfull.elem(0,4) = C15;
  Cfull.elem(0,5) = C16;

  Cfull.elem(1,2) = C23;
  Cfull.elem(1,3) = C24;
  Cfull.elem(1,4) = C25;
  Cfull.elem(1,5) = C26;

  Cfull.elem(2,3) = C34;
  Cfull.elem(2,4) = C35;
  Cfull.elem(2,5) = C36;

  Cfull.elem(3,4) = C45;
  Cfull.elem(3,5) = C46;

  Cfull.elem(4,5) = C56;

}



