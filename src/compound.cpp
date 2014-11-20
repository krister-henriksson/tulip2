



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

#include "elem-iacs.hpp"
#include "specs-fit-prop-pot.hpp"
#include "compound.hpp"
#include "physconst.hpp"
#include "exiterrors.hpp"

using namespace std;
using namespace utils;
using namespace exiterrors;
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
  Emix = false;
  B = false;
  Bp = false;
  Fmax = false;
  Pmax = false;
  displmax = false;
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
  Emix = 0;
  B = 0;
  Bp = 0;
  Fmax = 0;
  Pmax = 0;
  displmax = 0;
}

CompoundPropertiesUseUncertainties::CompoundPropertiesUseUncertainties()
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
  Emix = true;
  B = true;
  Bp = true;
  Fmax = true;
  Pmax = true;
  displmax = true;
}

CompoundPropertiesUseWeights::CompoundPropertiesUseWeights()
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
  Emix = false;
  B = false;
  Bp = false;
  Fmax = false;
  Pmax = false;
  displmax = false;
}

CompoundPropertiesUncertainties::CompoundPropertiesUncertainties()
  : C(6,6,1)
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
  Emix = 0.1;
  B = 0.1;
  Bp = 0.1;
  Fmax = 0.1;
  Pmax = 0.1;
  displmax = 0.1;
}

CompoundPropertiesWeights::CompoundPropertiesWeights()
  : C(6,6,0)
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
  Emix = 1;
  B = 1;
  Bp = 1;
  Fmax = 1;
  Pmax = 1;
  displmax = 1;
}




CompoundStructure::CompoundStructure()
  :
  elemnames(1, "none"),   // Single-atom basis, i.e. Bravais lattice
  pbc(3, false),
  origin(3,0),
  u1_vec(3,0), u2_vec(3,0), u3_vec(3,0),
  basis_elems(1, "none"),
  basis_vecs(1, Vector<double>(3, 0.0))
{

  filename = "none";
  name = "none";
  crystalname = "none";
  nelem = elemnames.size();
  csystem = "cubic";

  scalefactor = 1;

  use_int = false;

  // Orthonormal cubic system:
  u1_vec[0] = 1;
  u2_vec[1] = 1;
  u3_vec[2] = 1;

  nbasis = basis_elems.size();

  // Crystal is:
  //
  // r(n1,n2,n3,j) = n1 * u1_vec + n2 * u2_vec + n3 * u3_vec + basis_vecs[j]
  //
  // If pbc==false then n1=n2=n3=0 always, and 0<=j<nbasis.
  // Else n1,n2,n3 can vary, and, and 0<=j<nbasis.
}



// ###########################################################################
// ###########################################################################
//
// Create from model
//
// ###########################################################################
// ###########################################################################



void CompoundStructure::create_from_model(string name_in,
					  string elem1,
					  string elem2
					  ){
  
  int i;

  crystalname = name_in;
  name = crystalname + "-" + elem1;
  elemnames[0] = elem1;

  // finalize() will use 'ai' for the overall scaling of basis vectors:
  scalefactor = -1;


  pbc = Vector<bool>(3, true);

  basis_elems[0] = elem1;
  basis_vecs[0] = Vector<double>(3, 0.0);

  //use_u=false;
  //use_w=false;


  if (crystalname=="DIM1"){
    csystem = "none";
    pbc = Vector<bool>(3, false);

    nbasis = 2;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<nbasis; ++i) basis_elems[i] = elem1;
    i=0; basis_vecs[i] = Vector<double>(3, 0.0);
    i=1; basis_vecs[i] = Vector<double>(3, 0.0); basis_vecs[i][0] = 1.0; // bond axis in X direction

    return;
  }
  else if (crystalname=="DIM2"){
    nelem = 2;
    elemnames.resize(2);
    elemnames[1] = elem2;

    csystem = "none";
    pbc = Vector<bool>(3, false);

    nbasis = 2;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    i=0; basis_elems[i] = elem1;
    i=1; basis_elems[i] = elem2;

    i=0; basis_vecs[i] = Vector<double>(3, 0.0);
    i=1; basis_vecs[i] = Vector<double>(3, 0.0); basis_vecs[i][0] = 1.0; // bond axis in X direction

    return;
  }
  else if (crystalname=="SC"){
    return;
  }
  else if (crystalname=="BCC"){
    u1_vec[0] =  1.0; u1_vec[1] =  0.0; u1_vec[2] =  0.0;
    u2_vec[0] =  0.0; u2_vec[1] =  1.0; u2_vec[2] =  0.0;
    u3_vec[0] =  0.0; u3_vec[1] =  0.0; u3_vec[2] =  1.0;

    nbasis = 2;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    i=0; basis_elems[i] = elem1;
    i=1; basis_elems[i] = elem1;

    i=0; basis_vecs[i] = Vector<double>(3, 0.0);
    i=1; basis_vecs[i] = Vector<double>(3, 0.5);
    return;
  }
  else if (crystalname=="FCC"){
    nbasis = 4;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<nbasis; ++i) basis_elems[i] = elem1;
    for (i=0; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.0);

    i=1; basis_vecs[i][0] = 0.5; basis_vecs[i][1] = 0.5; basis_vecs[i][2] = 0.0;
    i=2; basis_vecs[i][0] = 0.5; basis_vecs[i][1] = 0.0; basis_vecs[i][2] = 0.5;
    i=3; basis_vecs[i][0] = 0.0; basis_vecs[i][1] = 0.5; basis_vecs[i][2] = 0.5;

    return;
  }
  else if (crystalname=="DIA"){
    nbasis = 8;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<nbasis; ++i) basis_elems[i] = elem1;
    for (i=0; i<4; ++i)      basis_vecs[i] = Vector<double>(3, 0.0);
    for (i=4; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.25);

    i=1; basis_vecs[i][0] += 0.5; basis_vecs[i][1] += 0.5; basis_vecs[i][2] += 0.0;
    i=2; basis_vecs[i][0] += 0.5; basis_vecs[i][1] += 0.0; basis_vecs[i][2] += 0.5;
    i=3; basis_vecs[i][0] += 0.0; basis_vecs[i][1] += 0.5; basis_vecs[i][2] += 0.5;

    i=5; basis_vecs[i][0] += 0.5; basis_vecs[i][1] += 0.5; basis_vecs[i][2] += 0.0;
    i=6; basis_vecs[i][0] += 0.5; basis_vecs[i][1] += 0.0; basis_vecs[i][2] += 0.5;
    i=7; basis_vecs[i][0] += 0.0; basis_vecs[i][1] += 0.5; basis_vecs[i][2] += 0.5;

    return;
  }
  else if (crystalname=="HCP"){
    csystem = "hexagonal";
    use_int = true;

    u2_vec[0] = 0.5; u2_vec[1] = 0.866025403784; u2_vec[2] = 0;

    nbasis = 2;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<nbasis; ++i) basis_elems[i] = elem1;
    for (i=0; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.0);

    i=0; basis_vecs[i] = Vector<double>(3, 0.0);
    i=1; basis_vecs[i][0] = 1.0/3.0; basis_vecs[i][1] = 1.0/3.0; basis_vecs[i][2] = 1.0/3.0;

    return;
  }
  else if (crystalname=="GRP"){ // graphene
    pbc[2] = false;
    csystem = "hexagonal";

    /*
      Primitive cell:
      a1 = a*xhat
      a2 = a/2*xhat + sqrt(3)*a/2*yhat
      a3 = 0.0
     */

    // Input lattice parameter a is called ahex, and is taken as r_nn !!!
    double ap=3.0;
    double bp=sqrt(3);

    u1_vec[0] = ap;  u1_vec[1] = 0.0; u1_vec[2] = 0.0;
    u2_vec[0] = 0.0; u2_vec[1] = bp;  u2_vec[2] = 0.0;
    u3_vec[0] = 0.0; u3_vec[1] = 0.0; u3_vec[2] = 1.0;

    nbasis = 4;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<nbasis; ++i) basis_elems[i] = elem1;
    for (i=0; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.0);

    use_int = true;

    i=0; basis_vecs[i]    = Vector<double>(3, 0.0);
    i=1; basis_vecs[i][0] = 1.0/3.0; basis_vecs[i][1] = 0;       basis_vecs[i][2] = 0;
    i=2; basis_vecs[i][0] = 1.0/2.0; basis_vecs[i][1] = 1.0/2.0; basis_vecs[i][2] = 0;
    i=3; basis_vecs[i][0] = 5.0/6.0; basis_vecs[i][1] = 1.0/2.0; basis_vecs[i][2] = 0;

    return;
  }
  else if (crystalname=="SH" || crystalname=="GRA"){
    csystem = "hexagonal";

    /*
      Primitive cell:
      a1 = a*xhat
      a2 = a/2*xhat + sqrt(3)*a/2*yhat
      a3 = c*zhat
     */

    // Input lattice parameter a is called ahex, and is taken as r_nn !!!
    double ap=3.0;
    double bp=sqrt(3);

    u1_vec[0] = ap;  u1_vec[1] = 0.0; u1_vec[2] = 0.0;
    u2_vec[0] = 0.0; u2_vec[1] = bp;  u2_vec[2] = 0.0;
    u3_vec[0] = 0.0; u3_vec[1] = 0.0; u3_vec[2] = 1.0;

    nbasis = 8;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<nbasis; ++i) basis_elems[i] = elem1;
    for (i=0; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.0);

    use_int = true;

    i=0; basis_vecs[i]    = Vector<double>(3, 0.0);
    i=1; basis_vecs[i][0] = 1.0/3.0; basis_vecs[i][1] = 0;       basis_vecs[i][2] = 0;
    i=2; basis_vecs[i][0] = 1.0/2.0; basis_vecs[i][1] = 1.0/2.0; basis_vecs[i][2] = 0;
    i=3; basis_vecs[i][0] = 5.0/6.0; basis_vecs[i][1] = 1.0/2.0; basis_vecs[i][2] = 0;

    i=4; basis_vecs[i][0] = 1.0/6.0; basis_vecs[i][1] = 1.0/2.0; basis_vecs[i][2] = 1.0/2.0;
    i=5; basis_vecs[i][0] = 1.0/3.0; basis_vecs[i][1] = 0;       basis_vecs[i][2] = 1.0/2.0;
    i=6; basis_vecs[i][0] = 2.0/3.0; basis_vecs[i][1] = 0;       basis_vecs[i][2] = 1.0/2.0;
    i=7; basis_vecs[i][0] = 5.0/6.0; basis_vecs[i][1] = 1.0/2.0; basis_vecs[i][2] = 1.0/2.0;

    return;
  }
  else if (crystalname=="CsCl"){
    name = name + "-" + elem2;
    nelem = 2;
    elemnames.resize(nelem);
    elemnames[1] = elem2;
 
    nbasis = 2;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    basis_elems[0] = elem1;
    basis_elems[1] = elem2;
    for (i=0; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.0);

    i=0; basis_vecs[i] = Vector<double>(3, 0.0);
    i=1; basis_vecs[i][0] = 0.5; basis_vecs[i][1] = 0.5; basis_vecs[i][2] = 0.5;

    return;
  }
  else if (crystalname=="NaCl"){
    name = name + "-" + elem2;
    nelem = 2;
    elemnames.resize(nelem);
    elemnames[1] = elem2;

    nbasis = 8;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<4; ++i)      basis_elems[i] = elem1;
    for (i=5; i<nbasis; ++i) basis_elems[i] = elem2;
    for (i=0; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.0);

    i=0; basis_vecs[i] = Vector<double>(3, 0.0);
    i=1; basis_vecs[i][0] = 0.5; basis_vecs[i][1] = 0.5; basis_vecs[i][2] = 0.0;
    i=2; basis_vecs[i][0] = 0.5; basis_vecs[i][1] = 0.0; basis_vecs[i][2] = 0.5;
    i=3; basis_vecs[i][0] = 0.0; basis_vecs[i][1] = 0.5; basis_vecs[i][2] = 0.5;

    i=4; basis_vecs[i][0] = 0.5; basis_vecs[i][1] = 0.0; basis_vecs[i][2] = 0.0;
    i=5; basis_vecs[i][0] = 0.0; basis_vecs[i][1] = 0.5; basis_vecs[i][2] = 0.0;
    i=6; basis_vecs[i][0] = 0.0; basis_vecs[i][1] = 0.0; basis_vecs[i][2] = 0.5;
    i=7; basis_vecs[i][0] = 0.5; basis_vecs[i][1] = 0.5; basis_vecs[i][2] = 0.5;

    return;
  }
  else if (crystalname=="ZnS"){
    name = name + "-" + elem2;
    nelem = 2;
    elemnames.resize(nelem);
    elemnames[1] = elem2;

    nbasis = 8;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<4; ++i)      basis_elems[i] = elem1;
    for (i=5; i<nbasis; ++i) basis_elems[i] = elem2;

    for (i=0; i<4; ++i)      basis_vecs[i] = Vector<double>(3, 0.0);
    for (i=4; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.25);

    i=1; basis_vecs[i][0] += 0.5; basis_vecs[i][1] += 0.5; basis_vecs[i][2] += 0.0;
    i=2; basis_vecs[i][0] += 0.5; basis_vecs[i][1] += 0.0; basis_vecs[i][2] += 0.5;
    i=3; basis_vecs[i][0] += 0.0; basis_vecs[i][1] += 0.5; basis_vecs[i][2] += 0.5;

    i=5; basis_vecs[i][0] += 0.5; basis_vecs[i][1] += 0.5; basis_vecs[i][2] += 0.0;
    i=6; basis_vecs[i][0] += 0.5; basis_vecs[i][1] += 0.0; basis_vecs[i][2] += 0.5;
    i=7; basis_vecs[i][0] += 0.0; basis_vecs[i][1] += 0.5; basis_vecs[i][2] += 0.5;

    return;
  }
  else {
    aborterror("Error: Compound type " + name + " not recognized. Exiting.");
    return;
  }

}





// ###########################################################################
// ###########################################################################
//
// Create from readin-file
//
// ###########################################################################
// ###########################################################################


void CompoundStructure::read_structure(void){
  ifstream fp;
  string line;
  vector<string> args;
  istringstream strbuf;
  int ns, i, j, k, tl;
  double td;


  /* .......................................................................
     Get basic structure data.
     ....................................................................... */
  
  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find/open file " +
	       filename + " for compound " +
	       name + ". Exiting.");
  

  int iline = 1;
  while (true){

    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    cout << line << endl;



    // iline = 1
    /* 1. Comment line: */

    // iline = 2
    /* 2. Overall scaling factor. */
    if (iline==2)
      strbuf.str(args[0]); strbuf >> scalefactor; strbuf.clear();

    // iline = 3...5
    /* 3-5.Primitive vectors: */
    if (iline==3 || iline==4 || iline==5){
      k = iline-3;

      if (ns<3)
	aborterror("ERROR: Too few coordinates for primitive vectors in file " +
		   filename + " for compound " +
		   name + ". Exiting");
      
      for (j=0; j<=2; j++){
	strbuf.str(args[j]); strbuf >> td; strbuf.clear();
	if      (k==0) u1_vec[j] = td;
	else if (k==1) u2_vec[j] = td;
	else if (k==2) u3_vec[j] = td;
      }

      if (ns>=4 && args[3]=="pbc")
	pbc[k] = true;
    }

    // iline=6
    /* 6. Basis vectors are in internal or direct format ? */
    if (iline==6){

      if (ns==0)
	aborterror("ERROR: No indicator for internal/direct format of basis vectors "
		   "given for compound " + name + ". Exiting.");
  
      if      (args[0][0]=='I' || args[0][0]=='i')
	use_int = true;
      else if (args[0][0]=='D' || args[0][0]=='d' || args[0][0]=='S' || args[0][0]=='s')
	use_int = false;
      else
	aborterror("ERROR: Unknown indicator " + args[0] +
		   " for internal/direct format of basis vectors for compound " + 
		   name + ". Exiting.");
    }


    // iline=7
    /* 7. Number of basis vectors: */
    if (iline==7){

      if (ns==0)
	aborterror("ERROR: No number of basis vectors given for compound " + name + ". Exiting.");
      
      strbuf.str(args[0]); strbuf >> tl; strbuf.clear();
      nbasis = tl;

      /* Allocate space for quantities depending on the number of basis vectors: */
      basis_elems.resize(tl);
      basis_vecs.resize(tl);
      for (i=0; i<nbasis; ++i)
	basis_vecs[i].resize(3);
    }

    // iline=8
    /* 8-(8+Nbasvec). Basis vectors: */
    if (iline >= 8 && iline < 8+nbasis){
      i = iline - 8;

      if (ns==0)
	aborterror("ERROR: Compound " + name + ": File ended prematurely when scanning "
		   "for basis vector " + tostring(i) + ". Exiting.");

      if (ns<4)
	aborterror("ERROR: Compound " + name + ": Too few fields for basis vector " + 
		   tostring(i) + ". Exiting.");

      //basis_elems[i] = args[0];
      strbuf.str(args[0]); strbuf >> basis_elems[i];   strbuf.clear();
      strbuf.str(args[1]); strbuf >> basis_vecs[i][0]; strbuf.clear();
      strbuf.str(args[2]); strbuf >> basis_vecs[i][1]; strbuf.clear();
      strbuf.str(args[3]); strbuf >> basis_vecs[i][2]; strbuf.clear();

      //cout << i << " " << basis_elems[i] << " " << basis_vecs[i] << endl;



      if (! fp)
	aborterror("Error while reading basis vectors for compound " + name + ". Exiting.");
    }

    iline++;

    if (!fp) break;
  }
  fp.close();
  fp.clear();

  finalize(scalefactor, scalefactor, scalefactor);

}



// ###########################################################################
// ###########################################################################
// ###########################################################################
// ###########################################################################

void CompoundStructure::finalize(const double ai,
				 const double bi,
				 const double ci){
  int i, k;
  double td;

  /* -----------------------------------------------------------------------------
     Get ready-to-use primitive and basis vectors.
     ----------------------------------------------------------------------------- */

  // ****************************************************************
  // Primitive vectors:
  // ****************************************************************
  for (k=0; k<3; ++k){
    u1_vec[k] *= ai;
    u2_vec[k] *= bi;
    u3_vec[k] *= ci;
  }


  // ****************************************************************
  // Basis vectors:
  // ****************************************************************
  if (! use_int){
    // Use non-internal format:
    td = (scalefactor < 0) ? ai : scalefactor;
    //td = ai;

    // Basis vectors are in non-internal format
    for (i=0; i<nbasis; ++i){
      for (k=0; k<3; ++k)
	basis_vecs[i][k] *= td;
    }
  }
  else {
    // Basis vectors are in internal format
    Vector< Vector<double> > bv;
    bv = basis_vecs;
    for (i=0; i<nbasis; ++i){
      for (k=0; k<3; ++k)
	basis_vecs[i][k] = 0.0
	  + bv[i][0] * u1_vec[k]
	  + bv[i][1] * u2_vec[k]
	  + bv[i][2] * u3_vec[k];
    }
  }

  /* -----------------------------------------------------------------------------
     Make all basis atoms be inside the cell, if possible
     ----------------------------------------------------------------------------- */





  Matrix<double> boxdir(3,3,0), Bravaismatrix_inv(3,3,0);
  Vector<double> tv(3,0), boxlen(3), drs(3,0), drc(3,0);

  tv=u1_vec; td=tv.normalize(); boxlen[0]=td; boxdir.col(0,tv);
  tv=u2_vec; td=tv.normalize(); boxlen[1]=td; boxdir.col(1,tv);
  tv=u3_vec; td=tv.normalize(); boxlen[2]=td; boxdir.col(2,tv);

  
  bool isCart = false;
  double eps = numeric_limits<double>::epsilon();;
  k = 0;
  // boxdir(0)
  tv[0] = boxdir.col(0)[0]; if (tv[0]<0) tv[0] *= -1;
  tv[1] = boxdir.col(0)[1]; if (tv[1]<0) tv[1] *= -1;
  tv[2] = boxdir.col(0)[2]; if (tv[2]<0) tv[2] *= -1;
  if (tv[1]<eps && tv[2]<eps) k++;
  // boxdir(1)
  tv[0] = boxdir.col(1)[0]; if (tv[0]<0) tv[0] *= -1;
  tv[1] = boxdir.col(1)[1]; if (tv[1]<0) tv[1] *= -1;
  tv[2] = boxdir.col(1)[2]; if (tv[2]<0) tv[2] *= -1;
  if (tv[0]<eps && tv[2]<eps) k++;
  // boxdir(2)
  tv[0] = boxdir.col(2)[0]; if (tv[0]<0) tv[0] *= -1;
  tv[1] = boxdir.col(2)[1]; if (tv[1]<0) tv[1] *= -1;
  tv[2] = boxdir.col(2)[2]; if (tv[2]<0) tv[2] *= -1;
  if (tv[0]<eps && tv[1]<eps) k++;
  if (k == 3) isCart = true;


  if (! isCart)
    boxdir.inverse( Bravaismatrix_inv );
  else {
    Bravaismatrix_inv.elem(0,0) = 1.0;
    Bravaismatrix_inv.elem(1,0) = 0.0;
    Bravaismatrix_inv.elem(2,0) = 0.0;

    Bravaismatrix_inv.elem(0,1) = 0.0;
    Bravaismatrix_inv.elem(1,1) = 1.0;
    Bravaismatrix_inv.elem(2,1) = 0.0;

    Bravaismatrix_inv.elem(0,2) = 0.0;
    Bravaismatrix_inv.elem(1,2) = 0.0;
    Bravaismatrix_inv.elem(2,2) = 1.0;
  }

  for (i=0; i<nbasis; ++i){
    // Get basis vector in skew system
    drc = basis_vecs[i];

    // Get distance in skew coordinate system, where periodics can be checked:
    drs[0] = drc[0];
    drs[1] = drc[1];
    drs[2] = drc[2];
    if (! isCart){
      drs[0] = Bravaismatrix_inv.elem(0,0) * drc[0]
	+ Bravaismatrix_inv.elem(0,1) * drc[1]
	+ Bravaismatrix_inv.elem(0,2) * drc[2];
      drs[1] = Bravaismatrix_inv.elem(1,0) * drc[0]
	+ Bravaismatrix_inv.elem(1,1) * drc[1]
	+ Bravaismatrix_inv.elem(1,2) * drc[2];
      drs[2] = Bravaismatrix_inv.elem(2,0) * drc[0]
	+ Bravaismatrix_inv.elem(2,1) * drc[1]
	+ Bravaismatrix_inv.elem(2,2) * drc[2];
    }

    // Periodics check:
    double pf1=0.0, pf2=1.0;
    while (pbc[0] && drs[0] >= pf2 * boxlen[0]) drs[0] -= boxlen[0];
    while (pbc[0] && drs[0] <  pf1 * boxlen[0]) drs[0] += boxlen[0];
    while (pbc[1] && drs[1] >= pf2 * boxlen[1]) drs[1] -= boxlen[1];
    while (pbc[1] && drs[1] <  pf1 * boxlen[1]) drs[1] += boxlen[1];
    while (pbc[2] && drs[2] >= pf2 * boxlen[2]) drs[2] -= boxlen[2];
    while (pbc[2] && drs[2] <  pf1 * boxlen[2]) drs[2] += boxlen[2];
    
    // Get distance in Cartesian coordinate system:
    drc[0] = drs[0];
    drc[1] = drs[1];
    drc[2] = drs[2];
    if (! isCart){
      drc[0] = boxdir.elem(0,0) * drs[0] + boxdir.elem(0,1) * drs[1] + boxdir.elem(0,2) * drs[2];
      drc[1] = boxdir.elem(1,0) * drs[0] + boxdir.elem(1,1) * drs[1] + boxdir.elem(1,2) * drs[2];
      drc[2] = boxdir.elem(2,0) * drs[0] + boxdir.elem(2,1) * drs[1] + boxdir.elem(2,2) * drs[2];
    }

    basis_vecs[i] = drc;
  }




  return;
}

// ###########################################################################
// ###########################################################################
// ###########################################################################
// ###########################################################################


#include "compound-strfit.cppinc"


// #####################################################################
// #####################################################################

#include "compound-listfit.cppinc"

