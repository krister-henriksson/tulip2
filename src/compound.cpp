



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

    // Lattice parameter a is taken as r_nn !!!

    u1_vec[0] = 3.0; u1_vec[1] = 0.0;       u1_vec[2] = 0.0;
    u2_vec[0] = 0.0; u2_vec[1] = sqrt(3.0); u2_vec[2] = 0.0;
    u3_vec[0] = 0.0; u3_vec[1] = 0.0;       u3_vec[2] = 1.0;

    nbasis = 4;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<nbasis; ++i) basis_elems[i] = elem1;
    for (i=0; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.0);

    use_int = true;

    i=0; basis_vecs[i] = Vector<double>(3, 0.0);
    i=1; basis_vecs[i][0] = 1.0/3.0; basis_vecs[i][1] = 0;       basis_vecs[i][2] = 0;
    i=2; basis_vecs[i][0] = 1.0/2.0; basis_vecs[i][1] = 1.0/2.0; basis_vecs[i][2] = 0;
    i=3; basis_vecs[i][0] = 5.0/6.0; basis_vecs[i][1] = 1.0/2.0; basis_vecs[i][2] = 0;

    return;
  }
  else if (crystalname=="GRA"){ //
    csystem = "hexagonal";

    // Lattice parameter a is taken as r_nn !!!

    u1_vec[0] = 3.0; u1_vec[1] = 0.0;       u1_vec[2] = 0.0;
    u2_vec[0] = 0.0; u2_vec[1] = sqrt(3.0); u2_vec[2] = 0.0;
    u3_vec[0] = 0.0; u3_vec[1] = 0.0;       u3_vec[2] = 1.0;

    nbasis = 8;
    basis_elems.resize(nbasis);
    basis_vecs.resize(nbasis);

    for (i=0; i<nbasis; ++i) basis_elems[i] = elem1;
    for (i=0; i<nbasis; ++i) basis_vecs[i] = Vector<double>(3, 0.0);

    use_int = true;

    i=0; basis_vecs[i] = Vector<double>(3, 0.0);
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


// #####################################################################
// #####################################################################





void CompoundStructure::finalize(const double ai,
				 const double bi,
				 const double ci){
  Vector< Vector<double> > bv;
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
    td = (scalefactor < 0) ? ai : scalefactor;
    td = ai;

    // Basis vectors are in non-internal format
    for (i=0; i<nbasis; ++i){
      for (k=0; k<3; ++k)
	basis_vecs[i][k] *= td;
    }
  }
  else {
    // Basis vectors are in internal format
    bv = basis_vecs;
    for (i=0; i<nbasis; ++i){
      basis_vecs[i][0] = 0;
      basis_vecs[i][1] = 0;
      basis_vecs[i][2] = 0;

      for (k=0; k<3; ++k)
	basis_vecs[i][k] += 
	  bv[i][0] * u1_vec[k]
	  + bv[i][1] * u2_vec[k]
	  + bv[i][2] * u3_vec[k];
    }
  }

  /*
  for (i=0; i<nbasis; ++i){
    cout << i << " " << basis_elems[i] << " " << basis_vecs[i] << endl;
  }
  */
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



// #####################################################################
// #####################################################################
// #####################################################################
// #####################################################################


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

  return N;
}



// #####################################################################
// #####################################################################



CompoundListFit::CompoundListFit(const Elements & el,
				 MDSettings & mds_specs_general,
				 string filename)
  : elem(el), ncompounds(0)
{

  ifstream fp;
  ofstream fpo;
  string line;
  vector<string> args, opts;
  istringstream strbuf;
  int ns, nlats, ilat, nelem;
  double td;
  Vector< Vector<double> > bv;
  int ot, it;
  bool reading_latinfo, elem_ok;
  int i, j, k, p;

  double eps = std::numeric_limits<double>::epsilon();



  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");


  nlats = 0;
  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;


    if (args[0]=="LAT") nlats++;


    if (! fp) break;
  }
  fp.close();
  fp.clear();

  cout << "Detected " << nlats << " compounds." << endl;

  compounds.resize(nlats);
  ncompounds = compounds.size();


  opts.resize(0);
  opts.push_back("a");
  opts.push_back("b");
  opts.push_back("c");
  opts.push_back("bpa");
  opts.push_back("cpa");
  opts.push_back("r0");
  opts.push_back("angle_ab");
  opts.push_back("angle_ac");
  opts.push_back("angle_bc");
  opts.push_back("Vatom");
  opts.push_back("Vat");
  opts.push_back("V0");
  opts.push_back("Ecoh");
  opts.push_back("Ec");
  opts.push_back("E0");
  opts.push_back("Eatom");
  opts.push_back("Eat");
  opts.push_back("Emix");
  opts.push_back("B");
  opts.push_back("B");
  opts.push_back("Bp");
  opts.push_back("dB/dP");
  opts.push_back("C");
  opts.push_back("Fmax");
  opts.push_back("F_max");
  opts.push_back("Pmax");
  opts.push_back("P_max");
  opts.push_back("displmax");
  opts.push_back("displ_max");



  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  ilat = 0;
  reading_latinfo = false;

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;


    if (args[0]=="LAT"){
      if (reading_latinfo){
	ilat++;
      }
      else {
	reading_latinfo = true;
	k = compounds.size();
	if (k<=ilat)
	  compounds.resize(k+1);
      }


      // Initialization:
      compounds[ilat].mds_specs = mds_specs_general;

      // Since we are reading this compound from a file, it is not a reference compound:
      compounds[ilat].mds_specs.is_ref_comp = false;


    }


    // *******************************************************************
    // Compound properties, needed for e.g. MD cell construction in
    // conjunction with the information in the compound's LAT file:
    // *******************************************************************

    else if (args[0]=="name"){
      compounds[ilat].name = args[1];
    }
    else if (args[0]=="crystalname"){
      compounds[ilat].crystalname = args[1];
    }
    else if (args[0]=="file"){
      compounds[ilat].filename = args[1];
    }

    else if (args[0]=="elements"){
      compounds[ilat].elemnames.resize(0);
      for (i=1; i<ns; ++i){
	if (args[i][0]=='#') break;
	compounds[ilat].elemnames.push_back(args[i]);
      }
      nelem = compounds[ilat].elemnames.size();
      compounds[ilat].nelem = nelem;


      // Check each supplied element, if it is known from before:
      for (i=0; i<nelem; ++i){
	elem_ok=false;
	for (j=0; j<elem.nelem(); ++j){
	  if (compounds[ilat].elemnames[i] == elem.idx2name(j)){
	    // found the element
	    elem_ok=true;
	    break;
	  }
	}
	if (! elem_ok){
	  aborterror("Element " + compounds[ilat].elemnames[i] + " was not found "
		     + "in the list of known element names. Exiting.");
	}
      }
    }

    else if (args[0]=="csystem"){
      if      ( args[1][0] == 'c' || args[1][0] == 'C' )
	compounds[ilat].csystem = "cubic";
      else if ( args[1][0] == 'h' || args[1][0] == 'H' )
	compounds[ilat].csystem = "hexagonal";
      else if ( args[1][0] == 'o' || args[1][0] == 'O' )
	compounds[ilat].csystem = "orthorombic";
      else {
	aborterror("ERROR: Unknown crystal system " + args[1] + ". Exiting.");
      }
    }



    // *******************************************************************
    // Compound-specific conditions for e.g. MD relaxation:
    // *******************************************************************

    else if (args[0]=="option"){

      if      (args[1]=="no_heating"){
	compounds[ilat].mds_specs.heating_allowed = false;
      }
      else if (args[1]=="fixed_geometry"){
	compounds[ilat].mds_specs.fixed_geometry = true;
      }
      else if (args[1]=="quench_always"){
	compounds[ilat].mds_specs.quench_always = true;
      }

    }

    else if (args[0][0]=='m' && args[0][1]=='d' && args[0][2]=='s'){

      if      (args[0]=="mds_skint"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.skint;
	strbuf.clear();
      }
      else if (args[0]=="mds_seed"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.seed;
	strbuf.clear();
      }
      else if (args[0]=="mds_ndump"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.ndump;
	strbuf.clear();
      }
      else if (args[0]=="mds_tstart"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.tstart;
	strbuf.clear();
      }
      else if (args[0]=="mds_tend"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.tend;
	strbuf.clear();
      }
      else if (args[0]=="mds_Tstart"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.Tstart;
	strbuf.clear();
      }
      else if (args[0]=="mds_dt"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.dt;
	strbuf.clear();
      }
      else if (args[0]=="mds_max_dt"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.max_dt;
	strbuf.clear();
      }
      else if (args[0]=="mds_max_dE"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.max_dE;
	strbuf.clear();
      }
      else if (args[0]=="mds_max_dr"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.max_dr;
	strbuf.clear();
      }
      else if (args[0]=="mds_btc_tau"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.btc_tau;
	if (compounds[ilat].mds_specs.btc_tau<0.0 ||
	    abs(compounds[ilat].mds_specs.btc_tau)<eps)
	  compounds[ilat].mds_specs.use_Tcontrol = false;
	else
	  compounds[ilat].mds_specs.use_Tcontrol = true;
	strbuf.clear();
      }

      else if (args[0]=="mds_btc_T0"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.btc_T0;
	strbuf.clear();
      }
      else if (args[0]=="mds_bpc_tau"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.bpc_tau;

	if (compounds[ilat].mds_specs.bpc_tau<0.0 ||
	    abs(compounds[ilat].mds_specs.bpc_tau)<eps)
	  compounds[ilat].mds_specs.use_Pcontrol = false;
	else
	  compounds[ilat].mds_specs.use_Pcontrol = true;
	strbuf.clear();

      }
      else if (args[0]=="mds_bpc_P0"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.bpc_P0;
	strbuf.clear();
      }
      else if (args[0]=="mds_bpc_scale"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.bpc_scale;
	strbuf.clear();
      }
      else if (args[0]=="mds_quench_tstart"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.quench_tstart;
	strbuf.clear();
      }
      else if (args[0]=="mds_quench_rate"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.quench_rate;
	strbuf.clear();
      }

      // OTHER OPTIONS POSSIBLE TO ADD HERE


    }


    // *******************************************************************
    // Compound properties to use as targets when fitting:
    // *******************************************************************

    else {
      // args[0] is being investigated.
      // It could be e.g. 'a', 'w_a', or 'u_a'


      ot = it = -1;
      for (unsigned int i=0; i<opts.size(); ++i){
	if      (args[0]==opts[i]){
	  it = i; ot = 1; break;
	}
	else if (args[0]==("w_" + opts[i])){
	  it = i; ot = 2; break;
	}
	else if (args[0]==("u_" + opts[i])){
	  it = i; ot = 3; break;
	}
      }

      string match="none";
      if (it>=0){
	match = opts[it];
	strbuf.str(args[1]); strbuf >> td; strbuf.clear();
      }
      
      if (match=="a"){
	if (ot==1){
	  compounds[ilat].prop_readin.a = td;
	  compounds[ilat].prop_use.a    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.a = td; compounds[ilat].use_w.a = true;  compounds[ilat].use_u.a = false; }
	else if (ot==3) { compounds[ilat].prop_u.a = td; compounds[ilat].use_w.a = false; compounds[ilat].use_u.a = true; }
      }
      else if (match=="b"){
	if (ot==1){
	  compounds[ilat].prop_readin.b = td;
	  compounds[ilat].prop_use.b    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.b = td; compounds[ilat].use_w.b = true;  compounds[ilat].use_u.b = false; }
	else if (ot==3) { compounds[ilat].prop_u.b = td; compounds[ilat].use_w.b = false; compounds[ilat].use_u.b = true; }
      }
      else if (match=="c"){
	if (ot==1){
	  compounds[ilat].prop_readin.c = td;
	  compounds[ilat].prop_use.c    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.c = td; compounds[ilat].use_w.c = true;  compounds[ilat].use_u.c = false; }
	else if (ot==3) { compounds[ilat].prop_u.c = td; compounds[ilat].use_w.c = false; compounds[ilat].use_u.c = true; }
      }
      else if (match=="bpa"){
	if (ot==1){
	  compounds[ilat].prop_readin.bpa = td;
	  compounds[ilat].prop_use.bpa    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.bpa = td; compounds[ilat].use_w.bpa = true;  compounds[ilat].use_u.bpa = false; }
	else if (ot==3) { compounds[ilat].prop_u.bpa = td; compounds[ilat].use_w.bpa = false; compounds[ilat].use_u.bpa = true; }
      }
      else if (match=="cpa"){
	if (ot==1){
	  compounds[ilat].prop_readin.cpa = td;
	  compounds[ilat].prop_use.cpa    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.cpa = td; compounds[ilat].use_w.cpa = true;  compounds[ilat].use_u.cpa = false; }
	else if (ot==3) { compounds[ilat].prop_u.cpa = td; compounds[ilat].use_w.cpa = false; compounds[ilat].use_u.cpa = true; }
      }
      else if (match=="r0"){
	if (ot==1){
	  compounds[ilat].prop_readin.r0 = td;
	  compounds[ilat].prop_use.r0    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.r0 = td; compounds[ilat].use_w.r0 = true;  compounds[ilat].use_u.r0 = false; }
	else if (ot==3) { compounds[ilat].prop_u.r0 = td; compounds[ilat].use_w.r0 = false; compounds[ilat].use_u.r0 = true; }
      }
      else if (match=="angle_ab"){
	if (ot==1){
	  compounds[ilat].prop_readin.angle_ab = td;
	  compounds[ilat].prop_use.angle_ab    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.angle_ab = td; compounds[ilat].use_w.angle_ab = true;  compounds[ilat].use_u.angle_ab = false; }
	else if (ot==3) { compounds[ilat].prop_u.angle_ab = td; compounds[ilat].use_w.angle_ab = false; compounds[ilat].use_u.angle_ab = true; }
      }
      else if (match=="angle_ac"){
	if (ot==1){
	  compounds[ilat].prop_readin.angle_ac = td;
	  compounds[ilat].prop_use.angle_ac    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.angle_ac = td; compounds[ilat].use_w.angle_ac = true;  compounds[ilat].use_u.angle_ac = false; }
	else if (ot==3) { compounds[ilat].prop_u.angle_ac = td; compounds[ilat].use_w.angle_ac = false; compounds[ilat].use_u.angle_ac = true; }
      }
      else if (match=="angle_bc"){
	if (ot==1){
	  compounds[ilat].prop_readin.angle_bc = td;
	  compounds[ilat].prop_use.angle_bc    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.angle_bc = td; compounds[ilat].use_w.angle_bc = true;  compounds[ilat].use_u.angle_bc = false; }
	else if (ot==3) { compounds[ilat].prop_u.angle_bc = td; compounds[ilat].use_w.angle_bc = false; compounds[ilat].use_u.angle_bc = true; }
      }
      else if (match=="V0" || match=="Vat" || match=="Vatom"){
	if (ot==1){
	  compounds[ilat].prop_readin.Vatom = td;
	  compounds[ilat].prop_use.Vatom    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Vatom = td; compounds[ilat].use_w.Vatom = true;  compounds[ilat].use_u.Vatom = false; }
	else if (ot==3) { compounds[ilat].prop_u.Vatom = td; compounds[ilat].use_w.Vatom = false; compounds[ilat].use_u.Vatom = true; }
      }
      else if (match=="Ec" || match=="Ecoh" || match=="Eat" || match=="Eatom"){
	if (ot==1){
	  compounds[ilat].prop_readin.Ecoh = td;
	  compounds[ilat].prop_use.Ecoh    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Ecoh = td; compounds[ilat].use_w.Ecoh = true;  compounds[ilat].use_u.Ecoh = false; }
	else if (ot==3) { compounds[ilat].prop_u.Ecoh = td; compounds[ilat].use_w.Ecoh = false; compounds[ilat].use_u.Ecoh = true; }
      }
      else if (match=="Emix"){
	if (ot==1){
	  compounds[ilat].prop_readin.Emix = td;
	  compounds[ilat].prop_use.Emix    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Emix = td; compounds[ilat].use_w.Emix = true;  compounds[ilat].use_u.Emix = false; }
	else if (ot==3) { compounds[ilat].prop_u.Emix = td; compounds[ilat].use_w.Emix = false; compounds[ilat].use_u.Emix = true; }
      }
      else if (match=="B"){
	if (ot==1){
	  compounds[ilat].prop_readin.B = td;
	  compounds[ilat].prop_use.B    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.B = td; compounds[ilat].use_w.B = true;  compounds[ilat].use_u.B = false; }
	else if (ot==3) { compounds[ilat].prop_u.B = td; compounds[ilat].use_w.B = false; compounds[ilat].use_u.B = true; }
      }
      else if (match=="Bp" || match=="dB/dP"){
	if (ot==1){
	  compounds[ilat].prop_readin.Bp = td;
	  compounds[ilat].prop_use.Bp    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Bp = td; compounds[ilat].use_w.Bp = true;  compounds[ilat].use_u.Bp = false; }
	else if (ot==3) { compounds[ilat].prop_u.Bp = td; compounds[ilat].use_w.Bp = false; compounds[ilat].use_u.Bp = true; }
      }
      else if (match[0]=='C'){
	int i = match[1];
	int j = match[2];

	if (i>=1 && i<=6 && j>=1 && j<=6){
	  if (ot==1){
	    compounds[ilat].prop_readin.C.elem(i-1,j-1) = td;
	    compounds[ilat].prop_use.C.elem(i-1,j-1)    = true;
	  }
	  else if (ot==2) { compounds[ilat].prop_w.C.elem(i-1,j-1) = td; compounds[ilat].use_w.C.elem(i-1,j-1) = true;  compounds[ilat].use_u.C.elem(i-1,j-1) = false; }
	  else if (ot==3) { compounds[ilat].prop_u.C.elem(i-1,j-1) = td; compounds[ilat].use_w.C.elem(i-1,j-1) = false; compounds[ilat].use_u.C.elem(i-1,j-1) = true; }
	}
      }
      else if (match=="Fmax" || match=="F_max"){
	if (ot==1){
	  compounds[ilat].prop_readin.Fmax = td;
	  compounds[ilat].prop_use.Fmax    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Fmax = td; compounds[ilat].use_w.Fmax = true;  compounds[ilat].use_u.Fmax = false; }
	else if (ot==3) { compounds[ilat].prop_u.Fmax = td; compounds[ilat].use_w.Fmax = false; compounds[ilat].use_u.Fmax = true; }
      }
      else if (match=="Pmax" || match=="P_max"){
	if (ot==1){
	  compounds[ilat].prop_readin.Pmax = td;
	  compounds[ilat].prop_use.Pmax    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Pmax = td; compounds[ilat].use_w.Pmax = true;  compounds[ilat].use_u.Pmax = false; }
	else if (ot==3) { compounds[ilat].prop_u.Pmax = td; compounds[ilat].use_w.Pmax = false; compounds[ilat].use_u.Pmax = true; }
      }
      else if (match=="displmax" || match=="displ_max"){
	if (ot==1){
	  compounds[ilat].prop_readin.displmax = td;
	  compounds[ilat].prop_use.displmax    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.displmax = td; compounds[ilat].use_w.displmax = true;  compounds[ilat].use_u.displmax = false; }
	else if (ot==3) { compounds[ilat].prop_u.displmax = td; compounds[ilat].use_w.displmax = false; compounds[ilat].use_u.displmax = true; }
      }

    }



    if (!fp) break; 
  }
  fp.clear();
  fp.close();

  cout << "Compounds information collected." << endl;

  /*
  if (compounds.size() != ilat)
    compounds.resize(ilat+1);
  */
  // ncompounds = compounds.size();



  cout << "Number of compounds read: " << compounds.size() << endl;
  cout << "Debugging compounds ..." << endl;
  /* -----------------------------------------------------------------------------
     Debugging of settings:
     ----------------------------------------------------------------------------- */
  for (ilat=0; ilat<compounds.size(); ++ilat){
    if (compounds[ilat].mds_specs.fixed_geometry){
      compounds[ilat].prop_use.displmax   = false;
      compounds[ilat].prop_use.a = false;
      compounds[ilat].prop_use.b = false;
      compounds[ilat].prop_use.c = false;
      compounds[ilat].prop_use.bpa = false;
      compounds[ilat].prop_use.cpa = false;
      compounds[ilat].prop_use.angle_ab = false;
      compounds[ilat].prop_use.angle_ac = false;
      compounds[ilat].prop_use.angle_bc = false;
      compounds[ilat].prop_use.r0 = false;
      compounds[ilat].prop_use.Vatom = false;
    }

    // Make sure we are not using the same properties in different
    // disguises, since this will lead to problems in some fitting
    // algorithms:
    compounds[ilat].check_and_fix_uses();

  }


  /* Some notes:

     - If a,b,c is mentioned in the geometry file they will be used. No checks or
     comparisons is made with the length of the primitive vectors in the LAT files.
   */


    
    
    

  /* -----------------------------------------------------------------------------
     Now get lattice data for the read-in structures.
     ----------------------------------------------------------------------------- */
  cout << "Reading structural info for compounds ... " << endl;
  for (ilat=0; ilat<compounds.size(); ++ilat){
    compounds[ilat].read_structure();
    // calls finalize() internally to insert a,b,c
  }



  /* ##############################################################################
     ##############################################################################
     Some debugging
     ##############################################################################
     ##############################################################################    
  */


  for (ilat=0; ilat<compounds.size(); ++ilat){

    if (compounds[ilat].name == "none")
      aborterror("Error: No compound name specified for compound " + tostring(ilat) + ". Exiting.");

    if (compounds[ilat].filename == "none")
      aborterror("Error: No file name specified for compound " + tostring(ilat) + ". Exiting.");

    if (compounds[ilat].nelem == 0)
      aborterror("Error: No elements specified for compound " + tostring(ilat) + ". Exiting.");

    if (compounds[ilat].nbasis == 0)
      aborterror("Error: No basis specified for compound " + tostring(ilat) + ". Exiting.");



    if (!compounds[ilat].pbc[0] ||
	!compounds[ilat].pbc[1] ||
	!compounds[ilat].pbc[2]){
      compounds[ilat].prop_use.B  = false;
      compounds[ilat].prop_use.Bp = false;
      compounds[ilat].prop_use.Pmax = false;
      for (k=0; k<6; ++k)
	for (p=0; p<6; ++p)
	  compounds[ilat].prop_use.C.elem(k,p) = false;
    }


    for (j=0; j<compounds[ilat].nbasis; ++j){
      elem_ok=false;
      for (i=0; i<compounds[ilat].nelem; ++i){
	if (compounds[ilat].basis_elems[j] == compounds[ilat].elemnames[i]){
	  // found the element
	  elem_ok=true;
	  break;
	}
      }
      if (! elem_ok){
	aborterror("Error: Compound " + tostring(ilat) + " basis vector " + tostring(j) + " "
		   "with element " + compounds[ilat].basis_elems[j] + " was not found "
		   + "in the list of known element names for this compound. Exiting.");
      }
    }
  }


  cout << "Compound information read-in completed." << endl;

  return;
}



int CompoundListFit::NData(){
  int N=0, i;

  for (i=0; i<compounds.size(); ++i)
    N += compounds[i].NData();

  return N;
}


