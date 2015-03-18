
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include <cmath>


#include "utils.hpp"
#include "utils-math.hpp"
#include "constants.hpp"
#include "utils-matrix.hpp"
#include "utils-string.hpp"
#include "utils-vector.hpp"
#include "utils-errors.hpp"
#include "utils-vecalg.hpp"

#include "compound.hpp"
#include "compoundfit.hpp"
#include "physconst.hpp"
#include "errors.hpp"

#include "lattice-simple.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;


#define USE_SPGLIB



#ifdef USE_SPGLIB

// Thanks to Morten for pointing out this construct for e.g. clang++:
extern "C" {
#include <spglib/spglib.h>
}

#endif


using namespace utils;
using namespace constants;
using namespace physconst;
using boost::format;


LatticeSimple::LatticeSimple(){
  nbasis   = 0;
  origin   = Vector3<double>(0.0);
  minpos   = Vector3<double>(0.0);
  avec     = Vector3<double>(0.0);
  bvec     = Vector3<double>(0.0);
  cvec     = Vector3<double>(0.0);
  pbc      = Vector3<bool>(true);
  csystem  = "unknown";
  csystem_sub = 0;
  csymaxis = "";
  pointgroup = "1";
  spacegroup = "1";
  spacegroup_number = 0;
}


void LatticeSimple::get_ipos(void){
  MatrixSq3<double> A(0), Ainv(0);
  Vector3<double> b(0), x(0);

  A.elem(0,0)=avec[0];
  A.elem(1,0)=avec[1];
  A.elem(2,0)=avec[2];
  A.elem(0,1)=bvec[0];
  A.elem(1,1)=bvec[1];
  A.elem(2,1)=bvec[2];
  A.elem(0,2)=cvec[0];
  A.elem(1,2)=cvec[1];
  A.elem(2,2)=cvec[2];

  for (int i=0; i<nbasis; ++i){
    b = pos[i];
    A.solve(b,x,Ainv);
    ipos[i] = x;
    for (int j=0; j<3; ++j){
      while (ipos[i][j]< 0.0) ipos[i][j]++;
      while (ipos[i][j]>=1.0) ipos[i][j]--;
    }
  }
}


void LatticeSimple::rotate(const MatrixSq3<double> & R){
  Vector3<double> tvec(0);
  tvec = R * avec; avec = tvec;
  tvec = R * bvec; bvec = tvec;
  tvec = R * cvec; cvec = tvec;
  for (int i=0; i<nbasis; ++i){
    tvec = R * pos[i]; pos[i] = tvec;
  }
  get_ipos();
}

void LatticeSimple::shift(const Vector3<double> & svec){
  for (int i=0; i<nbasis; ++i){
    pos[i][0] += svec[0];
    pos[i][1] += svec[1];
    pos[i][2] += svec[2];
  }
  get_ipos();
}




void LatticeSimple::dump_xyz(const std::string & filename){
  std::string dumpfile(filename);
  std::ofstream fout;
  fout.open(dumpfile.c_str());
  fout << nbasis << std::endl;
  fout << "" << std::endl;
  for (int i=0; i<nbasis; ++i){
    fout << "Au " 
	 << pos[i][0] << " "
	 << pos[i][1] << " "
	 << pos[i][2] << " " << i << " 1" << std::endl;
  }
  fout.close();
  fout.clear();
}




bool LatticeSimple::matrix_is_a_symm_op(const MatrixSq3<double> & R){
  int i,j,k,p,q, ts, tt;
  Vector3<double> bv, r, rs, rt, dr;
  bool found;
  int nfound=0;
  double tol = 1e-3;


  int imin=-1, imax=1;
  int jmin=-1, jmax=1;
  int kmin=-1, kmax=1;
  //cout << "name " << name << endl;

  for (q=0; q<nbasis; ++q){
    //cout << "q nbasis " << q << " " << nbasis << endl;

    r = origin + pos[q]; // transformed




    rs = R * r; // transformed
    //rs = R * pos[q]; // transformed
    ts = types[q];
    //cout << "R = " << R << endl;

    found = false;
    for (i=imin; i<=imax; ++i){
      for (j=jmin; j<=jmax; ++j){
	for (k=kmin; k<=kmax; ++k){
	  for (p=0; p<nbasis; ++p){
	    //cout << "p nbasis " << p << " " << nbasis << endl;
	    rt = origin + i*avec + j*bvec + k*cvec + pos[p];
	    //rtt = i*avec + j*bvec + k*cvec + pos[p];
	    tt = types[p];
	    dr = rs - rt;

	    if (tt == ts &&
		fp_is_small_tol(dr[0], tol) &&
		fp_is_small_tol(dr[1], tol) &&
		fp_is_small_tol(dr[2], tol)){
	      found=true;
	      nfound++;
	      /*
	      cout << "basis elem " << basis_elems[q] << " " << q
		   << " original basis elem " << basis_elems[p] << " " 
		   << basis_vecs[q] << " is symmetric with " << rtt
		   << "i j k p : "
		   << i << " " << j << " " << k << " " << p
		   << endl;
	      */
	      //i=1; j=1; k=1; p=nbasis;
	    }
	    if (found) break;
	  }
	  if (found) break;
	}
	if (found) break;
      }
      if (found) break;
    }
    
  }

  if (nfound==nbasis) return true;
  else return false;
}


// ########################################################################
// ########################################################################
// ########################################################################
// ########################################################################



void latsymm(Vector<CompoundStructureFit> & cmplist){


  /* Pages 280-3 of Nye:
     Laboratory system: x,y,z
     System used for tensors and matrices: X,Y,Z. These form a Cartesian system.

     Triclinic:
     - no special requirement for primitive vectors

     Monoclinic:
     - Y || y, sometimes: Z || y
     - alpha = gamma = 90 deg, y parallel to 2-fold axis

     Tetragonal, trigonal, hexagonal:
     - Z || z, X || x
     * tetragonal: a=b, alpha = beta = gamma = 90 deg, z parallel to 4-fold axes
     * trigonal: a=b, alpha = beta = 90 deg, gamma = 120 deg, z parallel to 3-fold axis
     * hexagonal: a=b, alpha = beta = 90 deg, gamma = 120 deg, z parallel to 6-fold axis

     Orthorombic and cubic:
     - X || x, Y || y, Z || z
     - alpha = beta = gamma = 90 deg
     * orthorombic: x,y,z parallel to 2-fold axes
     * cubic: x,y,z parallel to cube edges, body diagonals are the 3-fold axes
   */




  int ic, nc=cmplist.size();
  int nbasis;
  int i, j;
  int imin;
  double sum,summin;
  MatrixSq3<double> A, Ainv;
  Vector3<double> b, x;
  Vector3<double> u1_vec, u2_vec, u3_vec;
  std::string name;
  Vector< Vector3<double> > pos, ipos;
  CompoundStructureFit cmp;


#ifndef USE_SPGLIB


  double th;
  Vector3<double> origin(0);

  // Proper Cartesian rotations
  MatrixSq3<double> Rx60,  Ry60,  Rz60;
  MatrixSq3<double> Rx90,  Ry90,  Rz90;
  MatrixSq3<double> Rx120, Ry120, Rz120;
  MatrixSq3<double> Rx180, Ry180, Rz180;

  // Improper Cartesian rotations
  MatrixSq3<double> IRx60,  IRy60,  IRz60;
  MatrixSq3<double> IRx90,  IRy90,  IRz90;
  MatrixSq3<double> IRx120, IRy120, IRz120;
  MatrixSq3<double> IRx180, IRy180, IRz180;

  // Inversion:
  MatrixSq3<double> Inv;

  MatrixSq3<double> Mx, My, Mz; // same as IR*180(*)

  Vector3<double> dir_x, dir_y, dir_z;

  // Trigonal/hexagonal directions:
  Vector3<double> dir_h1, dir_h2, dir_h3;
  Vector3<double> dir_u1, dir_u2, dir_u3;
  // Planes:
  Vector3<double> dir_h1z, dir_h2z, dir_h3z;
  Vector3<double> dir_u1z, dir_u2z, dir_u3z;

  MatrixSq3<double> R180dh1, R180dh2, R180dh3;
  MatrixSq3<double> R180du1, R180du2, R180du3;

  MatrixSq3<double> IR180dh1, IR180dh2, IR180dh3;
  MatrixSq3<double> IR180du1, IR180du2, IR180du3;

  MatrixSq3<double> IR180dh1z, IR180dh2z, IR180dh3z;
  MatrixSq3<double> IR180du1z, IR180du2z, IR180du3z;

  // Cartesian directions:
  Vector3<double> dir_110, dir_bar110;

  // Cubic-cartesian directions:
  Vector3<double> dir_111, dir_bar111, dir_1bar11, dir_11bar1;

  MatrixSq3<double> R180d110, R180dbar110;
  MatrixSq3<double> IR180d110, IR180dbar110;

  MatrixSq3<double> R120d111, R120dbar111, R120d1bar11, R120d11bar1;
  MatrixSq3<double> IR120d111, IR120dbar111, IR120d1bar11, IR120d11bar1;


  std::cout << "Performing lattice symmetry analysis to obtain point group of compounds ..." << std::endl;


  get_rotation_matrix(Rx60, 2*PI/6, 0);
  get_rotation_matrix(Ry60, 2*PI/6, 1);
  get_rotation_matrix(Rz60, 2*PI/6, 2);
  get_rotation_matrix(Rx90, PI/2, 0);
  get_rotation_matrix(Ry90, PI/2, 1);
  get_rotation_matrix(Rz90, PI/2, 2);
  get_rotation_matrix(Rx120, 2*PI/3, 0);
  get_rotation_matrix(Ry120, 2*PI/3, 1);
  get_rotation_matrix(Rz120, 2*PI/3, 2);
  get_rotation_matrix(Rx180, PI, 0);
  get_rotation_matrix(Ry180, PI, 1);
  get_rotation_matrix(Rz180, PI, 2);

  get_improper_rotation_matrix(IRx60, 2*PI/6, 0);
  get_improper_rotation_matrix(IRy60, 2*PI/6, 1);
  get_improper_rotation_matrix(IRz60, 2*PI/6, 2);
  get_improper_rotation_matrix(IRx90, PI/2, 0);
  get_improper_rotation_matrix(IRy90, PI/2, 1);
  get_improper_rotation_matrix(IRz90, PI/2, 2);
  get_improper_rotation_matrix(IRx120, 2*PI/3, 0);
  get_improper_rotation_matrix(IRy120, 2*PI/3, 1);
  get_improper_rotation_matrix(IRz120, 2*PI/3, 2);
  get_improper_rotation_matrix(IRx180, PI, 0);
  get_improper_rotation_matrix(IRy180, PI, 1);
  get_improper_rotation_matrix(IRz180, PI, 2);

  Inv.elem(0,0) = Inv.elem(1,1) = Inv.elem(2,2) = -1.0; // Same as IR360.

  Mx = IRx180;
  My = IRy180;
  Mz = IRz180;


  dir_h1[0]=cos(0.0 * 2*PI);         dir_h1[1]=sin(0.0 * 2*PI);         dir_h1[2]=0;
  dir_h2[0]=cos(120.0/360.0 * 2*PI); dir_h2[1]=sin(120.0/360.0 * 2*PI); dir_h2[2]=0;
  dir_h3[0]=cos(240.0/360.0 * 2*PI); dir_h3[1]=sin(240.0/360.0 * 2*PI); dir_h3[2]=0;

  dir_u1[0]=cos( 30.0/360.0 * 2*PI); dir_u1[1]=sin( 30.0/360.0 * 2*PI); dir_u1[2]=0;
  dir_u2[0]=cos( 90.0/360.0 * 2*PI); dir_u2[1]=sin( 90.0/360.0 * 2*PI); dir_u2[2]=0;
  dir_u3[0]=cos(150.0/360.0 * 2*PI); dir_u3[1]=sin(150.0/360.0 * 2*PI); dir_u3[2]=0;

  dir_x[0]=1; dir_x[1]=0; dir_x[2]=0;
  dir_y[0]=0; dir_y[1]=1; dir_y[2]=0;
  dir_z[0]=0; dir_z[1]=0; dir_z[2]=1;

  vectorproduct(dir_h1, dir_z, dir_h1z);
  vectorproduct(dir_h2, dir_z, dir_h2z);
  vectorproduct(dir_h3, dir_z, dir_h3z);
  vectorproduct(dir_u1, dir_z, dir_u1z);
  vectorproduct(dir_u2, dir_z, dir_u2z);
  vectorproduct(dir_u3, dir_z, dir_u3z);

  dir_110[0]    =  1;
  dir_110[1]    =  1;
  dir_110[2]    =  0;
  dir_bar110[0] = -1;
  dir_bar110[1] =  1;
  dir_bar110[2] =  0;

  dir_111[0]    =  1;
  dir_111[1]    =  1;
  dir_111[2]    =  1;
  dir_bar111[0] = -1;
  dir_bar111[1] =  1;
  dir_bar111[2] =  1;
  dir_1bar11[0] =  1;
  dir_1bar11[1] = -1;
  dir_1bar11[2] =  1;
  dir_11bar1[0] =  1;
  dir_11bar1[1] =  1;
  dir_11bar1[2] = -1;


  th = PI;
  get_rotation_matrix_u(R180dh1, dir_h1, th);
  get_rotation_matrix_u(R180dh2, dir_h2, th);
  get_rotation_matrix_u(R180dh3, dir_h3, th);
  get_improper_rotation_matrix_u(IR180dh1, dir_h1, th);
  get_improper_rotation_matrix_u(IR180dh2, dir_h2, th);
  get_improper_rotation_matrix_u(IR180dh3, dir_h3, th);

  th = PI;
  get_rotation_matrix_u(R180du1, dir_u1, th);
  get_rotation_matrix_u(R180du2, dir_u2, th);
  get_rotation_matrix_u(R180du3, dir_u3, th);
  get_improper_rotation_matrix_u(IR180du1, dir_u1, th);
  get_improper_rotation_matrix_u(IR180du2, dir_u2, th);
  get_improper_rotation_matrix_u(IR180du3, dir_u3, th);

  th = PI;
  get_improper_rotation_matrix_u(IR180dh1z, dir_h1z, th);
  get_improper_rotation_matrix_u(IR180dh2z, dir_h2z, th);
  get_improper_rotation_matrix_u(IR180dh3z, dir_h3z, th);
  get_improper_rotation_matrix_u(IR180du1z, dir_u1z, th);
  get_improper_rotation_matrix_u(IR180du2z, dir_u2z, th);
  get_improper_rotation_matrix_u(IR180du3z, dir_u3z, th);

  th = PI;
  get_improper_rotation_matrix_u(IR180d110,    dir_110,    th);
  get_improper_rotation_matrix_u(IR180dbar110, dir_bar110, th);
  get_rotation_matrix_u(R180d110,   dir_110,   th);
  get_rotation_matrix_u(R180dbar110,dir_bar110,th);

  th = 2*PI/3;
  get_rotation_matrix_u(R120d111,    dir_111,    th);
  get_rotation_matrix_u(R120dbar111, dir_bar111, th);
  get_rotation_matrix_u(R120d1bar11, dir_1bar11, th);
  get_rotation_matrix_u(R120d11bar1, dir_11bar1, th);

  get_rotation_matrix_u(IR120d111,    dir_111,    th);
  get_rotation_matrix_u(IR120dbar111, dir_bar111, th);
  get_rotation_matrix_u(IR120d1bar11, dir_1bar11, th);
  get_rotation_matrix_u(IR120d11bar1, dir_11bar1, th);




  int i_cub=0, i_hex=0, i_trig=0, i_tetr=0, i_tric=0, i_mon=0, i_or=0;
  bool is_cubic, is_cubic1, is_cubic2;
  bool is_trig;
  bool is_hex;
  bool is_ort;
  bool is_tet;
  bool is_monoc;
  bool is_tric;
  int sinv;
  int s2x, s2y, s2z;
  int s3x, s3y, s3z;
  int s4x, s4y, s4z;
  int s6x, s6y, s6z;
  int smx, smy, smz;
  int sbar3x, sbar3y, sbar3z;
  int sbar4x, sbar4y, sbar4z;
  int sbar6x, sbar6y, sbar6z;
  int s3d111, s3dbar111, s3d1bar11, s3d11bar1;
  int sbar3d111, sbar3dbar111, sbar3d1bar11, sbar3d11bar1;
  int smd110, smdbar110;
  int s2d110, s2dbar110;
  int s2dh1, smdh1;
  int s2dh2, smdh2;
  int s2dh3, smdh3;
  int s2du1, smdu1;
  int s2du2, smdu2;
  int s2du3, smdu3;
  int smh1z, smu1z;
  int smh2z, smu2z;
  int smh3z, smu3z;
#endif



  for (ic=0; ic<nc; ++ic){
    // Local copy:
    cmp = cmplist[ic];
    cmp.csystem="unknown";
    cmp.csystem_sub=0;
    cmp.csymaxis = "";
    cmp.pointgroup = "1";
    cmp.spacegroup = "1";
    cmp.spacegroup_number = 0;

    if (!cmp.pbc[0] || !cmp.pbc[1] || !cmp.pbc[2]) continue;


    name = cmp.name;
    nbasis = cmp.nbasis;
    pos  = cmp.basis_vecs;
    ipos = cmp.basis_vecs;


    std::string dumpfile("symana-" + name + ".out");
    std::ofstream fout;
    fout.open(dumpfile.c_str());
    fout << "Symmetry analysis log for compound " << name << std::endl;


    //cout << "nbasis " << nbasis << endl;
    //cout << "types " << types << endl;


    // Find atom with smallest (most negative) coordinates:
    for (i=0; i<nbasis; ++i){
      sum = pos[i][0] + pos[i][1] + pos[i][2];
    
      if (i==0 || (i>0 && sum<summin)){
	imin   = i;
	summin = sum;
      }
    }

    // Find internal positions:
    A.elem(0,0)=cmp.u1_vec[0];
    A.elem(1,0)=cmp.u1_vec[1];
    A.elem(2,0)=cmp.u1_vec[2];
    A.elem(0,1)=cmp.u2_vec[0];
    A.elem(1,1)=cmp.u2_vec[1];
    A.elem(2,1)=cmp.u2_vec[2];
    A.elem(0,2)=cmp.u3_vec[0];
    A.elem(1,2)=cmp.u3_vec[1];
    A.elem(2,2)=cmp.u3_vec[2];
    for (i=0; i<nbasis; ++i){
      b = pos[i];
      A.solve(b,x,Ainv);
      ipos[i] = x;
      for (j=0; j<3; ++j){
	while (ipos[i][j]< 0.0) ipos[i][j]++;
	while (ipos[i][j]>=1.0) ipos[i][j]--;
      }
    }


    LatticeSimple lattice;

    lattice.nbasis   = nbasis;
    lattice.minpos   = pos[imin];
    lattice.origin   = Vector3<double>(0);
    lattice.avec  = cmp.u1_vec;
    lattice.bvec  = cmp.u2_vec;
    lattice.cvec  = cmp.u3_vec;
    lattice.pbc   = cmp.pbc;
    lattice.types = cmp.basis_types;
    lattice.pos   = pos;
    lattice.ipos  = ipos;
    /*
    cout << "lattice nbasis " << lattice.nbasis << endl;
    cout << "lattice origin " << lattice.origin << endl;
    cout << "lattice avec " << lattice.avec << endl;
    cout << "lattice bvec " << lattice.bvec << endl;
    cout << "lattice cvec " << lattice.cvec << endl;
    cout << "lattice pos " << lattice.pos << endl;
    cout << "lattice pbc " << lattice.pbc << endl;
    cout << "lattice types " << lattice.types << endl;
    */



#ifndef USE_SPGLIB


    sinv=0;
    s2x=0; s2y=0; s2z=0;
    s3x=0; s3y=0; s3z=0;
    s4x=0; s4y=0; s4z=0;
    s6x=0; s6y=0; s6z=0;
    smx=0; smy=0; smz=0;
    sbar3x=0; sbar3y=0; sbar3z=0;
    sbar4x=0; sbar4y=0; sbar4z=0;
    sbar6x=0; sbar6y=0; sbar6z=0;
    s3d111=0; s3dbar111=0; s3d1bar11=0; s3d11bar1=0;
    sbar3d111=0; sbar3dbar111=0; sbar3d1bar11=0; sbar3d11bar1=0;
    smd110=0; smdbar110=0;
    s2d110=0; s2dbar110=0;
    s2dh1=0; smdh1=0;
    s2dh2=0; smdh2=0;
    s2dh3=0; smdh3=0;
    s2du1=0; smdu1=0;
    s2du2=0; smdu2=0;
    s2du3=0; smdu3=0;
    smh1z=0; smu1z=0;
    smh2z=0; smu2z=0;
    smh3z=0; smu3z=0;



    LatticeSimple lattmp = lattice;
    /* lattmp.dump_xyz(name + "-tmp.xyz");
       lattmp.rotate(R120d11bar1);
       lattmp.dump_xyz(name + "-rot-tmp.xyz");
    */

    /* ################################################################
       Look for rotation symmetries using different shifts of atomic
       positions
       ################################################################
    */

    int ir, il, ish, nsh;
    Vector< Vector3<double> > shifts;
    Vector3<double> tv;

    shifts.push_back(tv);

    tv = -1.0 * lattice.avec;
    shifts.push_back(tv);

    tv = -1.0 * lattice.bvec;
    shifts.push_back(tv);

    tv = -1.0 * lattice.cvec;
    shifts.push_back(tv);

    tv = -0.5 * lattice.avec - 0.5 * lattice.bvec - 0.5 * lattice.cvec;
    shifts.push_back(tv);

    tv = -1.0 * lattice.minpos;
    shifts.push_back(tv);


    double CMpos[3];
    CMpos[0]=0;
    CMpos[1]=0;
    CMpos[2]=0;
    for (i=0; i<nbasis; ++i){
      for (j=0; j<3; ++j)
	CMpos[j] += pos[i][j];
    }
    for (i=0; i<3; ++i){
      CMpos[i] += 0.0*lattice.avec[i] + 0.0*lattice.bvec[i] + 0.0*lattice.cvec[i];
      CMpos[i] += 1.0*lattice.avec[i] + 0.0*lattice.bvec[i] + 0.0*lattice.cvec[i];
      CMpos[i] += 0.0*lattice.avec[i] + 1.0*lattice.bvec[i] + 0.0*lattice.cvec[i];
      CMpos[i] += 0.0*lattice.avec[i] + 0.0*lattice.bvec[i] + 1.0*lattice.cvec[i];
      CMpos[i] += 1.0*lattice.avec[i] + 1.0*lattice.bvec[i] + 0.0*lattice.cvec[i];
      CMpos[i] += 1.0*lattice.avec[i] + 0.0*lattice.bvec[i] + 1.0*lattice.cvec[i];
      CMpos[i] += 0.0*lattice.avec[i] + 1.0*lattice.bvec[i] + 1.0*lattice.cvec[i];
      CMpos[i] += 1.0*lattice.avec[i] + 1.0*lattice.bvec[i] + 1.0*lattice.cvec[i];
    }
    for (i=0; i<3; ++i){
      CMpos[i] /= nbasis + 8.0;
      tv[i] = -CMpos[i];
    }
    shifts.push_back(tv);

    nsh = shifts.size();


    /* ################################################################
       Look for rotation symmetries using different rotations
       of atomic positions
       ################################################################
    */

    for (ish=0; ish<nsh; ++ish){

      for (ir=0; ir<1; ++ir){

	LatticeSimple lattmp = lattice;
	/*
	lattmp.rotate(R120d11bar1);
	cout << lattmp.avec - lattice.avec << endl;
	cout << lattmp.bvec - lattice.bvec << endl;
	cout << lattmp.cvec - lattice.cvec << endl;
	for (int ijk=0; ijk<nbasis; ++ijk)
	  cout << lattmp.pos[ijk] - lattice.pos[ijk] << endl;

	exit(1);
	*/

	// Shift:
	lattmp.shift(shifts[ish]);

	// Rotation:
	MatrixSq3<double> Rdir;
	if (ir==0){
	  // Unit matrix:
	  for (int i1=0; i1<3; ++i1)
	    for (int i2=0; i2<3; ++i2)
	      Rdir.elem(i1,i2)=0.0;
	  Rdir.elem(0,0) = Rdir.elem(1,1) = Rdir.elem(2,2) = 1.0;
	}
	else if (ir==1){
	  // Make cvec coincide with z axis
	  get_matrix_for_rotation_to_coincide_with_axis(Rdir, dir_111, lattice.cvec);
	  lattmp.rotate(Rdir);
	}
	std::fout << "Using shift " << shifts[ish] << " and rotation " << Rdir << std::endl;





	if (lattmp.matrix_is_a_symm_op(Inv)){ sinv++; fout << "Compound " << name << " has symmetry Inv" << std::endl; }

	// ***********************************************************************
	// Proper rotations:
	// ***********************************************************************

	if (lattmp.matrix_is_a_symm_op(Rx90)){ s4x++; fout << "Compound " << name << " has symmetry Rot x 90deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Ry90)){ s4y++; fout << "Compound " << name << " has symmetry Rot y 90deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Rz90)){ s4z++; fout << "Compound " << name << " has symmetry Rot z 90deg" << std::endl; }

	if (lattmp.matrix_is_a_symm_op(Rx180)){ s2x++; fout << "Compound " << name << " has symmetry Rot x 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Ry180)){ s2y++; fout << "Compound " << name << " has symmetry Rot y 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Rz180)){ s2z++; fout << "Compound " << name << " has symmetry Rot z 180deg" << std::endl; }

	if (lattmp.matrix_is_a_symm_op(Rx60)){ s6x++; fout << "Compound " << name << " has symmetry Rot x 60deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Ry60)){ s6y++; fout << "Compound " << name << " has symmetry Rot y 60deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Rz60)){ s6z++; fout << "Compound " << name << " has symmetry Rot z 60deg" << std::endl; }

	if (lattmp.matrix_is_a_symm_op(Rx120)){ s3x++; fout << "Compound " << name << " has symmetry Rot x 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Ry120)){ s3y++; fout << "Compound " << name << " has symmetry Rot y 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Rz120)){ s3z++; fout << "Compound " << name << " has symmetry Rot z 120deg" << std::endl; }


	// ***********************************************************************
	// Improper rotations (Rotation + inversion in the origin):
	// ***********************************************************************
	/*
	if (lattmp.matrix_is_a_symm_op(Mx)){ smx++; fout << "Compound " << name << " has symmetry Mx (yz plane)" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(My)){ smy++; fout << "Compound " << name << " has symmetry My (xz plane)" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(Mz)){ smz++; fout << "Compound " << name << " has symmetry Mz (xy plane)" << std::endl; }
	*/
	if (lattmp.matrix_is_a_symm_op(IRx90)){ sbar4x++; fout << "Compound " << name << " has symmetry ImpRot x 90deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IRy90)){ sbar4y++; fout << "Compound " << name << " has symmetry ImpRot y 90deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IRz90)){ sbar4z++; fout << "Compound " << name << " has symmetry ImpRot z 90deg" << std::endl; }

	if (lattmp.matrix_is_a_symm_op(IRx180)){ smx++; fout << "Compound " << name << " has symmetry ImpRot x 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IRy180)){ smy++; fout << "Compound " << name << " has symmetry ImpRot y 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IRz180)){ smz++; fout << "Compound " << name << " has symmetry ImpRot z 180deg" << std::endl; }

	if (lattmp.matrix_is_a_symm_op(IRx60)){ sbar6x++; fout << "Compound " << name << " has symmetry ImpRot x 60deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IRy60)){ sbar6y++; fout << "Compound " << name << " has symmetry ImpRot y 60deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IRz60)){ sbar6z++; fout << "Compound " << name << " has symmetry ImpRot z 60deg" << std::endl; }

	if (lattmp.matrix_is_a_symm_op(IRx120)){ sbar3x++; fout << "Compound " << name << " has symmetry ImpRot x 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IRy120)){ sbar3y++; fout << "Compound " << name << " has symmetry ImpRot y 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IRz120)){ sbar3z++; fout << "Compound " << name << " has symmetry ImpRot z 120deg" << std::endl; }

  
	// ***********************************************************************
	// Rotations around cube space diagonals:
	// ***********************************************************************

	if (lattmp.matrix_is_a_symm_op(R120d111)){ s3d111++; fout << "Compound " << name << " has symmetry Rot (1,1,1) 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(R120dbar111)){ s3dbar111++; fout << "Compound " << name << " has symmetry Rot (-1,1,1) 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(R120d1bar11)){ s3d1bar11++; fout << "Compound " << name << " has symmetry Rot (1,-1,1) 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(R120d11bar1)){ s3d11bar1++; fout << "Compound " << name << " has symmetry Rot (1,1,-1) 120deg" << std::endl; }

	if (lattmp.matrix_is_a_symm_op(IR120d111)){ sbar3d111++; fout << "Compound " << name << " has symmetry Imp Rot (1,1,1) 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR120dbar111)){ sbar3dbar111++; fout << "Compound " << name << " has symmetry Imp Rot (-1,1,1) 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR120d1bar11)){ sbar3d1bar11++; fout << "Compound " << name << " has symmetry Imp Rot (1,-1,1) 120deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR120d11bar1)){ sbar3d11bar1++; fout << "Compound " << name << " has symmetry Imp Rot (1,1,-1) 120deg" << std::endl; }

	// ***********************************************************************
	// Improper rotations (reflections) around (in planes containing) cube side diagonals:
	// ***********************************************************************

	if (lattmp.matrix_is_a_symm_op(IR180d110)){ smd110++; fout << "Compound " << name << " has symmetry Imp Rot (1,1,0) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180dbar110)){ smdbar110++; fout << "Compound " << name << " has symmetry Imp Rot (-1,1,0) 180deg" << std::endl; }

	// ***********************************************************************
	// Rotations around cube side diagonals:
	// ***********************************************************************

	if (lattmp.matrix_is_a_symm_op(R180d110)){ s2d110++; fout << "Compound " << name << " has symmetry Rot (1,1,0) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(R180dbar110)){ s2dbar110++; fout << "Compound " << name << " has symmetry Rot (-1,1,0) 180deg" << std::endl; }

	// ***********************************************************************
	// Rotations relevant for trigonal/hexagonal systems:
	// ***********************************************************************
	// Directions: h1:=x, h2, h3, h4:=z
	// Angle between h1 and h2, h2 and h3, and h3 and h1 is 120 deg = 2*PI/3

	if (lattmp.matrix_is_a_symm_op(R180dh1)){ s2dh1++; fout << "Compound " << name << " has symmetry Rot (h1) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(R180dh2)){ s2dh2++; fout << "Compound " << name << " has symmetry Rot (h2) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(R180dh3)){ s2dh3++; fout << "Compound " << name << " has symmetry Rot (h3) 180deg" << std::endl; }

	// Axes orthogonal to h1,h2,h3:
	if (lattmp.matrix_is_a_symm_op(R180du1)){ s2du1++; fout << "Compound " << name << " has symmetry Rot (h1+30deg) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(R180du2)){ s2du2++; fout << "Compound " << name << " has symmetry Rot (h1+90deg) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(R180du3)){ s2du3++; fout << "Compound " << name << " has symmetry Rot (h1+150deg) 180deg" << std::endl; }



	// ***********************************************************************
	// Reflections in axes hi and ui, i=1,2,3.
	// Relevant for trigonal/hexagonal systems.
	// ***********************************************************************

	if (lattmp.matrix_is_a_symm_op(IR180dh1)){ smdh1++; fout << "Compound " << name << " has symmetry Imp Rot (h1) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180dh2)){ smdh2++; fout << "Compound " << name << " has symmetry Imp Rot (h2) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180dh3)){ smdh3++; fout << "Compound " << name << " has symmetry Imp Rot (h3) 180deg" << std::endl; }

	if (lattmp.matrix_is_a_symm_op(IR180du1)){ smdu1++; fout << "Compound " << name << " has symmetry Imp Rot (h1+30deg) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180du2)){ smdu2++; fout << "Compound " << name << " has symmetry Imp Rot (h1+90deg) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180du3)){ smdu3++; fout << "Compound " << name << " has symmetry Imp Rot (h1+150deg) 180deg" << std::endl; }



	// ***********************************************************************
	// Reflections in planes (hi, z) and (u1, z), i=1,2,3. Plane normal is e.g. h1 x z (vector product).
	// Relevant for trigonal/hexagonal systems.
	// ***********************************************************************
	if (lattmp.matrix_is_a_symm_op(IR180dh1z)){ smh1z++; fout << "Compound " << name << " has symmetry Imp Rot (h1 x z) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180dh2z)){ smh2z++; fout << "Compound " << name << " has symmetry Imp Rot (h2 x z) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180dh3z)){ smh3z++; fout << "Compound " << name << " has symmetry Imp Rot (h3 x z) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180du1z)){ smu1z++; fout << "Compound " << name << " has symmetry Imp Rot (u1 x z) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180du2z)){ smu2z++; fout << "Compound " << name << " has symmetry Imp Rot (u2 x z) 180deg" << std::endl; }
	if (lattmp.matrix_is_a_symm_op(IR180du3z)){ smu3z++; fout << "Compound " << name << " has symmetry Imp Rot (u3 x z) 180deg" << std::endl; }

      }

    }


    if (s6x && s3x) s3x=0; // 6-fold symmetry implies 3-fold symmetry
    if (s6y && s3y) s3y=0; // 6-fold symmetry implies 3-fold symmetry
    if (s6z && s3z) s3z=0; // 6-fold symmetry implies 3-fold symmetry

    if (s6x && s2x) s2x=0; // 6-fold symmetry implies 2-fold symmetry
    if (s6y && s2y) s2y=0; // 6-fold symmetry implies 2-fold symmetry
    if (s6z && s2z) s2z=0; // 6-fold symmetry implies 2-fold symmetry

    if (s4x && s2x) s2x=0; // 4-fold symmetry implies 2-fold symmetry
    if (s4y && s2y) s2y=0; // 4-fold symmetry implies 2-fold symmetry
    if (s4z && s2z) s2z=0; // 4-fold symmetry implies 2-fold symmetry

    if (sbar6x && s3x) s3x=0; // 6-fold symmetry implies 3-fold symmetry
    if (sbar6y && s3y) s3y=0; // 6-fold symmetry implies 3-fold symmetry
    if (sbar6z && s3z) s3z=0; // 6-fold symmetry implies 3-fold symmetry
    
    if (sbar4x && s2x) s2x=0; // 4-fold symmetry implies 2-fold symmetry
    if (sbar4y && s2y) s2y=0; // 4-fold symmetry implies 2-fold symmetry
    if (sbar4z && s2z) s2z=0; // 4-fold symmetry implies 2-fold symmetry
    


    // Summary:

    if (sinv) std::cout << "Compound " << name << " has symmetry Inv" << std::endl;
    if (s4x) std::cout << "Compound " << name << " has symmetry Rot x 90deg" << std::endl;
    if (s4y) std::cout << "Compound " << name << " has symmetry Rot y 90deg" << std::endl;
    if (s4z) std::cout << "Compound " << name << " has symmetry Rot z 90deg" << std::endl;
    if (s2x) std::cout << "Compound " << name << " has symmetry Rot x 180deg" << std::endl;
    if (s2y) std::cout << "Compound " << name << " has symmetry Rot y 180deg" << std::endl;
    if (s2z) std::cout << "Compound " << name << " has symmetry Rot z 180deg" << std::endl;
    if (s6x) std::cout << "Compound " << name << " has symmetry Rot x 60deg" << std::endl;
    if (s6y) std::cout << "Compound " << name << " has symmetry Rot y 60deg" << std::endl;
    if (s6z) std::cout << "Compound " << name << " has symmetry Rot z 60deg" << std::endl;
    if (s3x) std::cout << "Compound " << name << " has symmetry Rot x 120deg" << std::endl;
    if (s3y) std::cout << "Compound " << name << " has symmetry Rot y 120deg" << std::endl;
    if (s3z) std::cout << "Compound " << name << " has symmetry Rot z 120deg" << std::endl;
    if (smx) std::cout << "Compound " << name << " has symmetry Mx (yz plane)" << std::endl;
    if (smy) std::cout << "Compound " << name << " has symmetry My (xz plane)" << std::endl;
    if (smz) std::cout << "Compound " << name << " has symmetry Mz (xy plane)" << std::endl;
    if (sbar4x) std::cout << "Compound " << name << " has symmetry ImpRot x 90deg" << std::endl;
    if (sbar4y) std::cout << "Compound " << name << " has symmetry ImpRot y 90deg" << std::endl;
    if (sbar4z) std::cout << "Compound " << name << " has symmetry ImpRot z 90deg" << std::endl;
    if (smx) std::cout << "Compound " << name << " has symmetry ImpRot x 180deg" << std::endl;
    if (smy) std::cout << "Compound " << name << " has symmetry ImpRot y 180deg" << std::endl;
    if (smz) std::cout << "Compound " << name << " has symmetry ImpRot z 180deg" << std::endl;
    if (sbar6x) std::cout << "Compound " << name << " has symmetry ImpRot x 60deg" << std::endl;
    if (sbar6y) std::cout << "Compound " << name << " has symmetry ImpRot y 60deg" << std::endl;
    if (sbar6z) std::cout << "Compound " << name << " has symmetry ImpRot z 60deg" << std::endl;
    if (sbar3x) std::cout << "Compound " << name << " has symmetry ImpRot x 120deg" << std::endl;
    if (sbar3y) std::cout << "Compound " << name << " has symmetry ImpRot y 120deg" << std::endl;
    if (sbar3z) std::cout << "Compound " << name << " has symmetry ImpRot z 120deg" << std::endl;

    if (s3d111) std::cout << "Compound " << name << " has symmetry Rot (1,1,1) 120deg" << std::endl;
    if (s3dbar111) std::cout << "Compound " << name << " has symmetry Rot (-1,1,1) 120deg" << std::endl;
    if (s3d1bar11) std::cout << "Compound " << name << " has symmetry Rot (1,-1,1) 120deg" << std::endl;
    if (s3d11bar1) std::cout << "Compound " << name << " has symmetry Rot (1,1,-1) 120deg" << std::endl;

    if (sbar3d111) std::cout << "Compound " << name << " has symmetry Imp Rot (1,1,1) 120deg" << std::endl;
    if (sbar3dbar111) std::cout << "Compound " << name << " has symmetry Imp Rot (-1,1,1) 120deg" << std::endl;
    if (sbar3d1bar11) std::cout << "Compound " << name << " has symmetry Imp Rot (1,-1,1) 120deg" << std::endl;
    if (sbar3d11bar1) std::cout << "Compound " << name << " has symmetry Imp Rot (1,1,-1) 120deg" << std::endl;

    if (smd110) std::cout << "Compound " << name << " has symmetry Imp Rot (1,1,0) 180deg" << std::endl;
    if (smdbar110) std::cout << "Compound " << name << " has symmetry Imp Rot (-1,1,0) 180deg" << std::endl;

    if (s2d110) std::cout << "Compound " << name << " has symmetry Rot (1,1,0) 180deg" << std::endl;
    if (s2dbar110) std::cout << "Compound " << name << " has symmetry Rot (-1,1,0) 180deg" << std::endl;

    if (s2dh1) std::cout << "Compound " << name << " has symmetry Rot (h1) 180deg" << std::endl;
    if (s2dh2) std::cout << "Compound " << name << " has symmetry Rot (h2) 180deg" << std::endl;
    if (s2dh3) std::cout << "Compound " << name << " has symmetry Rot (h3) 180deg" << std::endl;

    if (s2du1) std::cout << "Compound " << name << " has symmetry Rot (h1+30deg) 180deg" << std::endl;
    if (s2du2) std::cout << "Compound " << name << " has symmetry Rot (h2+30deg) 180deg" << std::endl;
    if (s2du3) std::cout << "Compound " << name << " has symmetry Rot (h3+30deg) 180deg" << std::endl;

    if (smdh1) std::cout << "Compound " << name << " has symmetry Imp Rot (h1) 180deg" << std::endl;
    if (smdh2) std::cout << "Compound " << name << " has symmetry Imp Rot (h2) 180deg" << std::endl;
    if (smdh3) std::cout << "Compound " << name << " has symmetry Imp Rot (h3) 180deg" << std::endl;

    if (smdu1) std::cout << "Compound " << name << " has symmetry Imp Rot (h1+30deg) 180deg" << std::endl;
    if (smdu2) std::cout << "Compound " << name << " has symmetry Imp Rot (h2+30deg) 180deg" << std::endl;
    if (smdu3) std::cout << "Compound " << name << " has symmetry Imp Rot (h3+30deg) 180deg" << std::endl;

    if (smh1z) std::cout << "Compound " << name << " has symmetry Imp Rot (h1 x z) 180deg" << std::endl;
    if (smh2z) std::cout << "Compound " << name << " has symmetry Imp Rot (h2 x z) 180deg" << std::endl;
    if (smh3z) std::cout << "Compound " << name << " has symmetry Imp Rot (h3 x z) 180deg" << std::endl;

    if (smu1z) std::cout << "Compound " << name << " has symmetry Imp Rot (u1 x z) 180deg" << std::endl;
    if (smu2z) std::cout << "Compound " << name << " has symmetry Imp Rot (u2 x z) 180deg" << std::endl;
    if (smu3z) std::cout << "Compound " << name << " has symmetry Imp Rot (u3 x z) 180deg" << std::endl;




    // Now need to combine operations to figure out which crystal system each compound belongs to

    // Legend:
    // XY: axis X orthogonal to axis Y, rotation axes never parallel
    // Xm: axis X lies in mirror plane
    // X/m: axis X orthogonal to (any vector lying in the) mirror plane

    // Triclinic: only 1 and bar(1) symmetry
    // Monoclinic: 2, m, or 2/m,
    //   class 1 has symmetry axis in Cartesian Y direction
    //   class 2 has symmetry axis in Cartesian Z direction
    // Orthorombic: 222, mm2, or mmm
    // Trigonal:
    //   class 1: 3, bar(3)
    //   class 2: 32, 3m, bar(3) 2/m
    // Tetragonal:
    //   class 1: 4, bar(4), 4/m
    //   class 2: 422, 4mm, bar(4)2m, 4/m 2/m 2/m
    // Hexagonal: 6, bar(6), 6/m, 622, 6mm, bar(6)m2, 6/m 2/m 2/m
    // Cubic: 23, 2/m bar(3), 432, bar(4)3m, 4/m bar(3) 2/m



    is_cubic1 = false;
    is_cubic2 = false;

    if (s3d111 || s3dbar111 || s3d1bar11 || s3d11bar1)
      is_cubic1 = true;
    if (sbar3d111 || sbar3dbar111 || sbar3d1bar11 || sbar3d11bar1)
      is_cubic2 = true;

    is_cubic = is_cubic1 && is_cubic2;

    /*
    std::cout << "is_cubic1 " << is_cubic1 << std::endl;
    std::cout << "is_cubic2 " << is_cubic2 << std::endl;
    std::cout << "is_cubic  " << is_cubic  << std::endl;
    */

    is_trig = false;
    is_hex = false;
    is_ort = false;
    is_tet = false;
    is_monoc = false;
    is_tric = false;

    if (s3z || sbar3z) is_trig=true;
    if (s6z || sbar6z) is_hex=true;
    if (s4z || sbar4z) is_tet=true;
    /*
    if (is_cubic && is_trig) is_cubic=false;
    if (is_cubic && is_hex ) is_cubic=false;
    if (is_cubic && is_tet ) is_cubic=false;
    */
    if (is_trig && is_hex) is_trig=false;


    if (is_cubic==false){
      is_cubic1 = is_cubic2 = false;
    }


      
    cmp.csystem_sub = 0;




    // ##########################################################################
    // cubic
    // ##########################################################################
    if (is_cubic){

      if (is_cubic1){
	// .3.

	if ( (s4x && s4y && s4z) ||
	     (sbar4x && sbar4y && sbar4z) ){
	  // (4 or bar(4))3.
	  
	  if ( (s4x && s4y && s4z) &&
	       (s2d110 && s2dbar110) ){
	    // 432
	    cmp.csystem="cubic"; cmp.csystem_sub=3;
	    cmp.pointgroup = "432";
	  }
	  else if ( (sbar4x && sbar4y && sbar4z) &&
		    (smd110 && smdbar110) ){
	    // bar(4)3m
	    cmp.csystem="cubic"; cmp.csystem_sub=4;
	    cmp.pointgroup = "bar(4)3m";
	  }
	}
	else if (s2x && s2y && s2z){
	  // 23.
	  cmp.csystem="cubic"; cmp.csystem_sub=1;
	  cmp.pointgroup = "23";
	}
	else {
	  cmp.csystem="cubic"; cmp.csystem_sub=6;
	  cmp.pointgroup = ".3.";
	}
      }
      else if (is_cubic2){
	// .bar(3).
	
	if (smx && smy && smz){
	  // mbar(3).
	  
	  if ((s2x && s2y && s2z) &&
	      (smx && smy && smz)){
	    // 2/m bar(3)!!!!!
	    cmp.csystem="cubic"; cmp.csystem_sub=2;
	    cmp.pointgroup = "2/mbar(3)";
	  }
	  else if ( (s4x && s4y && s4z) &&
		    (smx && smy && smz) &&
		    (s2d110 && s2dbar110) &&
		    (smd110 && smdbar110) ){
	    // 4/m bar(3) 2/m
	    cmp.csystem="cubic"; cmp.csystem_sub=5;
	    cmp.pointgroup = "4/mbar(3)2/m";
	  }
	}
	else {
	  cmp.csystem="cubic"; cmp.csystem_sub=7;
	  cmp.pointgroup = ".bar(3).";
	}
      }
    }

    // ##########################################################################
    // tetragonal
    // ##########################################################################
    else if (is_tet){

      if (s4z && sbar4z==0){
	// 4..

	if (smz){

	  if ( (s2x && s2y) &&
	       (smx && smy) &&
	       (s2d110 && s2dbar110) &&
	       (smd110 && smdbar110) ){
	    // 4/m 2/m 2/m
	    cmp.csystem="tetragonal"; cmp.csystem_sub=7;
	    cmp.pointgroup = "4/m2/m2/m";
	  }
	  else {
	    // 4/m..
	    cmp.csystem="tetragonal"; cmp.csystem_sub=3;
	    cmp.pointgroup = "4/m";
	  }
	}
	else if ( (s2x && s2y) &&
		  (s2d110 && s2dbar110) ){
	  // 422
	  cmp.csystem="tetragonal"; cmp.csystem_sub=4;
	  cmp.pointgroup = "422";
	}
	else if ( (smx && smy) &&
		  (smd110 && smdbar110) ){
	  // 4mm
	  cmp.csystem="tetragonal"; cmp.csystem_sub=5;
	  cmp.pointgroup = "4mm";
	}
	else {
	  // 4..
	  cmp.csystem="tetragonal"; cmp.csystem_sub=1;
	  cmp.pointgroup = "4";
	}
      }
      else if (s4z==0 && sbar4z){
	// bar(4)..
	
	if ( (s2x && s2y) &&
	     (smd110 && smdbar110) ){
	  // bar(4)2m
	  cmp.csystem="tetragonal"; cmp.csystem_sub=6;
	  cmp.pointgroup = "bar(4)2m";
	}
	else {
	  // bar(4)
	  cmp.csystem="tetragonal"; cmp.csystem_sub=2;
	  cmp.pointgroup = "bar(4)";
	}
      }
    }

    // ##########################################################################
    // trigonal
    // ##########################################################################
    else if (is_trig){

      if (s3z && sbar3z==0){
	// 3..

	if (smdh1 && smdh2 && smdh3){
	  // 3m
	  cmp.csystem="trigonal"; cmp.csystem_sub=4;
	  cmp.pointgroup = "3m";
	}
	else if (s2dh1 && s2dh2 && s2dh3){
	  // 32
	  cmp.csystem="trigonal"; cmp.csystem_sub=3;
	  cmp.pointgroup = "32";
	}      
	else {
	  // 3
	  cmp.csystem="trigonal"; cmp.csystem_sub=1;
	  cmp.pointgroup = "3";
	}
      }
      else if (sbar3z && s3z==0){
	// bar(3)..

	if ( (s2dh1 && s2dh2 && s2dh3) &&
	     (smdh1 && smdh2 && smdh3) ){
	  // bar(3) 2/m
	  cmp.csystem="trigonal"; cmp.csystem_sub=5;
	  cmp.pointgroup = "bar(3)2/m";
	}
	else {
	  // bar(3)
	  cmp.csystem="trigonal"; cmp.csystem_sub=2;
	  cmp.pointgroup = "bar(3)";
	}
      }

    }

    // ##########################################################################
    // hexagonal
    // ##########################################################################
    else if (is_hex){

      if (s6z && sbar6z==0){
	// 6..

	if (smz){
	  // 6/m..

	  if ( (s2dh1 && s2dh2 && s2dh3) &&
	       (smdh1 && smdh2 && smdh3) &&
	       (s2du1 && s2du2 && s2du3) &&
	       (smdu1 && smdu2 && smdu3) ){
	    // 6/m 2/m 2/m
	    cmp.csystem="hexagonal"; cmp.csystem_sub=7;
	    cmp.pointgroup = "6/m2/m2/m";
	  }
	  else {
	    // 6/m
	    cmp.csystem="hexagonal"; cmp.csystem_sub=3;
	    cmp.pointgroup = "6/m";
	  }
	}
	else if ( (smdh1 && smdh2 && smdh3) &&
		  (smdu1 && smdu2 && smdu3) ){
	  // 6mm
	  cmp.csystem="hexagonal"; cmp.csystem_sub=5;
	  cmp.pointgroup = "6mm";
	}
	else if ( (s2dh1 && s2dh2 && s2dh3) &&
		  (s2du1 && s2du2 && s2du3) ){
	  // 622
	  cmp.csystem="hexagonal"; cmp.csystem_sub=4;
	  cmp.pointgroup = "622";
	}
	else {
	  // 6
	  cmp.csystem="hexagonal"; cmp.csystem_sub=1;
	  cmp.pointgroup = "6";
	}
	
      }
      else if (sbar6z && s6z==0){
	// bar(6)..

	if ( (smdh1 && smdh2 && smdh3) &&
	     (s2du1 && s2du2 && s2du3) ){
	  // bar(6)m2
	  cmp.csystem="hexagonal"; cmp.csystem_sub=6;
	  cmp.pointgroup = "bar(6)m2";
	}
	else {
	  // bar(6)
	  cmp.csystem="hexagonal"; cmp.csystem_sub=2;
	  cmp.pointgroup = "bar(6)";
	}

      }
    }

    // ##########################################################################
    // Cases: 2.., bar(2).., m..
    // ##########################################################################
    else if (s2x || smx || s2y || smy || s2z || smz){

      
      if ( (s2x && s2y && s2z) ||
	   (smx && smy && s2z) ||
	   (smx && smy && smz) ){
	// orthorombic
	
	if (s2x && s2y && s2z){
	  // 222
	  cmp.csystem="orthorombic"; cmp.csystem_sub=1;
	  cmp.pointgroup = "222";
	}
	else if (smx && smy && s2z){
	  // mm2
	  cmp.csystem="orthorombic"; cmp.csystem_sub=2;
	  cmp.pointgroup = "mm2";
	}
	else if ( (s2x && s2y && s2z) &&
		  (smx && smy && smz) ){
	  // 2/m 2/m 2/m
	  cmp.csystem="orthorombic"; cmp.csystem_sub=3;
	  cmp.pointgroup = "2/m2/m2/m";
	}

      }

      else if ( (s2z && s2y==0 && s2x==0) ||
		(s2y && s2z==0 && s2x==0) ||
		(s2x && s2y==0 && s2z==0) ||
		(smz && smy==0 && smx==0) ||
		(smy && smz==0 && smx==0) ||
		(smx && smy==0 && smz==0) ||
		(s2z && smz && s2y==0 && smy==0 && s2x==0 && smx==0) ||
		(s2y && smy && s2z==0 && smz==0 && s2x==0 && smx==0) ||
		(s2x && smx && s2y==0 && smy==0 && s2z==0 && smz==0) ){
	// monoclinic

	if (s2z && s2y==0 && s2x==0){ // 2
	  cmp.csystem="monoclinic"; cmp.csystem_sub=1;
	  cmp.csymaxis = "z";
	  cmp.pointgroup = "2";
	}
	else if (s2y && s2z==0 && s2x==0){ // 2
	  cmp.csystem="monoclinic"; cmp.csystem_sub=2;
	  cmp.csymaxis = "y";
	  cmp.pointgroup = "2";
	}
	else if (s2x && s2y==0 && s2z==0){ // 2
	  cmp.csystem="monoclinic"; cmp.csystem_sub=3;
	  cmp.pointgroup = "2";
	}
	else if (smz && smy==0 && smx==0){ // m
	  cmp.csystem="monoclinic"; cmp.csystem_sub=4;
	  cmp.csymaxis = "z";
	  cmp.pointgroup = "m";
	}
	else if (smy && smz==0 && smx==0){ // m
	  cmp.csystem="monoclinic"; cmp.csystem_sub=5;
	  cmp.csymaxis = "y";
	  cmp.pointgroup = "m";
	}
	else if (smx && smy==0 && smz==0){ // m
	  cmp.csystem="monoclinic"; cmp.csystem_sub=6;
	  cmp.pointgroup = "m";
	}
	else if (s2z && smz && s2y==0 && smy==0 && s2x==0 && smx==0){ // 2/m
	  cmp.csystem="monoclinic"; cmp.csystem_sub=7;
	  cmp.csymaxis = "z";
	  cmp.pointgroup = "2/m";
	}
	else if (s2y && smy && s2z==0 && smz==0 && s2x==0 && smx==0){ // 2/m
	  cmp.csystem="monoclinic"; cmp.csystem_sub=8;
	  cmp.csymaxis = "y";
	  cmp.pointgroup = "2/m";
	}
	else if (s2x && smx && s2y==0 && smy==0 && s2z==0 && smz==0){ // 2/m
	  cmp.csystem="monoclinic"; cmp.csystem_sub=9;
	  cmp.pointgroup = "2/m";
	}

      }
    }
    else if (sinv && cmp.csystem_sub<=0){
      // bar(1)
      cmp.csystem="triclinic"; cmp.csystem_sub=2;
      cmp.pointgroup = "bar(1)";
    }
    else if (cmp.csystem_sub<=0){
      cmp.csystem="triclinic"; cmp.csystem_sub=1;
      cmp.pointgroup = "1";
    }
#endif



#ifdef USE_SPGLIB
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################
    // spglib stuff
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################

    double lat[3][3];
    lat[0][0] = lattice.avec[0];
    lat[1][0] = lattice.avec[1];
    lat[2][0] = lattice.avec[2];
    lat[0][1] = lattice.bvec[0];
    lat[1][1] = lattice.bvec[1];
    lat[2][1] = lattice.bvec[2];
    lat[0][2] = lattice.cvec[0];
    lat[1][2] = lattice.cvec[1];
    lat[2][2] = lattice.cvec[2];

    int *stypes = new int [nbasis];
    
    typedef double (P3)[3];
    double (*position)[3]; // 'position' is a pointer to an object, object is of type '3-element-array of double'
    //position = new double [nbasis][3];
    position = new P3 [nbasis];

    for (i=0; i<nbasis; ++i){
      stypes[i] = lattice.types[i];
      position[i][0] = lattice.ipos[i][0];
      position[i][1] = lattice.ipos[i][1];
      position[i][2] = lattice.ipos[i][2];
    }

    lattice.csystem_sub = 1;
    cmp.csystem_sub = 1;

    int num_spg, num_atom = nbasis;
    char symbol[21];
    double tol = 1.0e-5;//sqrt( eps_d() );

    SpglibDataset *p_spg_dataset;
    p_spg_dataset = spg_get_dataset(lat, position, stypes, nbasis, 1e-5);
    std::cout << "spglib: spacegroup number   : " << p_spg_dataset->spacegroup_number << std::endl;
    fout << "spglib: spacegroup number   : " << p_spg_dataset->spacegroup_number << std::endl;

    lattice.spacegroup_number = p_spg_dataset->spacegroup_number;
    cmp.spacegroup_number = p_spg_dataset->spacegroup_number;

    // http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
    int sp=p_spg_dataset->spacegroup_number;
    if (sp<=0 || sp>230) aborterror("Spacegroup could not be found for compound " + name );

    std::string csystem;
    if      (sp>=1   && sp<=2  ) csystem="triclinic";
    else if (sp>=3   && sp<=15 ) csystem="monoclinic";
    else if (sp>=16  && sp<=74 ) csystem="orthorombic";
    else if (sp>=75  && sp<=142) csystem="tetragonal";
    else if (sp>=143 && sp<=167) csystem="trigonal";
    else if (sp>=168 && sp<=194) csystem="hexagonal";
    else if (sp>=195 && sp<=230) csystem="cubic";

    lattice.csystem = csystem;
    cmp.csystem = csystem;

    std::cout << "Compound " << name << " has crystal system " << csystem << std::endl;
    fout << "Compound " << name << " has crystal system " << csystem << std::endl;

    std::cout << "spglib: hall number         : " << p_spg_dataset->hall_number << std::endl;
    std::cout << "spglib: international symbol: " << p_spg_dataset->international_symbol << std::endl;
    std::cout << "spglib: hall symbol         : " << p_spg_dataset->hall_symbol << std::endl;
    std::cout << "spglib: setting             : " << p_spg_dataset->setting << std::endl;
    fout << "spglib: hall number         : " << p_spg_dataset->hall_number << std::endl;
    fout << "spglib: international symbol: " << p_spg_dataset->international_symbol << std::endl;
    fout << "spglib: hall symbol         : " << p_spg_dataset->hall_symbol << std::endl;
    fout << "spglib: setting             : " << p_spg_dataset->setting << std::endl;
    
    SpglibSpacegroupType spg_type = spg_get_spacegroup_type(p_spg_dataset->hall_number);
    std::cout << "spglib: international symbol (full) : " << spg_type.international_full << std::endl;
    std::cout << "spglib: international symbol (short): " << spg_type.international_short << std::endl;
    std::cout << "spglib: schoenflies                 : " << spg_type.schoenflies << std::endl;
    fout << "spglib: international symbol (full) : " << spg_type.international_full << std::endl;
    fout << "spglib: international symbol (short): " << spg_type.international_short << std::endl;
    fout << "spglib: schoenflies                 : " << spg_type.schoenflies << std::endl;

    lattice.spacegroup = spg_type.international_full;
    cmp.spacegroup = spg_type.international_full;

    
    // triclinic
    if (lattice.spacegroup_number>=  1 && lattice.spacegroup_number<=  1) lattice.pointgroup = "1";
    if (lattice.spacegroup_number>=  2 && lattice.spacegroup_number<=  2) lattice.pointgroup = "bar(3)";
    // monoclinic
    if (lattice.spacegroup_number>=  3 && lattice.spacegroup_number<=  5) lattice.pointgroup = "2";
    if (lattice.spacegroup_number>=  6 && lattice.spacegroup_number<=  9) lattice.pointgroup = "m";
    if (lattice.spacegroup_number>= 10 && lattice.spacegroup_number<= 15) lattice.pointgroup = "2/m";
    // orthorombic
    if (lattice.spacegroup_number>= 16 && lattice.spacegroup_number<= 24) lattice.pointgroup = "222";
    if (lattice.spacegroup_number>= 25 && lattice.spacegroup_number<= 46) lattice.pointgroup = "mm2";
    if (lattice.spacegroup_number>= 47 && lattice.spacegroup_number<= 74) lattice.pointgroup = "2/m2/m2/m";
    // tetragonal
    if (lattice.spacegroup_number>= 75 && lattice.spacegroup_number<= 80) lattice.pointgroup = "4";
    if (lattice.spacegroup_number>= 81 && lattice.spacegroup_number<= 82) lattice.pointgroup = "bar(4)";
    if (lattice.spacegroup_number>= 83 && lattice.spacegroup_number<= 88) lattice.pointgroup = "4/m";
    if (lattice.spacegroup_number>= 99 && lattice.spacegroup_number<=110) lattice.pointgroup = "4mm";
    if (lattice.spacegroup_number>=111 && lattice.spacegroup_number<=122) lattice.pointgroup = "bar(4)2m";
    if (lattice.spacegroup_number>= 89 && lattice.spacegroup_number<= 98) lattice.pointgroup = "422";
    if (lattice.spacegroup_number>=123 && lattice.spacegroup_number<=142) lattice.pointgroup = "4/m2/m2/m";
    // trigonal
    if (lattice.spacegroup_number>=143 && lattice.spacegroup_number<=146) lattice.pointgroup = "3";
    if (lattice.spacegroup_number>=147 && lattice.spacegroup_number<=148) lattice.pointgroup = "bar(3)";
    if (lattice.spacegroup_number>=149 && lattice.spacegroup_number<=155) lattice.pointgroup = "32";
    if (lattice.spacegroup_number>=162 && lattice.spacegroup_number<=167) lattice.pointgroup = "bar(3)2/m";
    if (lattice.spacegroup_number>=156 && lattice.spacegroup_number<=161) lattice.pointgroup = "3m";
    // hexagonal
    if (lattice.spacegroup_number>=168 && lattice.spacegroup_number<=173) lattice.pointgroup = "6";
    if (lattice.spacegroup_number>=174 && lattice.spacegroup_number<=174) lattice.pointgroup = "bar(6)";
    if (lattice.spacegroup_number>=175 && lattice.spacegroup_number<=176) lattice.pointgroup = "6/m";
    if (lattice.spacegroup_number>=177 && lattice.spacegroup_number<=182) lattice.pointgroup = "622";
    if (lattice.spacegroup_number>=183 && lattice.spacegroup_number<=186) lattice.pointgroup = "6mm";
    if (lattice.spacegroup_number>=187 && lattice.spacegroup_number<=190) lattice.pointgroup = "bar(6)m2";
    if (lattice.spacegroup_number>=191 && lattice.spacegroup_number<=194) lattice.pointgroup = "6/m2/m2/m";
    // cubic
    if (lattice.spacegroup_number>=195 && lattice.spacegroup_number<=199) lattice.pointgroup = "23";
    if (lattice.spacegroup_number>=200 && lattice.spacegroup_number<=206) lattice.pointgroup = "2/mbar(3)";
    if (lattice.spacegroup_number>=207 && lattice.spacegroup_number<=214) lattice.pointgroup = "432";
    if (lattice.spacegroup_number>=215 && lattice.spacegroup_number<=220) lattice.pointgroup = "bar(4)3m";
    if (lattice.spacegroup_number>=221 && lattice.spacegroup_number<=230) lattice.pointgroup = "4/mbar(3)2/m";

    cmp.pointgroup = lattice.pointgroup;



    char c = p_spg_dataset->setting[0];
    bool ax_is_small = fp_is_small_tol(lattice.avec[0], tol);
    bool ay_is_small = fp_is_small_tol(lattice.avec[1], tol);
    bool az_is_small = fp_is_small_tol(lattice.avec[2], tol);
    bool bx_is_small = fp_is_small_tol(lattice.bvec[0], tol);
    bool by_is_small = fp_is_small_tol(lattice.bvec[1], tol);
    bool bz_is_small = fp_is_small_tol(lattice.bvec[2], tol);
    bool cx_is_small = fp_is_small_tol(lattice.cvec[0], tol);
    bool cy_is_small = fp_is_small_tol(lattice.cvec[1], tol);
    bool cz_is_small = fp_is_small_tol(lattice.cvec[2], tol);
    if      (c=='a'){
      if      (ax_is_small && ay_is_small) cmp.csymaxis = "z";
      else if (ax_is_small && az_is_small) cmp.csymaxis = "y";
      else if (ay_is_small && az_is_small) cmp.csymaxis = "x";
      else cmp.csymaxis = "a";
    }
    else if (c=='b'){
      if      (bx_is_small && by_is_small) cmp.csymaxis = "z";
      else if (bx_is_small && bz_is_small) cmp.csymaxis = "y";
      else if (by_is_small && bz_is_small) cmp.csymaxis = "x";
      else cmp.csymaxis = "b";
    }
    else if (c=='c'){
      if      (cx_is_small && cy_is_small) cmp.csymaxis = "z";
      else if (cx_is_small && cz_is_small) cmp.csymaxis = "y";
      else if (cy_is_small && cz_is_small) cmp.csymaxis = "x";
    }


    if (cmp.csystem=="monoclinic"){
      if ( ! (cmp.csymaxis == "y" || cmp.csymaxis == "z") ){
	aborterror("ERROR: Compound " + name +
		   " should have main symmetry axis either 'y' or 'z'." +
		   " Now it is " + cmp.csymaxis);
      }
    }


    spg_free_dataset(p_spg_dataset);
    delete [] stypes;
    delete [] position;

    // ##########################################################################
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################

#endif


    std::cout << "crystal system " << cmp.csystem
	 << " subclass " << cmp.csystem_sub
	 << " point group " << cmp.pointgroup
	 << " symmetry axis " << cmp.csymaxis
	 << std::endl;
    fout << "crystal system " << cmp.csystem
	 << " subclass " << cmp.csystem_sub
	 << " point group " << cmp.pointgroup
	 << " symmetry axis " << cmp.csymaxis
	 << std::endl;


    fout.close();
    fout.clear();


    // Copy back:
    cmplist[ic] = cmp;

  }




}




