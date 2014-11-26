

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "atomsystem.hpp"
#include "exiterrors.hpp"
#include "constants.hpp"
#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-matrix3.hpp"
#include "utils-matrix.hpp"
#include "utils-string.hpp"
#include "utils-vector.hpp"

#include "compound.hpp"
#include "elem-iacs.hpp"
#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mdsettings.hpp"
#include "mtwister.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "potinfo.hpp"
#include "specs-fit-prop-pot.hpp"

using namespace std;
using namespace utils;
using namespace exiterrors;
using namespace constants;
using boost::format;
using std::ofstream;



MDSystem::MDSystem()
  : AtomSystem(),
    stresstensor_xyz(3,3,0.0),
    stresstensor_abc(3,3,0.0)
{
  name = "none";
  p_potinfo = 0;
  debug_creation = false;
  debug_mds = false;
  N[0]=0;
  N[1]=0;
  N[2]=0;
  vel.cap(100);
  acc.cap(100);
  frc.cap(100);
  frc_num.cap(100);
  virials.cap(100);
  Ep.cap(100);
  Ek.cap(100);
  dpos.cap(100);
  pos_int_tmp.cap(100);
  pos_int_ini.cap(100);
  pos_int_fin.cap(100);

  vel.resize(0);
  acc.resize(0);
  frc.resize(0);
  frc_num.resize(0);
  virials.resize(0);
  Ep.resize(0);
  Ek.resize(0);
  dpos.resize(0);
  pos_int_tmp.resize(0);
  pos_int_ini.resize(0);
  pos_int_fin.resize(0);

  rcut_max = dt = T = T_at_quench_start = V = P = Px = Py = Pz = 0;
  Ep_tot = Ek_tot = P_max = F_max = displ_max = 0;
  get_pot_force = get_pot_energy = false;

  iac_pure_ABOP = false;
  iac_pure_EAM  = false;
  sys_single_elem = false;

}






// ##############################################################################

/*
  If Ni> 0 make cell according to Ni numbers.

  If Ni<=0 and distmin>0 make cell so that boxlen in direction ai is >= distmin.
  If Ni<=0 and distmin<0 make small cell (N1=N2=N3=1).
*/

void MDSystem::create_from_structure(CompoundStructure & cmp,
				     double distmin // only makes sense for periodic directions
				     ){
  int i,j,k,p,iat;
  Vector<double> tv;
  Vector<double> origin(3, 0);



  pbc = cmp.pbc;

  if (N[0]<=0){
    N[0] = 1;
    if (pbc[0] && distmin>0){
      N[0] = distmin/cmp.u1_vec.magn();
      while (N[0] * cmp.u1_vec.magn() <= distmin) N[0]++;
    }
  }
  if (N[1]<=0){
    N[1] = 1;
    if (pbc[1] && distmin>0){
      N[1] = distmin/cmp.u2_vec.magn();
      while (N[1] * cmp.u2_vec.magn() <= distmin) N[1]++;
    }
  }
  if (N[2]<=0){
    N[2] = 1;
    if (pbc[2] && distmin>0){
      N[2] = distmin/cmp.u3_vec.magn();
      while (N[2] * cmp.u3_vec.magn() <= distmin) N[2]++;
    }
  }


  cout << "Creating MD system: Using N[0] N[1] N[2]  " << N[0] << " " << N[1] << " " << N[2] << endl;

  cout << "Creating MD system: scalefactor  " << cmp.scalefactor << endl;
  cout << "Creating MD system: internal lattice parameters  " 
       << cmp.lpa << " "
       << cmp.lpb << " "
       << cmp.lpc << endl;


  name = cmp.name;




  tv = cmp.u1_vec; tv.normalize(); set_boxdir(0, tv);
  tv = cmp.u2_vec; tv.normalize(); set_boxdir(1, tv);
  tv = cmp.u3_vec; tv.normalize(); set_boxdir(2, tv);

  boxlen[0] = N[0] * cmp.u1_vec.magn();
  boxlen[1] = N[1] * cmp.u2_vec.magn();
  boxlen[2] = N[2] * cmp.u3_vec.magn();

  update_box_geometry();

  cout << "Creating MD system: Box lengths: " << boxlen[0] << " " << boxlen[1] << " " << boxlen[2] << endl;

  // Default:
  origin[0] = -0.5 * boxlen[0];
  origin[1] = -0.5 * boxlen[1];
  origin[2] = -0.5 * boxlen[2];

  if (! cmp.use_readin_structure){
    cmp.origin_from_model(N[0], N[1], N[2]);
    origin = cmp.origin;
  }
  if (cmp.use_origin_spec){
    origin = cmp.origin;
  }


  cout << "Creating MD system: Origin: " << origin[0] << " " << origin[1] << " " << origin[2] << endl;







  /*
    cout << "boxlen 0: " << boxlen[0] << endl;
    cout << "boxlen 1: " << boxlen[1] << endl;
    cout << "boxlen 2: " << boxlen[2] << endl;

    cout << "boxdir 0: " << get_boxdir(0) << endl;
    cout << "boxdir 1: " << get_boxdir(1) << endl;
    cout << "boxdir 2: " << get_boxdir(2) << endl;

    cout << "N[0] N[1] N[2]: " << N[0] << " " << N[1] << " " << N[2] << endl;
  */

  // Create atom system:
  clear_all_atoms();


  for (i=0; i<N[0]; ++i){
    for (j=0; j<N[1]; ++j){
      for (k=0; k<N[2]; ++k){
	for (p=0; p<cmp.nbasis; ++p){
	  iat = add_atom();

	  pos[iat] = origin
	    + i * cmp.u1_vec + j * cmp.u2_vec + k * cmp.u3_vec
	    + cmp.basis_vecs[p];

	  // Atom type info will be filled in later. Now just put something here:
	  type[iat] = 0;
	  idx[iat] = iat;
	  matter[iat] = cmp.basis_elems[p];
	}
      }
    }
  }

  //  if (debug_creation)
  //cout << "Created system of " << pos.size() << " atoms." << endl;

  return;
}


















#include "mdsystem-relax.cppinc"




void MDSystem::transform_cell(const Matrix<double> & alpha_cart,
			      const double lowlim){
  int i, nat = natoms();
  Vector<double> boxlen_new(3,0), tv(3,0);
  Matrix<double> boxdir_new(3,3,0);

  double pf1=-0.5, pf2=0.5;
  if (lowlim>0){
    pf1 = 0.0;
    pf2 = 1.0;
  }

  update_box_geometry();

  if (pos_int_tmp.size()!=nat) pos_int_tmp.resize(nat);
  for (i=0; i<nat; ++i){
    if (pos_int_tmp[i].size()!=3) pos_int_tmp[i].resize(3);
  }


  // Get internal positions, with box lenghts scaled away.
#pragma omp parallel for schedule(static)
  for (i=0; i<nat; ++i){

    // Original:
    //pos_int_tmp[i] = Bravaismatrix_inv * pos[i];

    pos_int_tmp[i][0] = pos[i][0];
    pos_int_tmp[i][1] = pos[i][1];
    pos_int_tmp[i][2] = pos[i][2];
    if (! isCart){
      pos_int_tmp[i][0] = 0.0
	+ Bravaismatrix_inv.elem(0,0) * pos[i][0]
	+ Bravaismatrix_inv.elem(0,1) * pos[i][1]
	+ Bravaismatrix_inv.elem(0,2) * pos[i][2];
      pos_int_tmp[i][1] = 0.0
	+ Bravaismatrix_inv.elem(1,0) * pos[i][0]
	+ Bravaismatrix_inv.elem(1,1) * pos[i][1]
	+ Bravaismatrix_inv.elem(1,2) * pos[i][2];
      pos_int_tmp[i][2] = 0.0
	+ Bravaismatrix_inv.elem(2,0) * pos[i][0]
	+ Bravaismatrix_inv.elem(2,1) * pos[i][1]
	+ Bravaismatrix_inv.elem(2,2) * pos[i][2];
    }

    pos_int_tmp[i][0] /= boxlen[0];
    pos_int_tmp[i][1] /= boxlen[1];
    pos_int_tmp[i][2] /= boxlen[2];

  }

  // ***************************************************************
  // Transform box:
  // ***************************************************************

  // Multiply in the box lenghts so the direction vectors are complete:
  tv = boxlen[0] * boxdir.col(0); boxdir.col(0, tv);
  tv = boxlen[1] * boxdir.col(1); boxdir.col(1, tv);
  tv = boxlen[2] * boxdir.col(2); boxdir.col(2, tv);

  // Transform direction vectors:
  tv = alpha_cart * boxdir.col(0); boxdir_new.col(0, tv);
  tv = alpha_cart * boxdir.col(1); boxdir_new.col(1, tv);
  tv = alpha_cart * boxdir.col(2); boxdir_new.col(2, tv);

  // Normalize the direction vectors and retain the sizes as the new box lengths.
  tv = boxdir_new.col(0); boxlen_new[0] = tv.normalize(); boxdir_new.col(0, tv);
  tv = boxdir_new.col(1); boxlen_new[1] = tv.normalize(); boxdir_new.col(1, tv);
  tv = boxdir_new.col(2); boxlen_new[2] = tv.normalize(); boxdir_new.col(2, tv);

  boxlen = boxlen_new;
  boxdir = boxdir_new;

  update_box_geometry();

  // Get new atomic positions:
#pragma omp parallel for schedule(static)
  for (i=0; i<nat; ++i){

    pos_int_tmp[i][0] *= boxlen[0];
    pos_int_tmp[i][1] *= boxlen[1];
    pos_int_tmp[i][2] *= boxlen[2];

    while (pbc[0] && (pos_int_tmp[i][0] <  pf1 * boxlen[0])) pos_int_tmp[i][0] += boxlen[0];
    while (pbc[0] && (pos_int_tmp[i][0] >= pf2 * boxlen[0])) pos_int_tmp[i][0] -= boxlen[0];

    while (pbc[1] && (pos_int_tmp[i][1] <  pf1 * boxlen[1])) pos_int_tmp[i][1] += boxlen[1];
    while (pbc[1] && (pos_int_tmp[i][1] >= pf2 * boxlen[1])) pos_int_tmp[i][1] -= boxlen[1];

    while (pbc[2] && (pos_int_tmp[i][2] <  pf1 * boxlen[2])) pos_int_tmp[i][2] += boxlen[2];
    while (pbc[2] && (pos_int_tmp[i][2] >= pf2 * boxlen[2])) pos_int_tmp[i][2] -= boxlen[2];

    pos[i] = boxdir * pos_int_tmp[i];
  }


}





void MDSystem::handle_pbc_of_positions(const double lowlim){
  int i, nat = natoms();
  double pf1=-0.5, pf2=0.5;
  if (lowlim>0){
    pf1 = 0.0;
    pf2 = 1.0;
  }


#pragma omp parallel for schedule(static)
  for (i=0; i<nat; ++i){

    double drc0, drc1, drc2;
    double drs0, drs1, drs2;

    drc0 = pos[i][0];
    drc1 = pos[i][1];
    drc2 = pos[i][2];

    // Get distance in skew coordinate system, where periodics can be checked:
    drs0 = drc0;
    drs1 = drc1;
    drs2 = drc2;
    if (! isCart){
      drs0 = 0.0
	+ Bravaismatrix_inv.elem(0,0) * drc0
	+ Bravaismatrix_inv.elem(0,1) * drc1
	+ Bravaismatrix_inv.elem(0,2) * drc2;
      drs1 = 0.0
	+ Bravaismatrix_inv.elem(1,0) * drc0
	+ Bravaismatrix_inv.elem(1,1) * drc1
	+ Bravaismatrix_inv.elem(1,2) * drc2;
      drs2 = 0.0
	+ Bravaismatrix_inv.elem(2,0) * drc0
	+ Bravaismatrix_inv.elem(2,1) * drc1
	+ Bravaismatrix_inv.elem(2,2) * drc2;
    }
      
    // Periodics check:
    while (pbc[0] && drs0 <  pf1 * boxlen[0]) drs0 += boxlen[0];
    while (pbc[0] && drs0 >= pf2 * boxlen[0]) drs0 -= boxlen[0];
    while (pbc[1] && drs1 <  pf1 * boxlen[1]) drs1 += boxlen[1];
    while (pbc[1] && drs1 >= pf2 * boxlen[1]) drs1 -= boxlen[1];
    while (pbc[2] && drs2 <  pf1 * boxlen[2]) drs2 += boxlen[2];
    while (pbc[2] && drs2 >= pf2 * boxlen[2]) drs2 -= boxlen[2];

    // Get distance in Cartesian coordinate system:
    drc0 = drs0;
    drc1 = drs1;
    drc2 = drs2;
    if (! isCart){
      drc0 = boxdir.elem(0,0) * drs0 + boxdir.elem(0,1) * drs1 + boxdir.elem(0,2) * drs2;
      drc1 = boxdir.elem(1,0) * drs0 + boxdir.elem(1,1) * drs1 + boxdir.elem(1,2) * drs2;
      drc2 = boxdir.elem(2,0) * drs0 + boxdir.elem(2,1) * drs1 + boxdir.elem(2,2) * drs2;
    }
    pos[i][0] = drc0;
    pos[i][1] = drc1;
    pos[i][2] = drc2;
      
  }


}



void MDSystem::calc_closepacked_volume(){
  V = 0;
    
  double drsqmin=100,drsq=0;
  Vector<double> posi(3), posj(3), dr(3);
  int i,j,ij,counter=0, nat = natoms();



  for (i=0; i<nat; ++i){
    posi[0] = pos[i][0];
    posi[1] = pos[i][1];
    posi[2] = pos[i][2];
      
    for (ij=0; ij<neighborcollection[i].size(); ++ij){
      j = neighborcollection[i][ij];
      posj[0] = pos[j][0];
      posj[1] = pos[j][1];
      posj[2] = pos[j][2];

      if (isCart){
	drsq = 0.0;
	dr[0] = posi[0] - posj[0];
	while (pbc[0] && dr[0]<-0.5*boxlen[0]) dr[0] += boxlen[0];
	while (pbc[0] && dr[0]>=0.5*boxlen[0]) dr[0] -= boxlen[0];
	drsq += dr[0] * dr[0];
	dr[1] = posi[1] - posj[1];
	while (pbc[1] && dr[1]<-0.5*boxlen[1]) dr[1] += boxlen[1];
	while (pbc[1] && dr[1]>=0.5*boxlen[1]) dr[1] -= boxlen[1];
	drsq += dr[1] * dr[1];
	dr[2] = posi[2] - posj[2];
	while (pbc[2] && dr[2]<-0.5*boxlen[2]) dr[2] += boxlen[2];
	while (pbc[2] && dr[2]>=0.5*boxlen[2]) dr[2] -= boxlen[2];
	drsq += dr[2] * dr[2];
      }
      else {
	get_atom_distance_vec(posi, posj, dr);
	drsq = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      }
      if (counter==0 || (counter>0 && drsqmin<drsq)){
	drsqmin = drsq;
	counter++;
      }
    }

    V += 4*PI*drsqmin * sqrt(drsqmin)/3.0;
  }

  return;
}







#include "mdsystem-force.cppinc"








#if 0
  for (i=0; i<nat; ++i){
    // Get internal positions at end of simulation:
    int_pos_matrix.col(0, boxdir.col(0));
    int_pos_matrix.col(1, boxdir.col(1));
    int_pos_matrix.col(2, boxdir.col(2));

    int_pos_matrix.solve(pos[i], pos_int_fin[i], dummy_matrix);

    // Scale internal positions at start with current box lenghts:
    pos_int_ini[i][0] *= boxlen[0];
    pos_int_ini[i][1] *= boxlen[1];
    pos_int_ini[i][2] *= boxlen[2];


    for (j=0; j<3; ++j){
      tv1[j] = pos_int_fin[i][j] - pos_int_ini[i][j];
      if (pbc[j]){
	while (tv1[j] <= 0.0)       tv1[j] += boxlen[j];
	while (tv1[j] >  boxlen[j]) tv1[j] -= boxlen[j];
      }
    }
    tv2 = boxdir * tv1;
    td = tv2.magn();
    if (i==0 || (i>0 && (td>displ_max))) displ_max = td;
  }
#endif
