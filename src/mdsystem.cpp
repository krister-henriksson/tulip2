

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
#include "constants.hpp"
#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-matrix3.hpp"
#include "utils-matrix.hpp"
#include "utils-string.hpp"
#include "utils-vector.hpp"
#include "utils-errors.hpp"

#include "compound.hpp"
#include "elem-iacs.hpp"
//#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mdsettings.hpp"
#include "mtwister.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "potinfo.hpp"
//#include "specs-fit-prop-pot.hpp"
#include "errors.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;




using namespace utils;
using namespace constants;
using boost::format;




MDSystem::MDSystem()
  : AtomSystem(),
    stresstensor_xyz(0.0),
    stresstensor_abc(0.0),
    pos_CM(0),
    vel_CM(0)
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
  posint_ini.cap(100);
  posint_fin.cap(100);

  atom_is_fixed.cap(100);
  atom_freedir.cap(100);
  atom_freeplane.cap(100);

  atom_is_fixed.resize(0);
  atom_freedir.resize(0);
  atom_freeplane.resize(0);

  vel.resize(0);
  acc.resize(0);
  frc.resize(0);
  frc_num.resize(0);
  virials.resize(0);
  Ep.resize(0);
  Ek.resize(0);
  dpos.resize(0);
  pos_int_tmp.resize(0);
  posint_ini.resize(0);
  posint_fin.resize(0);

  rcut_max = dt = T = T_at_quench_start = P = Px = Py = Pz = 0;
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
  Vector3<double> tv;
  Vector3<double> origin(0);



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


  if (cmp.Ndesired[0] > 0 && cmp.Ndesired[0] > N[0]) N[0] = cmp.Ndesired[0];
  if (cmp.Neven_desired[0]){ if (N[0] % 2 != 0) N[0]++; }
  if (cmp.Nodd_desired[0] ){ if (N[0] % 2 == 0) N[0]++; }

  if (cmp.Ndesired[1] > 0 && cmp.Ndesired[1] > N[1]) N[1] = cmp.Ndesired[1];
  if (cmp.Neven_desired[1]){ if (N[1] % 2 != 0) N[1]++; }
  if (cmp.Nodd_desired[1] ){ if (N[1] % 2 == 0) N[1]++; }

  if (cmp.Ndesired[2] > 0 && cmp.Ndesired[2] > N[2]) N[2] = cmp.Ndesired[2];
  if (cmp.Neven_desired[2]){ if (N[2] % 2 != 0) N[2]++; }
  if (cmp.Nodd_desired[2] ){ if (N[2] % 2 == 0) N[2]++; }



  std::cout << "Creating MD system: name " << cmp.name << std::endl;
  std::cout << "Creating MD system: Using N[0] N[1] N[2]  " << N[0] << " " << N[1] << " " << N[2] << std::endl;

  std::cout << "Creating MD system: scalefactor  " << cmp.scalefactor << std::endl;
  std::cout << "Creating MD system: internal lattice parameters  " 
       << cmp.lpa << " "
       << cmp.lpb << " "
       << cmp.lpc << std::endl;


  name = cmp.name;




  tv = cmp.u1_vec; tv.normalize(); set_boxdir(0, tv);
  tv = cmp.u2_vec; tv.normalize(); set_boxdir(1, tv);
  tv = cmp.u3_vec; tv.normalize(); set_boxdir(2, tv);

  boxlen[0] = N[0] * cmp.u1_vec.magn();
  boxlen[1] = N[1] * cmp.u2_vec.magn();
  boxlen[2] = N[2] * cmp.u3_vec.magn();

  update_box_geometry();

  std::cout << "Creating MD system: Box lengths: " << boxlen[0] << " " << boxlen[1] << " " << boxlen[2] << std::endl;

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


  std::cout << "Creating MD system: Origin: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;







  /*
    std::cout << "boxlen 0: " << boxlen[0] << std::endl;
    std::cout << "boxlen 1: " << boxlen[1] << std::endl;
    std::cout << "boxlen 2: " << boxlen[2] << std::endl;

    std::cout << "boxdir 0: " << get_boxdir(0) << std::endl;
    std::cout << "boxdir 1: " << get_boxdir(1) << std::endl;
    std::cout << "boxdir 2: " << get_boxdir(2) << std::endl;

    std::cout << "N[0] N[1] N[2]: " << N[0] << " " << N[1] << " " << N[2] << std::endl;
  */

  // Create atom system:
  clear_all_atoms();


  int atom_counter=1;
  for (i=0; i<N[0]; ++i){
    for (j=0; j<N[1]; ++j){
      for (k=0; k<N[2]; ++k){
	for (p=0; p<cmp.nbasis; ++p){
	  iat = add_atom();

	  pos[iat] = origin
	    + i * cmp.u1_vec + j * cmp.u2_vec + k * cmp.u3_vec
	    + cmp.basis_vecs[p];

	  // Atom type info will be filled in later. Now just put something here:
	  sitetype[iat] = p;
	  type[iat] = 0;
	  itype[iat] = 0;
	  idx[iat] = atom_counter++;
	  matter[iat] = cmp.basis_elems[p];

	  atom_is_fixed[iat]  = cmp.basis_is_fixed[p];
	  //if (cmp.basis_freedir[p].size()==3)   atom_freedir[iat]   = cmp.basis_freedir[p];
	  //if (cmp.basis_freeplane[p].size()==3) atom_freeplane[iat] = cmp.basis_freeplane[p];
	  atom_freedir[iat]   = cmp.basis_freedir[p];
	  atom_freeplane[iat] = cmp.basis_freeplane[p];

	}
      }
    }
  }

  //  if (debug_creation)
  std::cout << "Created system of " << pos.size() << " atoms." << std::endl;

  return;
}




void MDSystem::transform_cell(const MatrixSq3<double> & alpha_cart,
			      const double lowlim){
  int i, nat = natoms();
  Vector3<double> boxlen_new(0), tv(0);
  MatrixSq3<double> boxdir_new(0);

  double pf1=-0.5, pf2=0.5;
  if (lowlim>0){
    pf1 = 0.0;
    pf2 = 1.0;
  }

  // Make sure box direction vectors are normalized and that the
  // Bravais matrix is up to date:
  update_box_geometry();

  if (pos_int_tmp.size()!=nat) pos_int_tmp.resize(nat);
  for (i=0; i<nat; ++i){
    pos_int_tmp[i] = Vector3<double>(0);
  }


  // Get internal positions, with box lenghts scaled away.

  // #pragma omp parallel for schedule(static)
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
  // #pragma omp parallel for schedule(static)
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



void MDSystem::get_CM_pos(void){
  handle_pbc_of_positions();
  get_CM_vec(pos_CM, 1);
}

void MDSystem::get_CM_vel(void){
  get_CM_vec(vel_CM, 2);
}

void MDSystem::get_CM_vec(Vector3<double> & vec, const int mode){
  vec[0] = 0.0;
  vec[1] = 0.0;
  vec[2] = 0.0;

  double mass_cm = 0.0;
  double mass;

  mass = elem.mass(matter[0]);

  if (mode==1){
    if (sys_single_elem){
      for (int i=0; i<natoms(); ++i){
	mass_cm += mass;
	vec[0] += mass * pos[i][0];
	vec[1] += mass * pos[i][1];
	vec[2] += mass * pos[i][2];
      }
    }
    else {
      for (int i=0; i<natoms(); ++i){
	mass = elem.mass( itype[i] ); //elem.name2idx(matter[i]);
	mass_cm += mass;
	vec[0] += mass * pos[i][0];
	vec[1] += mass * pos[i][1];
	vec[2] += mass * pos[i][2];
      }
    }
  }
  else if (mode==2){
    if (sys_single_elem){
      for (int i=0; i<natoms(); ++i){
	mass_cm += mass;
	vec[0] += mass * vel[i][0];
	vec[1] += mass * vel[i][1];
	vec[2] += mass * vel[i][2];
      }
    }
    else {
      for (int i=0; i<natoms(); ++i){
	mass = elem.mass( itype[i] ); //elem.name2idx(matter[i]);
	mass_cm += mass;
	vec[0] += mass * vel[i][0];
	vec[1] += mass * vel[i][1];
	vec[2] += mass * vel[i][2];
      }
    }
  }

  vec[0] /= mass_cm;
  vec[1] /= mass_cm;
  vec[2] /= mass_cm;
  
}

