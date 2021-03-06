



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

//#include "compound.hpp"
#include "elem-iacs.hpp"
//#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mdsettings.hpp"
//#include "mtwister.hpp"
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






void MDSystem::calc_forces_and_energies(){
  get_pot_force  = true;
  get_pot_energy = true;

  get_forces_and_energies_common();
  return;
}



double MDSystem::calc_potential_energy(){
  get_pot_force  = false;
  get_pot_energy = true;

  return get_forces_and_energies_common();
}



double MDSystem::get_forces_and_energies_common(){
  int i, j, k, t1, t2, t1_old;
  std::string s1, s2;
  int n_ABOP=0, n_EAM=0;
  int n_iacs=0, nat = natoms();


  /*
  iac_pure_abop
  iac_pure_eam
  sys_single_elem
*/

  int n2i = elem.name2idx( matter[0] );

  sys_single_elem=true;
  t1_old = n2i;

  // ######################################################################
  // ######################################################################
  // Clear all energies, forces, and virials:
  // ######################################################################
  // ######################################################################
  for (i=0; i<nat; ++i){
    Ep[i] = 0;
    Ek[i] = 0;
    frc[i][0] = 0;
    frc[i][1] = 0;
    frc[i][2] = 0;

    virials[i].elem(0,0) = 0.0;
    virials[i].elem(0,1) = 0.0;
    virials[i].elem(0,2) = 0.0;
    virials[i].elem(1,0) = 0.0;
    virials[i].elem(1,1) = 0.0;
    virials[i].elem(1,2) = 0.0;
    virials[i].elem(2,0) = 0.0;
    virials[i].elem(2,1) = 0.0;
    virials[i].elem(2,2) = 0.0;

    s1 = matter[i];
    t1 = elem.name2idx( matter[i] );

    if (t1 != t1_old) sys_single_elem=false;

    for (j=0; j<neighborcollection[i].size(); ++j){
      k  = neighborcollection[i][j];
      s2 = matter[k];
      t2 = elem.name2idx( matter[k] );

      //if (t1 != t2) sys_single_elem=false;

      n_iacs++;
      if      ((*p_potinfo).basepot(t1,t2)=="ABOP") n_ABOP++;
      else if ((*p_potinfo).basepot(t1,t2)=="EAM" ) n_EAM++;
      
      /*
      if      ((*p_potinfo).basepot(s1,s2)=="EAM")  use_EAM =true;
      else if ((*p_potinfo).basepot(s1,s2)=="ABOP") use_ABOP=true;
      */
    }
    t1_old=t1;
  }

  //cout << "get_forces_and_energies_common(): Made it here 01" << endl;

  if (n_ABOP> 0 && n_ABOP==n_iacs) iac_pure_ABOP=true;
  if (n_EAM > 0 && n_EAM ==n_iacs) iac_pure_EAM =true;


  Ep_tot = 0.0;
  if (n_EAM>0)  Ep_tot += force_EAM();
  //cout << "get_forces_and_energies_common(): Made it here 02" << endl;
  if (n_ABOP>0) Ep_tot += force_ABOP();
  //cout << "get_forces_and_energies_common(): Made it here 03" << endl;


  return Ep_tot;
}



