


#ifndef MDSYS_HPP
#define MDSYS_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-errors.hpp"

#include "atomsystem.hpp"


#include "mdsettings.hpp"
#include "potinfo.hpp"
#include "compound.hpp"
//#include "specs-fit-prop-pot.hpp"
#include "physconst.hpp"
#include "constants.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;


using namespace utils;
using namespace physconst;
using namespace constants;



//class CompoundStructure;



class MDSystem
  :
  public AtomSystem
{
public:
  // System-related properties:
  std::string name;
  MDSettings_common specs_common;
  MDSettings        specs;

  Elements             elem;
  Interactions         iacs;
  PotentialInformation * p_potinfo;

  bool debug_creation;
  bool debug_mds;

  int N[3];

  double rcut_max;
  double dt;
  double T;
  double T_at_quench_start;
  double P, Px, Py, Pz;
  MatrixSq3<double> stresstensor_xyz;
  MatrixSq3<double> stresstensor_abc;
  double Ep_tot;
  double Ek_tot;
  double P_max, F_max, displ_max;
  bool get_pot_force;
  bool get_pot_energy;

  bool sys_single_elem;
  // **********************************************
  // Optimization help variables, depend on potential database
  // **********************************************
  bool iac_pure_ABOP;
  bool iac_pure_EAM;




  // Atom-related properties:
  Vector< Vector3<double> > vel;
  Vector< Vector3<double> > acc;
  Vector< Vector3<double> > frc;
  Vector< Vector3<double> > frc_num;
  Vector< MatrixSq3<double> > virials;
  Vector<double> Ep;
  Vector<double> Ek;
  Vector< Vector3<double> > dpos;
  Vector< Vector3<double> > pos_int_tmp;
  Vector< Vector3<double> > pos_int_ini;
  Vector< Vector3<double> > pos_int_fin;
 

public:
  MDSystem();



  void create_from_structure(CompoundStructure & cmp,
			     double distmin);

  void transform_cell(const MatrixSq3<double> & alpha_cart,
		      const double lowlim=-1);

  void calc_forces_and_energies();

  double calc_potential_energy();

  double get_forces_and_energies_common();

  // **********************************************
  // List of force routines!
  // **********************************************
  double force_EAM();

  double force_ABOP();
  void   force_ABOP_perriot_K(int i, int j, double & Kij, Vector3<double> dposij,
			      Matrix<double> & rcut_all,
			      Matrix<CutoffScreeningPair> & rcs_all,
			      Matrix<std::string> & basepot_all,
			      Matrix<int> & basepot_vecidx_all);
  void   force_ABOP_perriot_K_frc(int i, int j, double & Kij, Vector3<double> dposij, double pref,
				  Matrix<double> & rcut_all,
				  Matrix<CutoffScreeningPair> & rcs_all,
				  Matrix<std::string> & basepot_all,
				  Matrix<int> & basepot_vecidx_all);
  
  void relax(void);

  void get_virials(int nat, MatrixSq3<double> & W,
		   double mux=1.0, double muy=1.0, double muz=1.0);

} ;









#endif


