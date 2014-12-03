


#ifndef MDSYS_HPP
#define MDSYS_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"

#include "atomsystem.hpp"
#include "mdsettings.hpp"
#include "potinfo.hpp"
//#include "compound.hpp"
#include "specs-fit-prop-pot.hpp"
#include "physconst.hpp"
#include "constants.hpp"
#include "exiterrors.hpp"

using namespace std;
using namespace utils;
using namespace exiterrors;
using namespace physconst;
using namespace constants;



class CompoundStructure;



class MDSystem
  :
  public AtomSystem
{
public:
  // System-related properties:
  string name;
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
  double V;
  double P, Px, Py, Pz;
  Matrix<double> stresstensor_xyz;
  Matrix<double> stresstensor_abc;
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
  Vector< Vector<double> > vel;
  Vector< Vector<double> > acc;
  Vector< Vector<double> > frc;
  Vector< Vector<double> > frc_num;
  Vector< Matrix<double> > virials;
  Vector<double> Ep;
  Vector<double> Ek;
  Vector< Vector<double> > dpos;
  Vector< Vector<double> > pos_int_tmp;
  Vector< Vector<double> > pos_int_ini;
  Vector< Vector<double> > pos_int_fin;
 

public:
  MDSystem();





  void create_from_structure(CompoundStructure & cmp,
			     double distmin);

  void handle_pbc_of_positions(const double lowlim=-1);

  void calc_closepacked_volume();

  void calc_forces_and_energies();

  double calc_potential_energy();

  double get_forces_and_energies_common();

  // **********************************************
  // List of force routines!
  // **********************************************
  double force_EAM();
  double force_ABOP();



  // MD:
  void relax(void);

  void transform_cell(const Matrix<double> & alpha_cart,
		      const double lowlim=-1);

  void get_virials(int nat, Matrix<double> & W,
		   double mux=1.0, double muy=1.0, double muz=1.0);

} ;









#endif


