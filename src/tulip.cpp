


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include <boost/format.hpp>

#include <cmath>

#include "chisq-basics.hpp"
#include "exiterrors.hpp"
#include "funcfit-basics.hpp"
#include "funcfit-conjgrad.hpp"
#include "funcfit-exceptions.hpp"
#include "funcfit-ls-gauss-newton.hpp"
#include "funcfit-ls-leve-marq.hpp"
#include "funcfit-ls-powelldogleg.hpp"
#include "funcfit-powell.hpp"
#include "funcfit-simplexfit.hpp"
#include "funcfit-diffevol.hpp"
#include "funcfit-partswarm.hpp"
#include "funcfit-beecolony.hpp"
#include "funcfit-gravsearch.hpp"
#include "funcfit-simann.hpp"
#include "constants.hpp"
#include "nr-f1dim.hpp"
#include "nr-golden.hpp"
#include "nr-linemethod.hpp"
#include "param.hpp"
#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-matrix3.hpp"
#include "utils-matrix.hpp"
#include "utils-string.hpp"
#include "utils-vector.hpp"

#include "atomsystem.hpp"
#include "compound.hpp"
#include "elem-iacs.hpp"
#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mtwister.hpp"
#include "param-pot.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "potinfo.hpp"
#include "specs-fit-prop-pot.hpp"
#include "latcalc.hpp"
#include "report.hpp"


#include <omp.h>
#include "omp-basics.hpp"


#include <ctime>


using namespace std;
using namespace utils;
using namespace funcfit;
using namespace exiterrors;
using namespace constants;
using namespace physconst;
using boost::format;




#define VERSION 2.0


int main(int argc, char *argv[]){
  // double eps = numeric_limits<double>::epsilon();
  int i,j,k,p,iref, ivec;
  bool run_quick, debug_forces, debug_pressure, report_mds_steps;
  Vector<bool> debug_fit_prop(5, false),debug_fit_pot(5,false);
  string arg, potfile, geomfile, specsfile;
  bool potfileOK, geomfileOK, specsfileOK;
  istringstream sstream;
  string tulip_info, dformat("%15.10f");

  OMP_Info omp_info;
  int omp_tid;
  // Default:
  int omp_nt_try = 1;

  double ts1 = omp_get_wtime();
  //clock_t clock1 = std::clock();


  tulip_info = "TULIP version " + tostring(VERSION) + " (c) Krister Henriksson 2006-2014";

  if (argc < 6){
    cout << tulip_info << endl;
    cout << "Purpose: Fit data to an interatomic potential." << endl;
    cout << "Usage:" << endl;
    cout << "     " << string(argv[0]) << " arguments [options]" << endl;
    cout << "Arguments:" << endl;
    cout << "     -pf file           Path to file containing potential information." << endl;
    cout << "     -gf file           Path to file containing geometry information." << endl;
    cout << "     -sf file           Path to file containing technical specifications about the calculations." << endl;
    cout << "" << endl;

    cout << "Options:" << endl;
    cout << "     -quick             Calculate properties of read-in geometries and quit. No relaxation" << endl;
    cout << "                        is performed on read-in compounds. Reference compounds are relaxed." << endl;
    cout << "     -dfitpropn         Show information about fitting of properties. Here 'n' must be" << endl;
    cout << "                        an integer. Supported: 0-4. 0: debug fitting method. 1-4: debug deeper." << endl;
    cout << "                        lying methods used by the fitting method. Default: not used" << endl;
    cout << "                        NOTE: 0 also shows some info about the initial Chi^2 object." << endl;
    cout << "     -dfitpotn          Show information about fitting of potentials. Here 'n' have a similar" << endl;
    cout << "                        role as for fitting of the properties." << endl;
    cout << "                        NOTE 1: 0 also shows some info about the initial Chi^2 object." << endl;
    cout << "                        NOTE 2: 'fitpot0' is always set to true, others are false by default." << endl;
    cout << "" << endl;

    cout << "     -dforces           Debug the forces. Default: not used" << endl;
    cout << "     -dpressure         Debug the pressure. Default: not used" << endl;
    cout << "     -dmdsprop          Debug MDS runs of the structures. Default: not used" << endl;

    cout << "     -dall              Activate all debugging options (top level only). Default: not used" << endl;
    
    cout << "" << endl;
    cout << "     -omp               Request maximal number of threads ("
	 << omp_info.nt_max() << ") for any OpenMP parts." << endl;
    cout << "     -omp_nt num        Request 'num' number of threads for any OpenMP parts. Default: "
	 << omp_info.nt_use() << endl;
    cout << "" << endl;

    return 0;
  }


  potfileOK = geomfileOK = specsfileOK = false;
  run_quick = false;
  debug_forces = false;
  debug_pressure = false;
  report_mds_steps = false;
  debug_fit_pot[0]  = true;


  /* ###############################################################################
     Parse the options.
     ############################################################################### */
  for (i=1; i<argc; i++){

    if (string(argv[i])=="-pf"){
      arg = string(argv[i+1]); sstream.str(arg);
      sstream >> potfile;
      sstream.clear();
      potfileOK = true; i++;
    }
    else if (string(argv[i])=="-gf"){
      arg = string(argv[i+1]); sstream.str(arg);
      sstream >> geomfile;
      sstream.clear();
      geomfileOK = true; i++;
    }
    else if (string(argv[i])=="-sf"){
      arg = string(argv[i+1]); sstream.str(arg);
      sstream >> specsfile;
      sstream.clear();
      specsfileOK = true; i++;
    }
    else if (string(argv[i])=="-quick"){
      run_quick = true;
    }
    else if (string(argv[i])=="-dfitprop0"){ debug_fit_prop[0] = true; }
    else if (string(argv[i])=="-dfitprop1"){ debug_fit_prop[1] = true; }
    else if (string(argv[i])=="-dfitprop2"){ debug_fit_prop[2] = true; }
    else if (string(argv[i])=="-dfitprop3"){ debug_fit_prop[3] = true; }
    else if (string(argv[i])=="-dfitprop4"){ debug_fit_prop[4] = true; }

    else if (string(argv[i])=="-dfitpot0"){ debug_fit_pot[0] = true; }
    else if (string(argv[i])=="-dfitpot1"){ debug_fit_pot[1] = true; }
    else if (string(argv[i])=="-dfitpot2"){ debug_fit_pot[2] = true; }
    else if (string(argv[i])=="-dfitpot3"){ debug_fit_pot[3] = true; }
    else if (string(argv[i])=="-dfitpot4"){ debug_fit_pot[4] = true; }

    else if (string(argv[i])=="-dforces"){
      debug_forces = true;
    }
    else if (string(argv[i])=="-dpressure"){
      debug_pressure = true;
    }

    else if (string(argv[i])=="-dmdsprop"){
      report_mds_steps = true;
    }
    else if (string(argv[i])=="-dall"){
      debug_fit_prop[0] = true;
      debug_fit_prop[1] = true;
      debug_fit_prop[2] = true;
      debug_fit_prop[3] = true;
      debug_fit_prop[4] = true;
      debug_fit_pot[0]  = true;
      debug_fit_pot[1]  = true;
      debug_fit_pot[2]  = true;
      debug_fit_pot[3]  = true;
      debug_fit_pot[4]  = true;
      debug_forces = true;
      debug_pressure = true;
      report_mds_steps = true;
    }
    else if (string(argv[i])=="-omp"){
      omp_nt_try = omp_info.nt_max();
    }
    else if (string(argv[i])=="-omp_nt"){
      arg = string(argv[i+1]); sstream.str(arg);
      sstream >> omp_nt_try;
      sstream.clear();
      i++;
    }
    
  }



  
  if (! potfileOK) aborterror("Error: Potential information file not specified. Exiting.");
  if (! geomfileOK) aborterror("Error: Geometry information file not specified. Exiting.");
  if (! specsfileOK) aborterror("Error: Specifications file not specified. Exiting.");

  cout << "Potential info file     : " << potfile << endl;
  cout << "Geometry info file      : " << geomfile << endl;
  cout << "Specifications info file: " << specsfile << endl;


  cout << "|||||||||||||||||||||||||||||||||||||||||||||||" << endl;

  omp_info.nt_use( omp_nt_try );
#pragma omp parallel private(omp_tid)
  {
    omp_tid = omp_info.tid();
    if (omp_tid==0){
      cout << "OpenMP threads used     : " << omp_info.nt_use() << endl;
      cout << "OpenMP threads max count: " << omp_info.nt_max() << endl;
    }
  }

  cout << "|||||||||||||||||||||||||||||||||||||||||||||||" << endl;




  cout << "**************************************************************" << endl;
  cout << "**************************************************************" << endl;
  cout << "" << endl;
  cout << tulip_info << endl;
  cout << "" << endl;
  cout << "**************************************************************" << endl;
  cout << "**************************************************************" << endl;

  cout << "" << endl;
  cout << "STARTING UP" << endl;
  cout << "  1. Constructing potential information object using potential and specifications files ..." << endl;

  PotentialInformationFit potinfo(potfile, specsfile);
  Elements elem = potinfo.elem;

  cout << "  2. Constructing list of compounds to use for fitting, using geometry information file ..." << endl;
  CompoundListFit complistfit(elem, potinfo.specs_prop.mds_specs, geomfile);
  int ncomp = complistfit.compounds.size();

  cout << "  3. Constructing parameter object ..." << endl;
  ParamPot param( &potinfo );

  cout << "" << endl;
  cout << "STARTUP COMPLETE" << endl;
  cout << "" << endl;


  potinfo.omp_info = omp_info;


  potinfo.specs_prop.mds_specs_common.debug_forces = debug_forces;
  potinfo.specs_prop.mds_specs_common.debug_pressure = debug_pressure;
  potinfo.specs_prop.mds_specs_common.report_step  = report_mds_steps;
  potinfo.specs_prop.mds_specs_common.quick_mode   = run_quick;



  potinfo.specs_prop.debug_fit_level0 = debug_fit_prop[0];
  potinfo.specs_prop.debug_fit_level1 = debug_fit_prop[1];
  potinfo.specs_prop.debug_fit_level2 = debug_fit_prop[2];
  potinfo.specs_prop.debug_fit_level3 = debug_fit_prop[3];
  potinfo.specs_prop.debug_fit_level4 = debug_fit_prop[4];

  potinfo.specs_pot.debug_fit_level0  = debug_fit_pot[0];
  potinfo.specs_pot.debug_fit_level1  = debug_fit_pot[1];
  potinfo.specs_pot.debug_fit_level2  = debug_fit_pot[2];
  potinfo.specs_pot.debug_fit_level3  = debug_fit_pot[3];
  potinfo.specs_pot.debug_fit_level4  = debug_fit_pot[4];


  // If true, then the quick_mode option must be set for all compounds
  // individually. For the reference compounds quick mode will not be used.



  // ################################################################################
  // ################################################################################
  //
  // Report on all settings ............
  //
  // ################################################################################
  // ################################################################################

  int nref;
  string s1, s2, s3;

  // ################################################################################  
  // Elements
  // ################################################################################
  cout << "" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << "Elements:" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  nref = potinfo.elem.nelem();
  for (i=0; i<nref; ++i){
    s1 = potinfo.elem.idx2name(i);
    cout << s1 << " atomtype " << potinfo.elem.atomtype(s1)
	 << " mass (amu)  " << potinfo.elem.mass(s1)
	 << " reference lattice " << potinfo.elem.reflat(s1) << endl;
  }

  // ################################################################################  
  // Interactions
  // ################################################################################  
  cout << "" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << "Interactions:" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  nref = potinfo.elem.nelem();
  for (i=0; i<nref; ++i){
    for (j=i; j<nref; ++j){
      s1 = potinfo.elem.idx2name(i);
      s2 = potinfo.elem.idx2name(j);
      
      cout << s1 << "-" << s2 << " iac " << potinfo.iacs.name(s1,s2)
	   << "(" << potinfo.basepot(s1,s2) << ")"
	   << " is fittable? " << potinfo.is_fittable(s1,s2)
	   << " uses reppot? " << potinfo.use_reppot(s1,s2);

      if (potinfo.use_reppot(s1,s2)){
	ivec = potinfo.reppot_vecidx(s1,s2);
	cout << " bermi= "  << potinfo.pot_Reppot[ivec].bfermi
	     << " rfermi= " << potinfo.pot_Reppot[ivec].rfermi;
      }

      cout << endl;
    }
  }

  cout << "" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << "FIXED interactions which are parametrized:" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  report_pot( &potinfo, false, true );

  cout << "" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << "FITTABLE interactions (parametrized):" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  report_pot( &potinfo, true, false );




  // ################################################################################
  // Specifications for compounds
  // ################################################################################
  cout << "" << endl;
  cout << "##########################################################################" << endl;
  cout << "Specifications for calculation of properties of compounds:" << endl;
  cout << "##########################################################################" << endl;
  cout << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << "General" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << endl;
  cout << "Random number seed                                 : " << potinfo.specs_prop.seed << endl;
  cout << "Fitting method                                     : " << potinfo.specs_prop.fitmet << endl;
  cout << "Min. iterations                                    : " << potinfo.specs_prop.nitermin << endl;
  cout << "Max. iterations                                    : " << potinfo.specs_prop.nitermax << endl;
  cout << "Tolerance for convergence of ChiSq                 : " << potinfo.specs_prop.functolabs << endl;
  cout << "Tolerance for convergence of ChiSq changes         : " << potinfo.specs_prop.functolrel << endl;
  cout << "Tolerance for convergence of ChiSq gradient        : " << potinfo.specs_prop.gradtolabs << endl;
  cout << "Tolerance for convergence of step length           : " << potinfo.specs_prop.steptolabs << endl;
  cout << "Tolerance for convergence of step length changes   : " << potinfo.specs_prop.steptolrel << endl;
  cout << "DOG-LEG: Initial trust region radius               : " << potinfo.specs_prop.dogleg_radius << endl;
  cout << "DOG-LEG: Smallest allowed trust region radius      : " << potinfo.specs_prop.dogleg_minradius << endl;
  cout << "SIMPLEX: Displacement when creating initial simplex: " << potinfo.specs_prop.simplex_delta << endl;
  cout << "Debug: level0                                      : " << potinfo.specs_prop.debug_fit_level0 << endl;
  cout << "Debug: level1                                      : " << potinfo.specs_prop.debug_fit_level1 << endl;
  cout << "Lattice tolerance                                  : " << potinfo.specs_prop.lattol << endl;
  cout << endl;
  cout << "Bulk modulus: Minimum strain (e.g. -0.01)          : " << potinfo.specs_prop.BM_fmin << endl;
  cout << "Bulk modulus: Maximum strain (e.g.  0.01)          : " << potinfo.specs_prop.BM_fmax << endl;
  cout << "Bulk modulus: Number of strain points (e.g. 10)    : " << potinfo.specs_prop.BM_Nf << endl;
  cout << "Bulk modulus: Multiplicative uncertainty factor (e.g. 0.10): " << potinfo.specs_prop.BM_ef << endl;
  cout << endl;
  cout << "Elastic moduli: Minimum strain (e.g. -0.01)        : " << potinfo.specs_prop.C_fmin << endl;
  cout << "Elastic moduli: Maximum strain (e.g.  0.01)        : " << potinfo.specs_prop.C_fmax << endl;
  cout << "Elastic moduli: Number of strain points (e.g. 10)  : " << potinfo.specs_prop.C_Nf << endl;
  cout << "Elastic moduli: Multiplicative uncertainty factor (e.g. 0.10): " << potinfo.specs_prop.C_ef << endl;
  cout << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << "Settings for MD simulations of read-in compounds" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << endl;
  cout << "Skin thickness of neighbor list                         : " << potinfo.specs_prop.mds_specs.skint << endl;
  cout << "Seed for random numbers                                 : " << potinfo.specs_prop.mds_specs.seed << endl;
  cout << "Starting time for simulation (usually 0.0)              : " << potinfo.specs_prop.mds_specs.tstart << endl;
  cout << "Ending time for simulation                              : " << potinfo.specs_prop.mds_specs.tend << endl;
  cout << "Dump info every n:th step, n is                         : " << potinfo.specs_prop.mds_specs.ndump << endl;
  cout << "Starting/desired temperature (e.g. 300.0)               : " << potinfo.specs_prop.mds_specs.Tstart << endl;
  cout << "Initial time step                                       : " << potinfo.specs_prop.mds_specs.dt << endl;
  cout << "  Maximum time step                                     : " << potinfo.specs_prop.mds_specs.max_dt << endl;
  cout << "  Maximum allowed energy change for an atom/timestep    : " << potinfo.specs_prop.mds_specs.max_dE << endl;
  cout << "  Maximum allowed distance traveled for an atom/timestep: " << potinfo.specs_prop.mds_specs.max_dr << endl;
  cout << "Use temperature control (Berendsen)?                    : " << potinfo.specs_prop.mds_specs.use_Tcontrol << endl;
  cout << "  T control: time constant (fs)                         : " << potinfo.specs_prop.mds_specs.btc_tau << endl;
  cout << "  T control: desired temperature                        : " << potinfo.specs_prop.mds_specs.btc_T0 << endl;
  cout << "Use quench?                                             : " << potinfo.specs_prop.mds_specs.use_quench << endl;
  cout << "  Quench: start time (fs)                               : " << potinfo.specs_prop.mds_specs.quench_tstart << endl;
  cout << "  Quench: rate (K/fs) (always pos. value)               : " << potinfo.specs_prop.mds_specs.quench_rate << endl;
  cout << "Use pressure control (Berendsen)?                       : " << potinfo.specs_prop.mds_specs.use_Pcontrol << endl;
  cout << "  P control: time constant (fs)                         : " << potinfo.specs_prop.mds_specs.bpc_tau << endl;
  cout << "  P control: desired pressure (GPa)                     : " << potinfo.specs_prop.mds_specs.bpc_P0 << endl; 
  cout << "  P control: scale                                      : " << potinfo.specs_prop.mds_specs.bpc_scale << endl; 
  cout << "    The scale is usually equal to the approximate bulk modulus (GPa)." << endl;
  cout << "" << endl;
  cout << "Options: heating allowed?                               : " << potinfo.specs_prop.mds_specs.heating_allowed << endl; 
  cout << "Options: fixed geometry?                                : " << potinfo.specs_prop.mds_specs.fixed_geometry << endl; 
  cout << "Options: quench always?                                 : " << potinfo.specs_prop.mds_specs.quench_always << endl; 
  cout << "  If true, all velocities zeroed at every time step." << endl;
  cout << "" << endl;
  cout << "Overruling options: report step?                                   : " << potinfo.specs_prop.mds_specs_common.report_step << endl;
  cout << "  If true, writes out physical info of the system at every time step." << endl;
  cout << "Overruling options: debug forces?                                  : " << potinfo.specs_prop.mds_specs_common.debug_forces << endl; 
  cout << "Overruling options: debug pressure?                                : " << potinfo.specs_prop.mds_specs_common.debug_pressure << endl; 
  cout << "Overruling options: quick mode?                                    : " << potinfo.specs_prop.mds_specs_common.quick_mode << endl; 
  cout << "  If true, no MD relaxation, only evaluation of potential energy." << endl;
  cout << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << "Settings for MD simulations of reference compounds" << endl;
  cout << "--------------------------------------------------------------------------" << endl;
  cout << endl;
  cout << "Skin thickness of neighbor list                         : " << potinfo.specs_prop.mds_specs_ref.skint << endl;
  cout << "Seed for random numbers                                 : " << potinfo.specs_prop.mds_specs_ref.seed << endl;
  cout << "Starting time for simulation (usually 0.0)              : " << potinfo.specs_prop.mds_specs_ref.tstart << endl;
  cout << "Ending time for simulation                              : " << potinfo.specs_prop.mds_specs_ref.tend << endl;
  cout << "Dump info every n:th step, n is                         : " << potinfo.specs_prop.mds_specs_ref.ndump << endl;
  cout << "Starting/desired temperature (e.g. 300.0)               : " << potinfo.specs_prop.mds_specs_ref.Tstart << endl;
  cout << "Initial time step                                       : " << potinfo.specs_prop.mds_specs_ref.dt << endl;
  cout << "  Maximum time step                                     : " << potinfo.specs_prop.mds_specs_ref.max_dt << endl;
  cout << "  Maximum allowed energy change for an atom/timestep    : " << potinfo.specs_prop.mds_specs_ref.max_dE << endl;
  cout << "  Maximum allowed distance traveled for an atom/timestep: " << potinfo.specs_prop.mds_specs_ref.max_dr << endl;
  cout << "Use temperature control (Berendsen)?                    : " << potinfo.specs_prop.mds_specs_ref.use_Tcontrol << endl;
  cout << "  T control: time constant (fs)                         : " << potinfo.specs_prop.mds_specs_ref.btc_tau << endl;
  cout << "  T control: desired temperature                        : " << potinfo.specs_prop.mds_specs_ref.btc_T0 << endl;
  cout << "Use quench?                                             : " << potinfo.specs_prop.mds_specs_ref.use_quench << endl;
  cout << "  Quench: start time (fs)                               : " << potinfo.specs_prop.mds_specs_ref.quench_tstart << endl;
  cout << "  Quench: rate (K/fs) (always pos. value)               : " << potinfo.specs_prop.mds_specs_ref.quench_rate << endl;
  cout << "Use pressure control (Berendsen)?                       : " << potinfo.specs_prop.mds_specs_ref.use_Pcontrol << endl;
  cout << "  P control: time constant (fs)                         : " << potinfo.specs_prop.mds_specs_ref.bpc_tau << endl;
  cout << "  P control: desired pressure (GPa)                     : " << potinfo.specs_prop.mds_specs_ref.bpc_P0 << endl; 
  cout << "  P control: scale                                      : " << potinfo.specs_prop.mds_specs_ref.bpc_scale << endl; 
  cout << "    The scale is usually equal to the approximate bulk modulus (GPa)." << endl;
  cout << "" << endl;
  cout << "Options: heating allowed?                               : " << potinfo.specs_prop.mds_specs_ref.heating_allowed << endl; 
  cout << "Options: fixed geometry?                                : " << potinfo.specs_prop.mds_specs_ref.fixed_geometry << endl; 
  cout << "Options: quench always?                                 : " << potinfo.specs_prop.mds_specs_ref.quench_always << endl; 
  cout << "  If true, all velocities zeroed at every time step." << endl;
  cout << "" << endl;
  cout << "Overruling options: report step?                                   : " << potinfo.specs_prop.mds_specs_common.report_step << endl;
  cout << "  If true, writes out physical info of the system at every time step." << endl;
  cout << "Overruling options: debug forces?                                  : " << potinfo.specs_prop.mds_specs_common.debug_forces << endl; 
  cout << "Overruling options: debug pressure?                                : " << potinfo.specs_prop.mds_specs_common.debug_pressure << endl; 
  cout << "Overruling options: quick mode?                                    : " << potinfo.specs_prop.mds_specs_common.quick_mode << endl; 
  cout << "  If true, no MD relaxation, only evaluation of potential energy." << endl;



  // ################################################################################  
  // Specifications for potentials
  // ################################################################################  
  cout << "" << endl;
  cout << "##########################################################################" << endl;
  cout << "Specifications for potential fitting:" << endl;
  cout << "##########################################################################" << endl;
  cout << "" << endl;
  cout << "Random number seed                                 : " << potinfo.specs_pot.seed << endl;
  cout << "Fitting method                                     : " << potinfo.specs_pot.fitmet << endl;
  cout << "Min. iterations                                    : " << potinfo.specs_pot.nitermin << endl;
  cout << "Max. iterations                                    : " << potinfo.specs_pot.nitermax << endl;
  cout << "Tolerance for convergence of ChiSq                 : " << potinfo.specs_pot.functolabs << endl;
  cout << "Tolerance for convergence of ChiSq changes         : " << potinfo.specs_pot.functolrel << endl;
  cout << "Tolerance for convergence of ChiSq gradient        : " << potinfo.specs_pot.gradtolabs << endl;
  cout << "Tolerance for convergence of step length           : " << potinfo.specs_pot.steptolabs << endl;
  cout << "Tolerance for convergence of step length changes   : " << potinfo.specs_pot.steptolrel << endl;
  cout << "DOG-LEG: Initial trust region radius               : " << potinfo.specs_pot.dogleg_radius << endl;
  cout << "DOG-LEG: Smallest allowed trust region radius      : " << potinfo.specs_pot.dogleg_minradius << endl;
  cout << "SIMPLEX: Displacement when creating initial simplex: " << potinfo.specs_pot.simplex_delta << endl;
  cout << "Debug: level0                                      : " << potinfo.specs_pot.debug_fit_level0 << endl;
  cout << "Debug: level1                                      : " << potinfo.specs_pot.debug_fit_level1 << endl;
  cout << "Penalty function: barrier                          : " << potinfo.specs_pot.barrier_scale << endl;











  // ################################################################################  
  // Compounds
  // ################################################################################  

  // Some adjustments ...
  for (i=0; i<ncomp; ++i){
    complistfit.compounds[i].check_and_fix_uses();
  }

  // Report:
  cout << "" << endl;
  cout << "##########################################################################" << endl;
  cout << "Compounds used for fitting:" << endl;
  cout << "##########################################################################" << endl;
  cout << "" << endl;
  
  for (i=0; i<ncomp; ++i){
    cout << "Compound " << i+1 << " of " << ncomp << endl;
    cout << "  Name                : " << complistfit.compounds[i].name << endl;
    cout << "  Read from file      : " << complistfit.compounds[i].filename << endl;
    cout << "  PBC                 : " << complistfit.compounds[i].pbc[0] << " "
	 << complistfit.compounds[i].pbc[1] << " "
	 << complistfit.compounds[i].pbc[2] << endl;
    cout << "  Crystal system      : " << complistfit.compounds[i].csystem << endl;
    cout << "  Scale factor        : " << complistfit.compounds[i].scalefactor << endl;
    cout << "  Use internal format?: " << complistfit.compounds[i].use_int << endl;
    cout << "  Direction vector 1  : " 
	 << format("%15.10f ") % complistfit.compounds[i].u1_vec[0]
	 << format("%15.10f ") % complistfit.compounds[i].u1_vec[1]
	 << format("%15.10f")  % complistfit.compounds[i].u1_vec[2] << endl;
    cout << "  Direction vector 2  : " 
	 << format("%15.10f ") % complistfit.compounds[i].u2_vec[0]
	 << format("%15.10f ") % complistfit.compounds[i].u2_vec[1]
	 << format("%15.10f")  % complistfit.compounds[i].u2_vec[2] << endl;
    cout << "  Direction vector 3  : " 
	 << format("%15.10f ") % complistfit.compounds[i].u3_vec[0]
	 << format("%15.10f ") % complistfit.compounds[i].u3_vec[1]
	 << format("%15.10f")  % complistfit.compounds[i].u3_vec[2] << endl;
    cout << "  Number of basis atoms: " << complistfit.compounds[i].nbasis << endl;


    string dumpfile(complistfit.compounds[i].filename + ".xyz");
    cout << "  Atom system will be written to file " << dumpfile << endl;
    ofstream fout;
    fout.open(dumpfile.c_str());
    fout << complistfit.compounds[i].basis_vecs.size() << endl;
    fout << "Frame number 0 time 0.0 fs boxsize "
	 << format("%15.10f ") % complistfit.compounds[i].u1_vec.magn()
	 << format("%15.10f ") % complistfit.compounds[i].u2_vec.magn()
	 << format("%15.10f ") % complistfit.compounds[i].u3_vec.magn()
	 << "   "
	 << format("%10.6f ") % complistfit.compounds[i].u1_vec[0]
	 << format("%10.6f ") % complistfit.compounds[i].u1_vec[1]
	 << format("%10.6f")  % complistfit.compounds[i].u1_vec[2]
	 << format("%10.6f ") % complistfit.compounds[i].u2_vec[0]
	 << format("%10.6f ") % complistfit.compounds[i].u2_vec[1]
	 << format("%10.6f")  % complistfit.compounds[i].u2_vec[2]
	 << format("%10.6f ") % complistfit.compounds[i].u3_vec[0]
	 << format("%10.6f ") % complistfit.compounds[i].u3_vec[1]
	 << format("%10.6f")  % complistfit.compounds[i].u3_vec[2] << endl;
    for (j=0; j<complistfit.compounds[i].basis_vecs.size(); ++j)
      fout << complistfit.compounds[i].basis_elems[j] << "  "
	   << format("%15.10f ") % complistfit.compounds[i].basis_vecs[j][0]
	   << format("%15.10f ") % complistfit.compounds[i].basis_vecs[j][1]
	   << format("%15.10f")  % complistfit.compounds[i].basis_vecs[j][2]
	   << j << endl;
    fout.close();
    fout.clear();


    cout << "  Compound-specific MD settings" << endl;
    cout << "  -----------------------------" << endl;
    cout << "Skin thickness of neighbor list                         : " << complistfit.compounds[i].mds_specs.skint << endl;
    cout << "Seed for random numbers                                 : " << complistfit.compounds[i].mds_specs.seed << endl;
    cout << "Starting time for simulation (usually 0.0)              : " << complistfit.compounds[i].mds_specs.tstart << endl;
    cout << "Ending time for simulation                              : " << complistfit.compounds[i].mds_specs.tend << endl;
    cout << "Dump info every n:th step, n is                         : " << complistfit.compounds[i].mds_specs.ndump << endl;
    cout << "Starting/desired temperature (e.g. 300.0)               : " << complistfit.compounds[i].mds_specs.Tstart << endl;
    cout << "Initial time step                                       : " << complistfit.compounds[i].mds_specs.dt << endl;
    cout << "  Maximum time step                                     : " << complistfit.compounds[i].mds_specs.max_dt << endl;
    cout << "  Maximum allowed energy change for an atom/timestep    : " << complistfit.compounds[i].mds_specs.max_dE << endl;
    cout << "  Maximum allowed distance traveled for an atom/timestep: " << complistfit.compounds[i].mds_specs.max_dr << endl;
    cout << "Use temperature control (Berendsen)?                    : " << complistfit.compounds[i].mds_specs.use_Tcontrol << endl;
    cout << "  T control: time constant (fs)                         : " << complistfit.compounds[i].mds_specs.btc_tau << endl;
    cout << "  T control: desired temperature                        : " << complistfit.compounds[i].mds_specs.btc_T0 << endl;
    cout << "Use quench?                                             : " << complistfit.compounds[i].mds_specs.use_quench << endl;
    cout << "  Quench: start time (fs)                               : " << complistfit.compounds[i].mds_specs.quench_tstart << endl;
    cout << "  Quench: rate (K/fs) (always pos. value)               : " << complistfit.compounds[i].mds_specs.quench_rate << endl;
    cout << "Use pressure control (Berendsen)?                       : " << complistfit.compounds[i].mds_specs.use_Pcontrol << endl;
    cout << "  P control: time constant (fs)                         : " << complistfit.compounds[i].mds_specs.bpc_tau << endl;
    cout << "  P control: desired pressure (GPa)                     : " << complistfit.compounds[i].mds_specs.bpc_P0 << endl; 
    cout << "  P control: scale                                      : " << complistfit.compounds[i].mds_specs.bpc_scale << endl; 
    cout << "    The scale is usually equal to the approximate bulk modulus (GPa)." << endl;
    cout << "" << endl;
    cout << "Options: heating allowed?                               : " << complistfit.compounds[i].mds_specs.heating_allowed << endl; 
    cout << "Options: fixed geometry?                                : " << complistfit.compounds[i].mds_specs.fixed_geometry << endl; 
    cout << "Options: quench always?                                 : " << complistfit.compounds[i].mds_specs.quench_always << endl; 
    cout << "  If true, all velocities zeroed at every time step." << endl;

    cout << "--------------------------------------------------------------------------" << endl;
  }




  cout << "INFO: Number of free fitting prameters: " << param.NXfree() << endl;
  cout << "INFO: Number of data points           : " << complistfit.NData() << endl;
  cout << "--------------------------------------------------------------------------" << endl;


  Vector<CompoundStructureFit> DX;
  Vector<double> DY, DUY, DWY;



  // DataX:
  DX.resize(0);
  for (i=0; i<ncomp; ++i){

    DX.push_back(complistfit.compounds[i]);
  }


  // DataY: 
  // WeightsDataY: 

  DY.resize(0);
  DUY.resize(0);
  DWY.resize(0);

  for (i=0; i<ncomp; ++i){
    
    if (complistfit.compounds[i].prop_use.a){ DY.push_back(complistfit.compounds[i].prop_readin.a);
      if (complistfit.compounds[i].use_u.a){
	DUY.push_back(complistfit.compounds[i].prop_u.a);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.a);
      }
    }
    if (complistfit.compounds[i].prop_use.b){ DY.push_back(complistfit.compounds[i].prop_readin.b);
      if (complistfit.compounds[i].use_u.b){
	DUY.push_back(complistfit.compounds[i].prop_u.b);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.b);
      }
    }

    if (complistfit.compounds[i].prop_use.c){ DY.push_back(complistfit.compounds[i].prop_readin.c);
      if (complistfit.compounds[i].use_u.c){
	DUY.push_back(complistfit.compounds[i].prop_u.c);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.c);
      }
    }

    if (complistfit.compounds[i].prop_use.bpa){ DY.push_back(complistfit.compounds[i].prop_readin.bpa);
      if (complistfit.compounds[i].use_u.bpa){
	DUY.push_back(complistfit.compounds[i].prop_u.bpa);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.bpa);
      }
    }

    if (complistfit.compounds[i].prop_use.cpa){ DY.push_back(complistfit.compounds[i].prop_readin.cpa);
      if (complistfit.compounds[i].use_u.cpa){
	DUY.push_back(complistfit.compounds[i].prop_u.cpa);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.cpa);
      }
    }

    if (complistfit.compounds[i].prop_use.r0){ DY.push_back(complistfit.compounds[i].prop_readin.r0);
      if (complistfit.compounds[i].use_u.r0){
	DUY.push_back(complistfit.compounds[i].prop_u.r0);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.r0);
      }
    }

    if (complistfit.compounds[i].prop_use.angle_ab){ DY.push_back(complistfit.compounds[i].prop_readin.angle_ab);
      if (complistfit.compounds[i].use_u.angle_ab){
	DUY.push_back(complistfit.compounds[i].prop_u.angle_ab);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.angle_ab);
      }
    }

    if (complistfit.compounds[i].prop_use.angle_ac){ DY.push_back(complistfit.compounds[i].prop_readin.angle_ac);
      if (complistfit.compounds[i].use_u.angle_ac){
	DUY.push_back(complistfit.compounds[i].prop_u.angle_ac);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.angle_ac);
      }
    }

    if (complistfit.compounds[i].prop_use.angle_bc){ DY.push_back(complistfit.compounds[i].prop_readin.angle_bc);
      if (complistfit.compounds[i].use_u.angle_bc){
	DUY.push_back(complistfit.compounds[i].prop_u.angle_bc);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.angle_bc);
      }
    }

    if (complistfit.compounds[i].prop_use.Vatom){ DY.push_back(complistfit.compounds[i].prop_readin.Vatom);
      if (complistfit.compounds[i].use_u.Vatom){
	DUY.push_back(complistfit.compounds[i].prop_u.Vatom);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.Vatom);
      }
    }

    if (complistfit.compounds[i].prop_use.Ecoh){ DY.push_back(complistfit.compounds[i].prop_readin.Ecoh);
      if (complistfit.compounds[i].use_u.Ecoh){
	DUY.push_back(complistfit.compounds[i].prop_u.Ecoh);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.Ecoh);
      }
    }

    if (complistfit.compounds[i].prop_use.Emix){ DY.push_back(complistfit.compounds[i].prop_readin.Emix);
      if (complistfit.compounds[i].use_u.Emix){
	DUY.push_back(complistfit.compounds[i].prop_u.Emix);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.Emix);
      }
    }

    if (complistfit.compounds[i].prop_use.B){ DY.push_back(complistfit.compounds[i].prop_readin.B);
      if (complistfit.compounds[i].use_u.B){
	DUY.push_back(complistfit.compounds[i].prop_u.B);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.B);
      }
    }

    if (complistfit.compounds[i].prop_use.Bp){ DY.push_back(complistfit.compounds[i].prop_readin.Bp);
      if (complistfit.compounds[i].use_u.Bp){
	DUY.push_back(complistfit.compounds[i].prop_u.Bp);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.Bp);
      }
    }

    for (k=0; k<6; ++k){
      for (p=0; p<6; ++p){
	if (complistfit.compounds[i].prop_use.C.elem(k,p)){ DY.push_back(complistfit.compounds[i].prop_readin.C.elem(k,p));
	  if (complistfit.compounds[i].use_u.C.elem(k,p)){
	    DUY.push_back(complistfit.compounds[i].prop_u.C.elem(k,p));
	    DWY.push_back(-1.0);
	  }
	  else {
	    DUY.push_back(-1.0);
	    DWY.push_back(complistfit.compounds[i].prop_w.C.elem(k,p));
	  }
	}
      }
    }

    if (complistfit.compounds[i].prop_use.Fmax){ DY.push_back(complistfit.compounds[i].prop_readin.Fmax);
      if (complistfit.compounds[i].use_u.Fmax){
	DUY.push_back(complistfit.compounds[i].prop_u.Fmax);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.Fmax);
      }
    }

    if (complistfit.compounds[i].prop_use.Pmax){ DY.push_back(complistfit.compounds[i].prop_readin.Pmax);
      if (complistfit.compounds[i].use_u.Pmax){
	DUY.push_back(complistfit.compounds[i].prop_u.Pmax);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.Pmax);
      }
    }

    if (complistfit.compounds[i].prop_use.displmax){ DY.push_back(complistfit.compounds[i].prop_readin.displmax);
      if (complistfit.compounds[i].use_u.displmax){
	DUY.push_back(complistfit.compounds[i].prop_u.displmax);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.displmax);
      }
    }


  }


  report_prop( complistfit.compounds, true );








  // ################################################################################
  // Get properties (at least cohesive energies) of reference structures !!!
  // ################################################################################

  cout << "" << endl;
  cout << "##########################################################################" << endl;
  cout << "Reference compounds" << endl;
  cout << "##########################################################################" << endl;
  cout << "" << endl;

  nref = potinfo.elem.nelem();

  potinfo.Ecoh_ref.resize(nref);
  for (iref=0; iref<nref; ++iref) potinfo.Ecoh_ref[iref] = 0;



  for (iref=0; iref<nref; ++iref){
    string sref = potinfo.elem.idx2name(iref);
    string latref = potinfo.elem.reflat(sref);

    if (latref=="none"){
      cout << "Elemental combination " << sref << "-" << sref << " is to be fitted."
	   << " No reference lattice calculation possible. Continuing." << endl;
      continue;
    }


    Vector<CompoundStructureFit> cmpref(1, CompoundStructureFit());;
    double a,b,c;

    // ###############################################################
    // 1. Create from generic model:
    // ###############################################################
    cmpref[0].create_from_model(latref, sref, sref);

    // ###############################################################
    // 2. Use specific compound data:
    // ###############################################################
    a = potinfo.elem.reflat_a(sref);
    b = potinfo.elem.reflat_b(sref);
    c = potinfo.elem.reflat_c(sref);
    cmpref[0].finalize(a,b,c);




    cout << "Reference compound " << iref+1 << " of " << nref << endl;
    cout << "  Name                : " << cmpref[0].name << endl;
    cout << "  PBC                 : " << cmpref[0].pbc[0] << " "
	 << cmpref[0].pbc[1] << " "
	 << cmpref[0].pbc[2] << endl;
    cout << "  Crystal system      : " << cmpref[0].csystem << endl;
    cout << "  Scale factor        : " << cmpref[0].scalefactor << endl;
    cout << "  Lattice parameters provided" << endl;
    cout << "                     a: " << a << endl;
    cout << "                     b: " << b << endl;
    cout << "                     c: " << c << endl;
    cout << "                => c/a: " << c/a << endl;
    cout << "  Use internal format?: " << cmpref[0].use_int << endl;
    cout << "  Direction vector 1  : " 
	 << format("%15.10f ") % cmpref[0].u1_vec[0]
	 << format("%15.10f ") % cmpref[0].u1_vec[1]
	 << format("%15.10f")  % cmpref[0].u1_vec[2] << endl;
    cout << "  Direction vector 2  : " 
	 << format("%15.10f ") % cmpref[0].u2_vec[0]
	 << format("%15.10f ") % cmpref[0].u2_vec[1]
	 << format("%15.10f")  % cmpref[0].u2_vec[2] << endl;
    cout << "  Direction vector 3  : " 
	 << format("%15.10f ") % cmpref[0].u3_vec[0]
	 << format("%15.10f ") % cmpref[0].u3_vec[1]
	 << format("%15.10f")  % cmpref[0].u3_vec[2] << endl;
    cout << "  Number of basis atoms: " << cmpref[0].nbasis << endl;






    // Set read-in properties and indicate which ones we want printed later:
    cmpref[0].prop_readin.a = a;
    cmpref[0].prop_use.a = true;
    cmpref[0].prop_readin.b = b;
    cmpref[0].prop_use.b = true;
    cmpref[0].prop_readin.c = c;
    cmpref[0].prop_use.c = true;

    cmpref[0].prop_use.Ecoh = true;

    cmpref[0].prop_use.B = true;
    cmpref[0].prop_use.Bp = true;

    bool anyC = true;
    for (i=0; i<6; ++i){
      for (j=0; j<6; ++j){
	cmpref[0].prop_use.C.elem(i,j) = true;
      }
    }

    cmpref[0].prop_use.Fmax = true;
    cmpref[0].prop_use.Pmax = true;
    cmpref[0].prop_use.displmax = true;

    if (!cmpref[0].pbc[0] || !cmpref[0].pbc[1] || !cmpref[0].pbc[2]){
      if (cmpref[0].nbasis==2)
	cmpref[0].prop_use.r0 = true;

      cmpref[0].prop_use.B  = false;
      cmpref[0].prop_use.Bp = false;
      cmpref[0].prop_use.Pmax = false;
      anyC = false;
      for (k=0; k<6; ++k)
	for (p=0; p<6; ++p)
	  cmpref[0].prop_use.C.elem(k,p) = false;
    }




    // MD settings for reference compounds:
    cmpref[0].mds_specs = potinfo.specs_prop.mds_specs_ref;





    // Calculate properties:
    latcalc(param, cmpref);

    // Report:
    cout << "Properties of reference compound for element " << sref << ":" << endl;
    cout << "  Lattice parameter a     : " << format("%15.10f") % cmpref[0].prop_pred.a << endl;
    cout << "  Lattice parameter b     : " << format("%15.10f") % cmpref[0].prop_pred.b << endl;
    cout << "  Lattice parameter c     : " << format("%15.10f") % cmpref[0].prop_pred.c << endl;
    cout << "  Cohesive energy   Ecoh  : " << format("%15.10f") % cmpref[0].prop_pred.Ecoh
	 << " for type " << potinfo.elem.name2idx(sref) << endl;

    potinfo.Ecoh_ref[ potinfo.elem.name2idx(sref) ] = cmpref[0].prop_pred.Ecoh;

    if (cmpref[0].prop_use.r0)
      cout << "  Dimer bond lenght r0    : " << format("%15.10f") % cmpref[0].prop_pred.r0 << endl;

    if (cmpref[0].prop_use.B)
      cout << "  Bulk modulus      B     : " << format("%15.10f") % cmpref[0].prop_pred.B << endl;
    if (cmpref[0].prop_use.Bp)
      cout << "  Pressure derivative of B: " << format("%15.10f") % cmpref[0].prop_pred.Bp << endl;

    if (anyC){
      if (cmpref[0].csystem=="cubic"){
	cout << "  Elastic constant     C11: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,0) << endl;
	cout << "  Elastic constant     C12: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,1) << endl;
	cout << "  Elastic constant     C44: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(3,3) << endl;
      }
      else if (cmpref[0].csystem=="hexagonal"){
	cout << "  Elastic constant     C11: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,0) << endl;
	cout << "  Elastic constant     C12: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,1) << endl;
	cout << "  Elastic constant     C13: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,2) << endl;
	cout << "  Elastic constant     C33: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(2,2) << endl;
	cout << "  Elastic constant     C44: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(3,3) << endl;
      }
      else if (cmpref[0].csystem=="orthorombic"){
	cout << "  Elastic constant     C11: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,0) << endl;
	cout << "  Elastic constant     C12: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,1) << endl;
	cout << "  Elastic constant     C13: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,2) << endl;
	cout << "  Elastic constant     C22: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(1,1) << endl;
	cout << "  Elastic constant     C23: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(1,2) << endl;
	cout << "  Elastic constant     C33: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(2,2) << endl;
	cout << "  Elastic constant     C44: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(3,3) << endl;
	cout << "  Elastic constant     C55: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(4,4) << endl;
	cout << "  Elastic constant     C66: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(5,5) << endl;
      }
    }
    cout << "  Maximum force           : " << format("%15.10f") % cmpref[0].prop_pred.Fmax << endl;
    if (cmpref[0].prop_use.Pmax)
      cout << "  Maximum pressure        : " << format("%15.10f") % cmpref[0].prop_pred.Pmax << endl;
    cout << "  Maximum displacement    : " << format("%15.10f") % cmpref[0].prop_pred.displmax << endl;


    cout << "--------------------------------------------------------------------------" << endl;




  }




  exit(1);



  /*

  if (debug_forces){
    cout << "Warning: Debugging of forces done, performing quick and dirty exiting from main routine. Bye." << endl;
    return EXIT_SUCCESS;
  }
  */




  // ################################################################################
  // Potential fitting
  // ################################################################################




  // #############################################################
  // Form merit function object
  // #############################################################

  cout << "Setting up merit function ..." << endl;


  ChiSqFunc<ParamPot, CompoundStructureFit, double> cs;
  Vector<double> Xopt;
  Cond_Conv  cond_conv;
  Cond_Debug cond_debug;
  Cond_Print cond_print;


  cs.Param() = param;
  cs.DataX() = DX;
  cs.DataY() = DY;
  cs.DataUncertaintyY()  = DUY;
  cs.DataWeightY()       = DWY;
  cs.ModelFuncPointer()  = latcalc;
  cs.ReportFuncPointer() = report_pot_prop;

  cs.barrier_scale() = potinfo.specs_pot.barrier_scale;

  cs.finalize_setup();







  cond_conv.functolabs = potinfo.specs_pot.functolabs;
  cond_conv.functolrel = potinfo.specs_pot.functolrel;
  cond_conv.gradtolabs = potinfo.specs_pot.gradtolabs;
  cond_conv.steptolabs = potinfo.specs_pot.steptolabs;
  cond_conv.steptolrel = potinfo.specs_pot.steptolrel;
  cond_conv.nitermin   = potinfo.specs_pot.nitermin;
  cond_conv.nitermax   = potinfo.specs_pot.nitermax;

  cond_conv.report_conv = potinfo.specs_pot.report_conv;
  cond_conv.prefix_report_conv = "POTFIT conv: ";


  cond_debug.debug_fit_level0 = potinfo.specs_pot.debug_fit_level0;
  cond_debug.debug_fit_level1 = potinfo.specs_pot.debug_fit_level1;
  cond_debug.debug_fit_level2 = potinfo.specs_pot.debug_fit_level2;
  cond_debug.debug_fit_level3 = potinfo.specs_pot.debug_fit_level3;
  cond_debug.debug_fit_level4 = potinfo.specs_pot.debug_fit_level4;
  cond_debug.prefix_debug_fit_level0 = "POTFIT debug0: ";
  cond_debug.prefix_debug_fit_level1 = "POTFIT debug1: ";
  cond_debug.prefix_debug_fit_level2 = "POTFIT debug2: ";
  cond_debug.prefix_debug_fit_level3 = "POTFIT debug3: ";
  cond_debug.prefix_debug_fit_level4 = "POTFIT debug4: ";


  cond_print.report_iter = true;
  cond_print.report_error= true;
  cond_print.report_warn = true;
  cond_print.prefix_report_iter  = "POTFIT iter: ";
  cond_print.prefix_report_warn  = "POTFIT warn: ";
  cond_print.prefix_report_error = "POTFIT error: ";



  if (cond_debug.debug_fit_level0){
    cond_conv.report_conv  = true;
    cond_print.report_iter = true;
    cond_print.report_warn = true;
    cond_print.report_error= true;
    cs.debug();
  }


  if (run_quick){
    cond_conv.nitermin   = 0;
    cond_conv.nitermax   = 0;
  }





  // Vector<double> latcalc(ParamPot & parpot, Vector<CompoundStructureFit> & cmpvec);

  cout << endl;
  cout << "******************************************************" << endl;
  cout << "******************************************************" << endl;
  cout << "**                                                  **" << endl;
  cout << "**          Starting potential fitting ...          **" << endl;
  cout << "**                                                  **" << endl;
  cout << "******************************************************" << endl;
  cout << "******************************************************" << endl;
  cout << endl;
  cout << "INFO: Number of free fitting prameters: " << param.NXfree() << endl;
  cout << "INFO: Number of data points           : " << complistfit.NData() << endl;
  cout << endl;


  int seed = potinfo.specs_pot.seed;

  /* ###############################################################################
     Fit function.
     ############################################################################### */
  if (potinfo.specs_pot.fitmet=="CG"){
    // Conjugate Gradients
    cout << "Using conjugate gradients method." << endl;
    ConjGrad< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="PM"){
    // Powell's method
    cout << "Using Powell's method." << endl;
    Powell< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="GN"){
    // Gauss-Newton
    cout << "Using Gauss-Newton method." << endl;
    GaussNewton< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="LM"){
    // Levenberg-Marquardt
    cout << "Using Levenberg-Marquardt method." << endl;
    LeveMarq< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="DL"){
    // Powell dog-leg
    cout << "Using Powell dog-leg method." << endl;
    PowellDogLeg< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       potinfo.specs_pot.dogleg_radius,
		       potinfo.specs_pot.dogleg_minradius,
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="SM"){
    // Simplex method
    cout << "Using Simplex method." << endl;
    SimplexFit< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    /*
    Vector<double> X_displ( cs.Param().X().size(), 
			    potinfo.specs_pot.simplex_delta );
    */
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
  }

  else if (potinfo.specs_pot.fitmet=="DE"){
    // Differential evolution
    cs.barrier_scale() = 0.0;
    DiffEvol< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="PS"){
    // Particle Swarm
    cs.barrier_scale() = 0.0;
    PartSwarm< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="BC"){
    // Bee colony
    cs.barrier_scale() = 0.0;
    BeeColony< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="GS"){
    // Gravitational Search
    cs.barrier_scale() = 0.0;
    GravSearch< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
  }
  else if (potinfo.specs_pot.fitmet=="SA"){
    // Simulated Annealing
    cs.barrier_scale() = 0.0;
    SimAnn< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);

    int N = cs.Param().X().size();
    Vector<double> Xd(N);
    double td;
    for (int i=0; i<N; ++i){
      td = cs.Param().X(i); if (td<0) td *= -1.0;
      Xd[i] = td;
    }

    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       Xd,
		       seed,
		       cond_conv, cond_debug, cond_print);
  }
  else
    aborterror("Error: Unknown fitting method " + potinfo.specs_pot.fitmet + ". " +
	       "Exiting.");






  /* ###############################################################################
     Report on fitted parameters.
     ############################################################################### */

  //cs.Param().X() = Xopt;
  //cout << "Made it here" << endl;
  //cout << "Xopt - cs.Param().X() is " << endl;
  //cout << Xopt - cs.Param().X() << endl;

  cout << endl;
  cout << "******************************************************" << endl;
  cout << "******************************************************" << endl;
  cout << "**                                                  **" << endl;
  cout << "**          Potential fitting terminated.           **" << endl;
  cout << "**                                                  **" << endl;
  cout << "******************************************************" << endl;
  cout << "******************************************************" << endl;
  cout << endl;


  // Calculate properties using optimized parameters:
  cs(Xopt);
  // Report:
  report_pot_prop( cs.Param(), cs.DataX(), cs.DataY(), cs.ModelDataY() );


  double ts2 = omp_get_wtime();
  //clock_t clock2 = std::clock();
  
  cout << "**********************************************************************" << endl;
  cout << "Total time (seconds) spent in program: "
       << ts2 - ts1
       << endl;
  cout << "**********************************************************************" << endl;




  return EXIT_SUCCESS;
}





#if 0
  Matrix<double> H;
  H = cs.Hessian(Xopt);
  int NXf = cs.Param().NXfree();
  Vector<double> Xfopt(cs.free_parameters(Xopt)), dXfopt(NXf);
  double re;

  for (i=0; i!=NXf; ++i)
    dXfopt[i] = sqrt(2.0 / H.elem(i,i));

  cout << "####################################################################" << endl;
  cout << "Fitted parameters and their uncertainties:" << endl;
  cout << "Parameter    Value      Uncertainty    Rel. uncertainty:" << endl;
  cout << "####################################################################" << endl;
  for (i=0; i<NXf; i++){
    re = 0.0;
    if (abs(Xfopt[i])>eps) re = abs(dXfopt[i] / Xfopt[i]);
    cout << "x" << i << " " << Xfopt[i] << " " << dXfopt[i] << " " << re << endl;
  }
  cout << "####################################################################" << endl;
#endif
