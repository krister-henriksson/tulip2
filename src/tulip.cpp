


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include <boost/format.hpp>

#include <cmath>

#include "chisq-basics.hpp"
#include "funcfit-basics.hpp"
#include "funcfit-conjgrad.hpp"
#include "funcfit-errors.hpp"
#include "funcfit-gauss-newton.hpp"
#include "funcfit-leve-marq.hpp"
#include "funcfit-powelldogleg.hpp"
#include "funcfit-powell.hpp"
#include "funcfit-simplexfit.hpp"
#include "funcfit-diffevol.hpp"
#include "funcfit-partswarm.hpp"
#include "funcfit-beecolony.hpp"
#include "funcfit-gravsearch.hpp"
#include "funcfit-simann.hpp"
#include "funcfit-moldyn.hpp"
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
#include "utils-errors.hpp"

#include "atomsystem.hpp"
#include "compound.hpp"
#include "compoundfit-list.hpp"
#include "elem-iacs.hpp"
#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mtwister.hpp"
#include "param-pot.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "potinfo.hpp"
#include "specs-fit-prop-pot.hpp"
#include "get-comp-prop.hpp"
#include "report.hpp"
#include "lattice-simple.hpp"
#include "get-ini-fit-data.hpp"


#include <omp.h>
#include "omp-basics.hpp"


#include <ctime>


#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;


using namespace utils;
using namespace funcfit;
using namespace constants;
using namespace physconst;
using boost::format;




#define VERSION 2.0


int main(int argc, char *argv[]){
  // double eps = numeric_limits<double>::epsilon();
  int i,j,k,p,iref, ivec;
  bool run_refonly, use_relonly, ini_fit_data_mode, debug_forces, debug_pressure, report_mds_steps;
  Vector<bool> debug_fit_prop(5, false),debug_fit_pot(5,false);
  std::string arg, potfile(""), geomfile(""), specsfile("");
  bool potfileOK, geomfileOK, specsfileOK;
  std::istringstream sstream;
  std::string tulip_info, dformat("%15.10f");

  OMP_Info omp_info;
  int omp_tid;
  // Default:
  int omp_nt_try = 1;

  double ts1 = omp_get_wtime();
  //clock_t clock1 = std::clock();

  bool u_d_xyz_fmt=true;
  std::string d_xyz_fmt="extxyz";




  tulip_info = "TULIP version " + tostring(VERSION) + " (c) Krister Henriksson 2013-";


  potfileOK = geomfileOK = specsfileOK = false;
  run_refonly = false;
  use_relonly = false;
  ini_fit_data_mode = false;
  debug_forces = false;
  debug_pressure = false;
  report_mds_steps = false;
  debug_fit_pot[0]  = true;
  


  /* ###############################################################################
     Parse the options.
     ############################################################################### */
  for (i=1; i<argc; i++){

    if (string(argv[i])=="-pf"){
      arg = std::string(argv[i+1]); sstream.str(arg); i++;
      sstream >> potfile;
      sstream.clear();
      potfileOK = true;
    }
    else if (string(argv[i])=="-gf"){
      arg = std::string(argv[i+1]); sstream.str(arg); i++;
      sstream >> geomfile;
      sstream.clear();
      geomfileOK = true;
    }
    else if (string(argv[i])=="-sf"){
      arg = std::string(argv[i+1]); sstream.str(arg); i++;
      sstream >> specsfile;
      sstream.clear();
      specsfileOK = true;
    }
    else if (string(argv[i])=="-ro"){
      run_refonly = true;
    }
    else if (string(argv[i])=="-nof"){
      use_relonly = true;
    }
    else if (string(argv[i])=="-xyz"){
      d_xyz_fmt   = "xyz";
      u_d_xyz_fmt = true;      
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
    else if (string(argv[i])=="-mif"){
      ini_fit_data_mode = true;
    }
    else if (string(argv[i])=="-omp"){
      omp_nt_try = omp_info.nt_max();
    }
    else if (string(argv[i])=="-omp_nt"){
      arg = std::string(argv[i+1]); sstream.str(arg); i++;
      sstream >> omp_nt_try;
      sstream.clear();
    }
    
  }



  if (potfileOK == false && geomfileOK == false){

    std::cout << tulip_info << std::endl;
    std::cout << "Purpose: Fit data to an interatomic potential." << std::endl;
    std::cout << "Usage:" << std::endl;
    std::cout << "     " << std::string(argv[0]) << " arguments [options]" << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "     -pf file           Path to file containing potential information." << std::endl;
    std::cout << "     -gf file           Path to file containing geometry information." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Options:" << std::endl;
    std::cout << "     -sf file           Path to file containing technical specifications about the calculations." << std::endl;

    std::cout << "     -ro                Only calculate properties of reference compounds, then exit. Default: not used." << std::endl;
    std::cout << "     -nof               Only calculate properties of reference and read-in compounds, then exit. Default: not used." << std::endl;
    std::cout << "     -xyz               Use traditional XYZ format when writing XYZ files. Default: not used." << std::endl;
    std::cout << "                        The extended XYZ format (http://jrkermode.co.uk/quippy/io.html#extendedxyz)" << std::endl;
    std::cout << "                        is used by default." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "     -dfitpropn         Show information about fitting of properties. Here 'n' must be" << std::endl;
    std::cout << "                        an integer. Supported: 0-4. 0: debug fitting method. 1-4: debug deeper." << std::endl;
    std::cout << "                        lying methods used by the fitting method. Default: not used" << std::endl;
    std::cout << "                        NOTE: 0 also shows some info about the initial Chi^2 object." << std::endl;
    std::cout << "     -dfitpotn          Show information about fitting of potentials. Here 'n' have a similar" << std::endl;
    std::cout << "                        role as for fitting of the properties." << std::endl;
    std::cout << "                        NOTE 1: 0 also shows some info about the initial Chi^2 object." << std::endl;
    std::cout << "                        NOTE 2: 'fitpot0' is always set to true, others are false by default." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "     -dforces           Debug the forces. Default: not used" << std::endl;
    std::cout << "     -dpressure         Debug the pressure. Default: not used" << std::endl;
    std::cout << "     -dmdsprop          Debug MDS runs of the structures. Default: not used" << std::endl;

    std::cout << "     -dall              Activate all debugging options (top level only). Default: not used" << std::endl;

    std::cout << "" << std::endl;
    std::cout << "     -mif               Suggest an initial fit and exit. Default: not used" << std::endl;
    std::cout << "                        ABOP: Put D0=0.0 for the binary interaction you want to fit. Keep all other" << std::endl;
    std::cout << "                        parametrizations at their normal values." << std::endl;
    std::cout << "" << std::endl;

#if 0
    std::cout << "" << std::endl;
    std::cout << "     -omp               Request maximal number of threads ("
	 << omp_info.nt_max() << ") for any OpenMP parts." << std::endl;
    std::cout << "     -omp_nt num        Request 'num' number of threads for any OpenMP parts. Default: "
	 << omp_info.nt_use() << std::endl;
    std::cout << "" << std::endl;
#endif

    return 0;
  }





  
  if (! potfileOK) aborterror("Error: Potential information file not specified. Exiting.");
  if (! geomfileOK) aborterror("Error: Geometry information file not specified. Exiting.");
  //  if (! specsfileOK) aborterror("Error: Specifications file not specified. Exiting.");

  std::cout << "Potential info file     : " << potfile << std::endl;
  std::cout << "Geometry info file      : " << geomfile << std::endl;
  std::cout << "Specifications info file: " << specsfile << std::endl;


  std::cout << "|||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;

  omp_info.nt_use( omp_nt_try );
#pragma omp parallel private(omp_tid)
  {
    omp_tid = omp_info.tid();
    if (omp_tid==0){
      std::cout << "OpenMP threads used     : " << omp_info.nt_use() << std::endl;
      std::cout << "OpenMP threads max count: " << omp_info.nt_max() << std::endl;
    }
  }

  std::cout << "|||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;




  std::cout << "**************************************************************" << std::endl;
  std::cout << "**************************************************************" << std::endl;
  std::cout << "" << std::endl;
  std::cout << tulip_info << std::endl;
  std::cout << "" << std::endl;
  std::cout << "**************************************************************" << std::endl;
  std::cout << "**************************************************************" << std::endl;

  std::cout << "" << std::endl;
  std::cout << "STARTING UP" << std::endl;
  std::cout << "  1. Constructing potential information object using potential and specifications files ..." << std::endl;

  PotentialInformationFit potinfo(potfile, specsfile);
  Elements elem = potinfo.elem;

  std::cout << "  2. Constructing list of compounds to use for fitting, using geometry information file ..." << std::endl;
  CompoundListFit complistfit(elem, potinfo.specs_prop.mds_specs, geomfile);
  int ncomp = complistfit.compounds.size();

  std::cout << "  3. Constructing parameter object ..." << std::endl;
  ParamPot param( &potinfo );

  std::cout << "" << std::endl;
  std::cout << "STARTUP COMPLETE" << std::endl;
  std::cout << "" << std::endl;


  potinfo.omp_info = omp_info;


  potinfo.specs_prop.mds_specs_common.debug_forces         = debug_forces;
  potinfo.specs_prop.mds_specs_common.debug_pressure       = debug_pressure;
  potinfo.specs_prop.mds_specs_common.report_step          = report_mds_steps;
  potinfo.specs_prop.mds_specs_common.use_def_dump_xyz_fmt = u_d_xyz_fmt;
  potinfo.specs_prop.mds_specs_common.def_dump_xyz_fmt     = d_xyz_fmt;




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





  // ################################################################################
  // ################################################################################
  //
  // Report on all settings ............
  //
  // ################################################################################
  // ################################################################################

  int nref;
  std::string s1, s2, s3;

  // ################################################################################  
  // Elements
  // ################################################################################
  std::cout << "" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  std::cout << "Elements:" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  nref = potinfo.elem.nelem();
  for (i=0; i<nref; ++i){
    s1 = potinfo.elem.idx2name(i);
    std::cout << s1
	 << " element index " << potinfo.elem.name2idx(s1)
	 << " atom type " << potinfo.elem.atomtype(s1)
	 << " mass (amu)  " << potinfo.elem.mass(s1)
	 << " reference lattice " << potinfo.elem.reflat(s1) << std::endl;
  }

  // ################################################################################  
  // Interactions
  // ################################################################################  
  std::cout << "" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  std::cout << "Interactions:" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  nref = potinfo.elem.nelem();
  for (i=0; i<nref; ++i){
    for (j=i; j<nref; ++j){
      s1 = potinfo.elem.idx2name(i);
      s2 = potinfo.elem.idx2name(j);
      
      std::cout << s1 << "-" << s2 << " iac " << potinfo.iacs.name(s1,s2)
	   << "(" << potinfo.basepot(s1,s2) << ")"
	   << " is fittable? " << potinfo.is_fittable(s1,s2)
	   << " uses reppot? " << potinfo.use_reppot(s1,s2);

      /*
      if (potinfo.use_reppot(s1,s2)){
	ivec = potinfo.reppot_vecidx(s1,s2);
	cout << " bermi= "  << potinfo.pot_Reppot[ivec].bfermi
	     << " rfermi= " << potinfo.pot_Reppot[ivec].rfermi;
      }
      */


      std::cout << std::endl;
    }
  }

  std::cout << "" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  std::cout << "FIXED interactions which are parametrized:" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  report_pot( &potinfo, false, true );

  std::cout << "" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  std::cout << "FITTABLE interactions (parametrized):" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  report_pot( &potinfo, true, false );




  // ################################################################################
  // Specifications for compounds
  // ################################################################################
  std::cout << "" << std::endl;
  std::cout << "##########################################################################" << std::endl;
  std::cout << "Specifications for calculation of properties of compounds:" << std::endl;
  std::cout << "##########################################################################" << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  std::cout << "General" << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;

  std::cout << "Random number seed                                 : " << potinfo.specs_prop.seed << std::endl;
  std::cout << "Fitting method                                     : " << potinfo.specs_prop.fitmet << std::endl;
  std::cout << "Min. iterations                                    : " << potinfo.specs_prop.nitermin << std::endl;
  std::cout << "Max. iterations                                    : " << potinfo.specs_prop.nitermax << std::endl;
  std::cout << "Restart at iteration (if positive)                 : " << potinfo.specs_prop.niterrestart << std::endl;
  std::cout << "Tolerance (abs.) for convergence of ChiSq                 : " << potinfo.specs_prop.functolabs << std::endl;
  std::cout << "Tolerance (rel.) for convergence of ChiSq changes         : " << potinfo.specs_prop.functolrel << std::endl;
  std::cout << "Tolerance (abs.) for convergence of ChiSq gradient        : " << potinfo.specs_prop.gradtolabs << std::endl;
  std::cout << "Tolerance (abs.) for convergence of step length           : " << potinfo.specs_prop.steptolabs << std::endl;
  std::cout << "Tolerance (rel.) for convergence of step length changes   : " << potinfo.specs_prop.steptolrel << std::endl;
  std::cout << "DOG-LEG: Initial trust region radius               : " << potinfo.specs_prop.dogleg_radius << std::endl;
  std::cout << "DOG-LEG: Smallest allowed trust region radius      : " << potinfo.specs_prop.dogleg_minradius << std::endl;
  std::cout << "SIMPLEX: Displacement when creating initial simplex: " << potinfo.specs_prop.simplex_delta << std::endl;
  std::cout << "SIMANN: Relative displacement when creating trial  : " << potinfo.specs_prop.simann_delta_rel << std::endl;
  std::cout << "MolDynFit: min_dx                                  : " << potinfo.specs_prop.moldyn_min_dx << std::endl;
  std::cout << "MolDynFit: max_dx                                  : " << potinfo.specs_prop.moldyn_max_dx << std::endl;
  std::cout << "Penalty function: barrier                          : " << potinfo.specs_prop.barrier_scale << std::endl;
  std::cout << "Penalty function: use barrier rescaling?           : " << potinfo.specs_prop.use_barrier_rescaling << std::endl;
  std::cout << "Normalize data using initial values?               : " << potinfo.specs_prop.use_data_scales << std::endl;
  std::cout << "Debug: level0                                      : " << potinfo.specs_prop.debug_fit_level0 << std::endl;
  std::cout << "Debug: level1                                      : " << potinfo.specs_prop.debug_fit_level1 << std::endl;
  std::cout << "Debug: level2                                      : " << potinfo.specs_prop.debug_fit_level2 << std::endl;
  std::cout << "Debug: level3                                      : " << potinfo.specs_prop.debug_fit_level3 << std::endl;
  std::cout << "Debug: level4                                      : " << potinfo.specs_prop.debug_fit_level4<< std::endl;
  std::cout << "Lattice tolerance                                  : " << potinfo.specs_prop.lattol << std::endl;

  std::cout << std::endl;
  std::cout << "Bulk modulus: Relax strained systems?              : " << potinfo.specs_prop.BM_rel_sys << std::endl;
  std::cout << "Bulk modulus: Minimum strain (e.g. -0.01)          : " << potinfo.specs_prop.BM_fmin << std::endl;
  std::cout << "Bulk modulus: Maximum strain (e.g.  0.01)          : " << potinfo.specs_prop.BM_fmax << std::endl;
  std::cout << "Bulk modulus: Number of strain points (e.g. 10)    : " << potinfo.specs_prop.BM_Nf << std::endl;
  std::cout << "Bulk modulus: Multiplicative uncertainty factor (e.g. 0.10): " << potinfo.specs_prop.BM_ef << std::endl;
  std::cout << std::endl;
  std::cout << "Elastic moduli: Relax strained systems?            : " << potinfo.specs_prop.C_rel_sys << std::endl;
  std::cout << "Elastic moduli: Minimum strain (e.g. -0.01)        : " << potinfo.specs_prop.C_fmin << std::endl;
  std::cout << "Elastic moduli: Maximum strain (e.g.  0.01)        : " << potinfo.specs_prop.C_fmax << std::endl;
  std::cout << "Elastic moduli: Number of strain points (e.g. 10)  : " << potinfo.specs_prop.C_Nf << std::endl;
  std::cout << "Elastic moduli: Multiplicative uncertainty factor (e.g. 0.10): " << potinfo.specs_prop.C_ef << std::endl;


  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "Settings for MD simulations of read-in compounds" << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;

  std::cout << "Skin thickness of neighbor list                         : " << potinfo.specs_prop.mds_specs.skint << std::endl;
  std::cout << "Seed for random numbers                                 : " << potinfo.specs_prop.mds_specs.seed << std::endl;
  std::cout << "Starting time for simulation (usually 0.0)              : " << potinfo.specs_prop.mds_specs.tstart << std::endl;
  std::cout << "Ending time for simulation                              : " << potinfo.specs_prop.mds_specs.tend << std::endl;
  std::cout << "Dump info every n:th step, n is                         : " << potinfo.specs_prop.mds_specs.ndump << std::endl;
  std::cout << "Starting/desired temperature (e.g. 300.0)               : " << potinfo.specs_prop.mds_specs.Tstart << std::endl;
  std::cout << "Initial time step                                       : " << potinfo.specs_prop.mds_specs.dt << std::endl;
  std::cout << "  Maximum time step                                     : " << potinfo.specs_prop.mds_specs.max_dt << std::endl;
  std::cout << "  Maximum allowed energy change for an atom/timestep    : " << potinfo.specs_prop.mds_specs.max_dE << std::endl;
  std::cout << "  Maximum allowed distance traveled for an atom/timestep: " << potinfo.specs_prop.mds_specs.max_dr << std::endl;
  std::cout << "Use temperature control (Berendsen)?                    : " << potinfo.specs_prop.mds_specs.use_Tcontrol << std::endl;
  std::cout << "  T control: time constant (fs)                         : " << potinfo.specs_prop.mds_specs.btc_tau << std::endl;
  std::cout << "  T control: desired temperature                        : " << potinfo.specs_prop.mds_specs.btc_T0 << std::endl;
  std::cout << "Use quench?                                             : " << potinfo.specs_prop.mds_specs.use_quench << std::endl;
  std::cout << "  Quench: start time (fs)                               : " << potinfo.specs_prop.mds_specs.quench_tstart << std::endl;
  std::cout << "  Quench: rate (K/fs) (always pos. value)               : " << potinfo.specs_prop.mds_specs.quench_rate << std::endl;
  std::cout << "Use pressure control (Berendsen)?                       : " << potinfo.specs_prop.mds_specs.use_Pcontrol << std::endl;
  std::cout << "  P control: time constant (fs)                         : " << potinfo.specs_prop.mds_specs.bpc_tau << std::endl;
  std::cout << "  P control: desired pressure (GPa)                     : " << potinfo.specs_prop.mds_specs.bpc_P0 << std::endl; 
  std::cout << "  P control: scale                                      : " << potinfo.specs_prop.mds_specs.bpc_scale << std::endl; 
  std::cout << "    The scale is usually equal to the approximate bulk modulus (GPa)." << std::endl;
  std::cout << "Terminate program is T > T_limit?                       : " << potinfo.specs_prop.mds_specs.use_error_T_gt << std::endl; 
  std::cout << "  T_limit                                               : " << potinfo.specs_prop.mds_specs.error_T_gt << std::endl; 
  std::cout << "Terminate program is dt < dt_limit?                     : " << potinfo.specs_prop.mds_specs.use_error_dt_lt << std::endl; 
  std::cout << "  dt_limit                                              : " << potinfo.specs_prop.mds_specs.error_dt_lt << std::endl; 
  std::cout << "Terminate program is boxlen[1/2/3] > bl_limit?          : " << potinfo.specs_prop.mds_specs.use_error_boxlen_gt << std::endl; 
  std::cout << "  bl_limit                                              : " << potinfo.specs_prop.mds_specs.error_boxlen_gt << std::endl; 
  std::cout << "" << std::endl;

  std::cout << "Options: external relaxation?                           : " << potinfo.specs_prop.mds_specs.ext_relax << std::endl; 
  std::cout << "Options: quench always?                                 : " << potinfo.specs_prop.mds_specs.quench_always << std::endl; 
  std::cout << "  If true, all velocities zeroed at every time step." << std::endl;
  std::cout << "" << std::endl;

  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "Settings for MD simulations of reference compounds" << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;

  std::cout << "Skin thickness of neighbor list                         : " << potinfo.specs_prop.mds_specs_ref.skint << std::endl;
  std::cout << "Seed for random numbers                                 : " << potinfo.specs_prop.mds_specs_ref.seed << std::endl;
  std::cout << "Starting time for simulation (usually 0.0)              : " << potinfo.specs_prop.mds_specs_ref.tstart << std::endl;
  std::cout << "Ending time for simulation                              : " << potinfo.specs_prop.mds_specs_ref.tend << std::endl;
  std::cout << "Dump info every n:th step, n is                         : " << potinfo.specs_prop.mds_specs_ref.ndump << std::endl;
  std::cout << "Starting/desired temperature (e.g. 300.0)               : " << potinfo.specs_prop.mds_specs_ref.Tstart << std::endl;
  std::cout << "Initial time step                                       : " << potinfo.specs_prop.mds_specs_ref.dt << std::endl;
  std::cout << "  Maximum time step                                     : " << potinfo.specs_prop.mds_specs_ref.max_dt << std::endl;
  std::cout << "  Maximum allowed energy change for an atom/timestep    : " << potinfo.specs_prop.mds_specs_ref.max_dE << std::endl;
  std::cout << "  Maximum allowed distance traveled for an atom/timestep: " << potinfo.specs_prop.mds_specs_ref.max_dr << std::endl;
  std::cout << "Use temperature control (Berendsen)?                    : " << potinfo.specs_prop.mds_specs_ref.use_Tcontrol << std::endl;
  std::cout << "  T control: time constant (fs)                         : " << potinfo.specs_prop.mds_specs_ref.btc_tau << std::endl;
  std::cout << "  T control: desired temperature                        : " << potinfo.specs_prop.mds_specs_ref.btc_T0 << std::endl;
  std::cout << "Use quench?                                             : " << potinfo.specs_prop.mds_specs_ref.use_quench << std::endl;
  std::cout << "  Quench: start time (fs)                               : " << potinfo.specs_prop.mds_specs_ref.quench_tstart << std::endl;
  std::cout << "  Quench: rate (K/fs) (always pos. value)               : " << potinfo.specs_prop.mds_specs_ref.quench_rate << std::endl;
  std::cout << "Use pressure control (Berendsen)?                       : " << potinfo.specs_prop.mds_specs_ref.use_Pcontrol << std::endl;
  std::cout << "  P control: time constant (fs)                         : " << potinfo.specs_prop.mds_specs_ref.bpc_tau << std::endl;
  std::cout << "  P control: desired pressure (GPa)                     : " << potinfo.specs_prop.mds_specs_ref.bpc_P0 << std::endl; 
  std::cout << "  P control: scale                                      : " << potinfo.specs_prop.mds_specs_ref.bpc_scale << std::endl; 
  std::cout << "    The scale is usually equal to the approximate bulk modulus (GPa)." << std::endl;
  std::cout << "Terminate program is T > T_limit?                       : " << potinfo.specs_prop.mds_specs_ref.use_error_T_gt << std::endl; 
  std::cout << "  T_limit                                               : " << potinfo.specs_prop.mds_specs_ref.error_T_gt << std::endl; 
  std::cout << "Terminate program is dt < dt_limit?                     : " << potinfo.specs_prop.mds_specs_ref.use_error_dt_lt << std::endl; 
  std::cout << "  dt_limit                                              : " << potinfo.specs_prop.mds_specs_ref.error_dt_lt << std::endl; 
  std::cout << "Terminate program is boxlen[1/2/3] > bl_limit?          : " << potinfo.specs_prop.mds_specs_ref.use_error_boxlen_gt << std::endl; 
  std::cout << "  bl_limit                                              : " << potinfo.specs_prop.mds_specs_ref.error_boxlen_gt << std::endl; 

  std::cout << "" << std::endl;
  std::cout << "Options: external relaxation?                           : " << potinfo.specs_prop.mds_specs_ref.ext_relax << std::endl; 
  std::cout << "Options: quench always?                                 : " << potinfo.specs_prop.mds_specs_ref.quench_always << std::endl; 
  std::cout << "  If true, all velocities zeroed at every time step." << std::endl;
  std::cout << "" << std::endl;

  std::cout << std::endl;
  std::cout << "**************************************************************************" << std::endl;
  std::cout << "Common settings for MD simulations of any compound" << std::endl;
  std::cout << "**************************************************************************" << std::endl;
  std::cout << std::endl;

  std::cout << "Report step?                                   : " << potinfo.specs_prop.mds_specs_common.report_step << std::endl;
  std::cout << "  If true, writes out physical info of the system at every time step." << std::endl;
  std::cout << "Debug forces?                                  : " << potinfo.specs_prop.mds_specs_common.debug_forces << std::endl; 
  std::cout << "Debug pressure?                                : " << potinfo.specs_prop.mds_specs_common.debug_pressure << std::endl; 
  std::cout << "Use default XYZ format?                        : " << potinfo.specs_prop.mds_specs_common.use_def_dump_xyz_fmt << std::endl;
  std::cout << "Default XYZ format?                            : " << potinfo.specs_prop.mds_specs_common.def_dump_xyz_fmt << std::endl;
  std::cout << "" << std::endl;

  std::cout << std::endl;
  std::cout << "//////////////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << "Special settings for various types of MDS cases" << std::endl;
  std::cout << "//////////////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << std::endl;

  std::cout << "'''''''''''''''''''' Elastic properties ''''''''''''''''''''" << std::endl;
  std::cout << "  quench_always = true" << std::endl;
  std::cout << std::endl;



  // ################################################################################  
  // Specifications for potentials
  // ################################################################################  
  std::cout << "" << std::endl;
  std::cout << "##########################################################################" << std::endl;
  std::cout << "Specifications for potential fitting:" << std::endl;
  std::cout << "##########################################################################" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Random number seed                                 : " << potinfo.specs_pot.seed << std::endl;
  std::cout << "Fitting method                                     : " << potinfo.specs_pot.fitmet << std::endl;
  std::cout << "Min. iterations                                    : " << potinfo.specs_pot.nitermin << std::endl;
  std::cout << "Max. iterations                                    : " << potinfo.specs_pot.nitermax << std::endl;
  std::cout << "Restart at iteration (if positive)                 : " << potinfo.specs_pot.niterrestart << std::endl;
  std::cout << "Tolerance (abs.) for convergence of ChiSq                 : " << potinfo.specs_pot.functolabs << std::endl;
  std::cout << "Tolerance (rel.) for convergence of ChiSq changes         : " << potinfo.specs_pot.functolrel << std::endl;
  std::cout << "Tolerance (abs.) for convergence of ChiSq gradient        : " << potinfo.specs_pot.gradtolabs << std::endl;
  std::cout << "Tolerance (abs.) for convergence of step length           : " << potinfo.specs_pot.steptolabs << std::endl;
  std::cout << "Tolerance (rel.) for convergence of step length changes   : " << potinfo.specs_pot.steptolrel << std::endl;
  std::cout << "DOG-LEG: Initial trust region radius               : " << potinfo.specs_pot.dogleg_radius << std::endl;
  std::cout << "DOG-LEG: Smallest allowed trust region radius      : " << potinfo.specs_pot.dogleg_minradius << std::endl;
  std::cout << "SIMPLEX: Displacement when creating initial simplex: " << potinfo.specs_pot.simplex_delta << std::endl;
  std::cout << "SIMANN: Relative displacement when creating trial  : " << potinfo.specs_pot.simann_delta_rel << std::endl;
  std::cout << "MolDynFit: min_dx                                  : " << potinfo.specs_pot.moldyn_min_dx << std::endl;
  std::cout << "MolDynFit: max_dx                                  : " << potinfo.specs_pot.moldyn_max_dx << std::endl;
  std::cout << "Penalty function: barrier                          : " << potinfo.specs_pot.barrier_scale << std::endl;
  std::cout << "Penalty function: use barrier rescaling?           : " << potinfo.specs_pot.use_barrier_rescaling << std::endl;
  std::cout << "Normalize data using initial values?               : " << potinfo.specs_pot.use_data_scales << std::endl;
  std::cout << "Debug: level0                                      : " << potinfo.specs_pot.debug_fit_level0 << std::endl;
  std::cout << "Debug: level1                                      : " << potinfo.specs_pot.debug_fit_level1 << std::endl;
  std::cout << "Debug: level2                                      : " << potinfo.specs_pot.debug_fit_level2 << std::endl;
  std::cout << "Debug: level3                                      : " << potinfo.specs_pot.debug_fit_level3 << std::endl;
  std::cout << "Debug: level4                                      : " << potinfo.specs_pot.debug_fit_level4 << std::endl;






  /*
  std::cout << myomp_get_chunksize(sizeof(double)) << std::endl;
  exit(1);
  */




  // ################################################################################  
  // Compounds
  // ################################################################################  

  // Some adjustments ...
  for (i=0; i<ncomp; ++i){
    complistfit.compounds[i].check_and_fix_uses();
  }

  // Report:
  std::cout << "" << std::endl;
  std::cout << "##########################################################################" << std::endl;
  std::cout << "Compounds used for fitting:" << std::endl;
  std::cout << "##########################################################################" << std::endl;
  std::cout << "" << std::endl;
  
  for (i=0; i<ncomp; ++i){
    std::cout << "Compound " << i+1 << " of " << ncomp << std::endl;
    std::cout << "  Name                : " << complistfit.compounds[i].name << std::endl;
    std::cout << "  Read from file      : " << complistfit.compounds[i].filename << std::endl;
    std::cout << "  PBC                 : " << complistfit.compounds[i].pbc[0] << " "
	 << complistfit.compounds[i].pbc[1] << " "
	 << complistfit.compounds[i].pbc[2] << std::endl;
    std::cout << "  Crystal system      : " << complistfit.compounds[i].csystem << std::endl;
    std::cout << "  Scale factor        : " << complistfit.compounds[i].scalefactor << std::endl;
    std::cout << "  Use internal format?: " << complistfit.compounds[i].use_int << std::endl;
    std::cout << "  Direction vector 1  : " 
	 << format("%15.10f ") % complistfit.compounds[i].u1_vec[0]
	 << format("%15.10f ") % complistfit.compounds[i].u1_vec[1]
	 << format("%15.10f")  % complistfit.compounds[i].u1_vec[2] << std::endl;
    std::cout << "  Direction vector 2  : " 
	 << format("%15.10f ") % complistfit.compounds[i].u2_vec[0]
	 << format("%15.10f ") % complistfit.compounds[i].u2_vec[1]
	 << format("%15.10f")  % complistfit.compounds[i].u2_vec[2] << std::endl;
    std::cout << "  Direction vector 3  : " 
	 << format("%15.10f ") % complistfit.compounds[i].u3_vec[0]
	 << format("%15.10f ") % complistfit.compounds[i].u3_vec[1]
	 << format("%15.10f")  % complistfit.compounds[i].u3_vec[2] << std::endl;
    std::cout << "  Number of basis atoms: " << complistfit.compounds[i].nbasis << std::endl;

    std::cout << "" << std::endl;
    std::cout << "Desired N values for compound construction      : "
	 << complistfit.compounds[i].Ndesired[0] << " "
	 << complistfit.compounds[i].Ndesired[1] << " "
	 << complistfit.compounds[i].Ndesired[2] << std::endl;
    std::cout << "Desired even N values for compound construction?: "
	 << complistfit.compounds[i].Neven_desired[0] << " "
	 << complistfit.compounds[i].Neven_desired[1] << " "
	 << complistfit.compounds[i].Neven_desired[2] << std::endl;
    std::cout << "Desired odd N values for compound construction? : "
	 << complistfit.compounds[i].Nodd_desired[0] << " "
	 << complistfit.compounds[i].Nodd_desired[1] << " "
	 << complistfit.compounds[i].Nodd_desired[2] << std::endl;

    std::cout << "" << std::endl;
    std::cout << "Compound is a reference for changes in Ecoh (Ecoh_delta) ? : "
	 << complistfit.compounds[i].Ecoh_delta_refcomp << std::endl;

    // constraints
    int nic=0;
    for (int ic=0; ic<complistfit.compounds[i].nbasis; ++ic){
      if (complistfit.compounds[i].basis_is_fixed[ic]) nic++;
      if (complistfit.compounds[i].basis_freedir[ic].size()==3) nic++;
      if (complistfit.compounds[i].basis_freeplane[ic].size()==3) nic++;
    }
    std::cout << "" << std::endl;
    std::cout << "Number of constrained atoms: " << nic << std::endl;
    for (int ic=0; ic<complistfit.compounds[i].nbasis; ++ic){
      if (complistfit.compounds[i].basis_is_fixed[ic])
	std::cout << "  Atom " << ic
		  << " is fixed (0=false, 1=true)   " << complistfit.compounds[i].basis_is_fixed[ic] << std::endl;
      if (complistfit.compounds[i].basis_freedir[ic].size()==3){
	std::cout << "  Atom " << ic
		  << " has a free direction         " << complistfit.compounds[i].basis_freedir[ic] << std::endl;
      }
      if (complistfit.compounds[i].basis_freeplane[ic].size()==3){
	std::cout << "  Atom " << ic
		  << " has a free plane with normal " << complistfit.compounds[i].basis_freeplane[ic] << std::endl;
      }
    }


    std::cout << "" << std::endl;
    std::string dumpfile(complistfit.compounds[i].filename + ".xyz");
    std::cout << "  Atom system will be written to file of format Extended XYZ for debugging purposes: " << dumpfile << std::endl;
    std::ofstream fout;
    fout.open(dumpfile.c_str());
    fout << complistfit.compounds[i].basis_vecs.size() << std::endl;
    fout << "Lattice=\""
	 << format("%10.6f ") % complistfit.compounds[i].u1_vec[0]
	 << format("%10.6f ") % complistfit.compounds[i].u1_vec[1]
	 << format("%10.6f")  % complistfit.compounds[i].u1_vec[2]
	 << format("%10.6f ") % complistfit.compounds[i].u2_vec[0]
	 << format("%10.6f ") % complistfit.compounds[i].u2_vec[1]
	 << format("%10.6f")  % complistfit.compounds[i].u2_vec[2]
	 << format("%10.6f ") % complistfit.compounds[i].u3_vec[0]
	 << format("%10.6f ") % complistfit.compounds[i].u3_vec[1]
	 << format("%10.6f")  % complistfit.compounds[i].u3_vec[2]
	 << "\" Properties=species:S:1:pos:R:3:index:I:1 Time=0.0" << std::endl;
    for (j=0; j<complistfit.compounds[i].basis_vecs.size(); ++j){
      fout << complistfit.compounds[i].basis_elems[j]
	   << "  "
	   << format("%15.10f ") % complistfit.compounds[i].basis_vecs[j][0]
	   << format("%15.10f ") % complistfit.compounds[i].basis_vecs[j][1]
	   << format("%15.10f")  % complistfit.compounds[i].basis_vecs[j][2]
	   << "  "
	   << potinfo.elem.name2idx( complistfit.compounds[i].basis_elems[j] )
	   << "  "
	   << j << std::endl;
    }
    fout.close();
    fout.clear();



    std::cout << "  Compound-specific MD settings" << std::endl;
    std::cout << "  -----------------------------" << std::endl;
    std::cout << "Skin thickness of neighbor list                         : " << complistfit.compounds[i].mds_specs.skint << std::endl;
    std::cout << "Seed for random numbers                                 : " << complistfit.compounds[i].mds_specs.seed << std::endl;
    std::cout << "Starting time for simulation (usually 0.0)              : " << complistfit.compounds[i].mds_specs.tstart << std::endl;
    std::cout << "Ending time for simulation                              : " << complistfit.compounds[i].mds_specs.tend << std::endl;
    std::cout << "Dump info every n:th step, n is                         : " << complistfit.compounds[i].mds_specs.ndump << std::endl;
    std::cout << "Starting temperature (e.g. 300.0)                       : " << complistfit.compounds[i].mds_specs.Tstart << std::endl;
    std::cout << "Initial time step                                       : " << complistfit.compounds[i].mds_specs.dt << std::endl;
    std::cout << "  Maximum time step                                     : " << complistfit.compounds[i].mds_specs.max_dt << std::endl;
    std::cout << "  Maximum allowed energy change for an atom/timestep    : " << complistfit.compounds[i].mds_specs.max_dE << std::endl;
    std::cout << "  Maximum allowed distance traveled for an atom/timestep: " << complistfit.compounds[i].mds_specs.max_dr << std::endl;
    std::cout << "Use temperature control (Berendsen)?                    : " << complistfit.compounds[i].mds_specs.use_Tcontrol << std::endl;
    std::cout << "  T control: time constant (fs)                         : " << complistfit.compounds[i].mds_specs.btc_tau << std::endl;
    std::cout << "  T control: desired temperature                        : " << complistfit.compounds[i].mds_specs.btc_T0 << std::endl;
    std::cout << "Use quench?                                             : " << complistfit.compounds[i].mds_specs.use_quench << std::endl;
    std::cout << "  Quench: start time (fs)                               : " << complistfit.compounds[i].mds_specs.quench_tstart << std::endl;
    std::cout << "  Quench: rate (K/fs) (always pos. value)               : " << complistfit.compounds[i].mds_specs.quench_rate << std::endl;
    std::cout << "Use pressure control (Berendsen)?                       : " << complistfit.compounds[i].mds_specs.use_Pcontrol << std::endl;
    std::cout << "  P control: time constant (fs)                         : " << complistfit.compounds[i].mds_specs.bpc_tau << std::endl;
    std::cout << "  P control: desired pressure (GPa)                     : " << complistfit.compounds[i].mds_specs.bpc_P0 << std::endl; 
    std::cout << "  P control: scale                                      : " << complistfit.compounds[i].mds_specs.bpc_scale << std::endl; 
    std::cout << "    The scale is usually equal to the approximate bulk modulus (GPa)." << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Options: external relaxation?                           : " << complistfit.compounds[i].mds_specs.ext_relax << std::endl; 
    std::cout << "Options: quench always?                                 : " << complistfit.compounds[i].mds_specs.quench_always << std::endl; 
    std::cout << "  If true, all velocities zeroed at every time step." << std::endl;

    std::cout << "--------------------------------------------------------------------------" << std::endl;
  }



  std::cout << "INFO: Number of free fitting prameters: " << param.NXfree() << std::endl;
  std::cout << "INFO: Number of data points           : " << complistfit.NData() << std::endl;
  std::cout << "--------------------------------------------------------------------------" << std::endl;


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

    if (complistfit.compounds[i].prop_use.Ecoh_delta){ DY.push_back(complistfit.compounds[i].prop_readin.Ecoh_delta);
      if (complistfit.compounds[i].use_u.Ecoh_delta){
	DUY.push_back(complistfit.compounds[i].prop_u.Ecoh_delta);
	DWY.push_back(-1.0);
      }
      else {
	DUY.push_back(-1.0);
	DWY.push_back(complistfit.compounds[i].prop_w.Ecoh_delta);
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

    // Forces !!!
    if (complistfit.compounds[i].prop_use.frc){
      int nb = complistfit.compounds[i].basis_elems.size();
      for (int iat=0; iat<nb; ++iat){
	for (int k=0; k<3; ++k){

	  DY.push_back( complistfit.compounds[i].prop_readin.frc[iat][k] );

	  if (complistfit.compounds[i].use_u.frc){
	    DUY.push_back( complistfit.compounds[i].prop_u.frc[iat][k] );
	    DWY.push_back(-1.0);
	  }
	  else {
	    DUY.push_back(-1.0);
	    DWY.push_back( complistfit.compounds[i].prop_w.frc[iat][k] );
	  }
	}
      }

    }





  }


  report_prop( complistfit.compounds, std::cout, true );








  // ################################################################################
  // Get properties (at least cohesive energies) of reference structures !!!
  // ################################################################################

  std::cout << "" << std::endl;
  std::cout << "##########################################################################" << std::endl;
  std::cout << "Reference compounds" << std::endl;
  std::cout << "##########################################################################" << std::endl;
  std::cout << "" << std::endl;

  nref = potinfo.elem.nelem();

  potinfo.Ecoh_ref.resize(nref);
  for (iref=0; iref<nref; ++iref) potinfo.Ecoh_ref[iref] = 0;



  for (iref=0; iref<nref; ++iref){
    std::string sref = potinfo.elem.idx2name(iref);
    std::string latref = potinfo.elem.reflat(sref);

    if (latref=="skip"){
      std::cout << "*********************************************************************************" << std::endl;
      std::cout << "*********************************************************************************" << std::endl;
      std::cout << "WARNING: Skipping calculation of properties of the reference lattice "
		<< "for pure element " << sref << "." << std::endl;
      std::cout << "*********************************************************************************" << std::endl;
      std::cout << "*********************************************************************************" << std::endl;
      continue;
    }



    if (potinfo.is_fittable(sref,sref)){
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cout << "Elemental combination " << sref << "-" << sref << " is to be fitted."
	   << " No reference lattice calculation possible. Continuing." << std::endl;
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      continue;
    }


    Vector<CompoundStructureFit> cmpref(1, CompoundStructureFit());;
    double a,b,c;

    // ###############################################################
    // Create from generic model:
    // ###############################################################
    a = potinfo.elem.reflat_a(sref);
    b = potinfo.elem.reflat_b(sref);
    c = potinfo.elem.reflat_c(sref);
    cmpref[0].create_from_model(elem,latref,sref,sref, a,b,c);





    std::cout << "Reference compound " << iref+1 << " of " << nref << std::endl;
    std::cout << "  Name                : " << cmpref[0].name << std::endl;
    std::cout << "  PBC                 : " << cmpref[0].pbc[0] << " "
	 << cmpref[0].pbc[1] << " "
	 << cmpref[0].pbc[2] << std::endl;
    std::cout << "  Crystal system      : " << cmpref[0].csystem << std::endl;
    std::cout << "  Scale factor        : " << cmpref[0].scalefactor << std::endl;
    std::cout << "  Lattice parameters provided" << std::endl;
    std::cout << "                     a: " << a << std::endl;
    std::cout << "                     b: " << b << std::endl;
    std::cout << "                     c: " << c << std::endl;
    std::cout << "                => c/a: " << c/a << std::endl;
    std::cout << "  Use internal format?: " << cmpref[0].use_int << std::endl;
    std::cout << "  Direction vector 1  : " 
	 << format("%15.10f ") % cmpref[0].u1_vec[0]
	 << format("%15.10f ") % cmpref[0].u1_vec[1]
	 << format("%15.10f")  % cmpref[0].u1_vec[2] << std::endl;
    std::cout << "  Direction vector 2  : " 
	 << format("%15.10f ") % cmpref[0].u2_vec[0]
	 << format("%15.10f ") % cmpref[0].u2_vec[1]
	 << format("%15.10f")  % cmpref[0].u2_vec[2] << std::endl;
    std::cout << "  Direction vector 3  : " 
	 << format("%15.10f ") % cmpref[0].u3_vec[0]
	 << format("%15.10f ") % cmpref[0].u3_vec[1]
	 << format("%15.10f")  % cmpref[0].u3_vec[2] << std::endl;
    std::cout << "  Number of basis atoms: " << cmpref[0].nbasis << std::endl;






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

    cmpref[0].prop_use.Fmax = false;
    cmpref[0].prop_use.Pmax = false;
    cmpref[0].prop_use.displmax = false;
    cmpref[0].prop_use.frc = false;

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


    latsymm(cmpref);

    //cmpref[0].check_crystal_symm();
    cmpref[0].check_and_fix_Cij();



    // Calculate properties:
    get_comp_prop(param, cmpref);

    // Report:
    std::cout << "Properties of reference compound for element " << sref << ":" << std::endl;
    std::cout << "  Lattice parameter a     : " << format("%15.10f") % cmpref[0].prop_pred.a << std::endl;
    std::cout << "  Lattice parameter b     : " << format("%15.10f") % cmpref[0].prop_pred.b << std::endl;
    std::cout << "  Lattice parameter c     : " << format("%15.10f") % cmpref[0].prop_pred.c << std::endl;
    std::cout << "  Cohesive energy   Ecoh  : " << format("%15.10f") % cmpref[0].prop_pred.Ecoh
	 << " for element index " << potinfo.elem.name2idx(sref)
	 << " having atom type " << potinfo.elem.atomtype(sref)
	 << std::endl;

    potinfo.Ecoh_ref[ potinfo.elem.name2idx(sref) ] = cmpref[0].prop_pred.Ecoh;

    if (cmpref[0].prop_use.r0)
      std::cout << "  Dimer bond lenght r0    : " << format("%15.10f") % cmpref[0].prop_pred.r0 << std::endl;

    if (cmpref[0].prop_use.B)
      std::cout << "  Bulk modulus      B     : " << format("%15.10f") % cmpref[0].prop_pred.B << std::endl;
    if (cmpref[0].prop_use.Bp)
      std::cout << "  Pressure derivative of B: " << format("%15.10f") % cmpref[0].prop_pred.Bp << std::endl;

    if (anyC){
      if (cmpref[0].csystem=="cubic"){
	cout << "  Elastic constant     C11: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,0) << std::endl;
	cout << "  Elastic constant     C12: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,1) << std::endl;
	cout << "  Elastic constant     C44: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(3,3) << std::endl;
      }
      else if (cmpref[0].csystem=="hexagonal"){
	cout << "  Elastic constant     C11: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,0) << std::endl;
	cout << "  Elastic constant     C12: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,1) << std::endl;
	cout << "  Elastic constant     C13: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,2) << std::endl;
	cout << "  Elastic constant     C33: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(2,2) << std::endl;
	cout << "  Elastic constant     C44: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(3,3) << std::endl;
      }
      else if (cmpref[0].csystem=="orthorombic"){
	cout << "  Elastic constant     C11: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,0) << std::endl;
	cout << "  Elastic constant     C12: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,1) << std::endl;
	cout << "  Elastic constant     C13: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(0,2) << std::endl;
	cout << "  Elastic constant     C22: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(1,1) << std::endl;
	cout << "  Elastic constant     C23: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(1,2) << std::endl;
	cout << "  Elastic constant     C33: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(2,2) << std::endl;
	cout << "  Elastic constant     C44: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(3,3) << std::endl;
	cout << "  Elastic constant     C55: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(4,4) << std::endl;
	cout << "  Elastic constant     C66: " << format("%15.10f") % cmpref[0].prop_pred.C.elem(5,5) << std::endl;
      }
    }
    std::cout << "  Maximum force           : " << format("%15.10f") % cmpref[0].prop_pred.Fmax << std::endl;
    std::cout << "  Maximum pressure        : " << format("%15.10f") % cmpref[0].prop_pred.Pmax << std::endl;
    std::cout << "  Maximum displacement    : " << format("%15.10f") % cmpref[0].prop_pred.displmax << std::endl;


    std::cout << "--------------------------------------------------------------------------" << std::endl;




  }






  if (run_refonly){
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "Done with reference compounds. Quitting." << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    return EXIT_SUCCESS;
  }







  // ***********************************************************************************
  // ***********************************************************************************
  // ***********************************************************************************
  // ***********************************************************************************
  if (ini_fit_data_mode){

    std::cout << "Getting some info on the specified compounds ..." << std::endl;
    get_ini_fit_data(param, DX);


  }
  else if (use_relonly){

    // ################################################################################
    // A
    // ################################################################################

    //report_pot_prop( param, DX, DY, MDY );

    std::cout << "'Relax only' option used." << std::endl;
    std::cout << "1. Report on potential parameters:" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    report_pot( &potinfo, true, true );
    std::cout << "2. Relaxing read-in compounds and obtaining their properties ..." << std::endl;
    get_comp_prop(param, DX);
    std::cout << "----------------------------------------------------------------" << std::endl;
    report_prop( DX );
    std::cout << "Done." << std::endl;

  }
  else {

    // ################################################################################
    // B
    // ################################################################################
    // Potential fitting
    // ################################################################################


    // #############################################################
    // Form merit function object
    // #############################################################

    std::cout << "Setting up merit function ..." << std::endl;


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
    cs.ModelFuncPointer()  = get_comp_prop;
    cs.ReportFuncPointer() = report_pot_prop;

    cs.barrier_scale() = potinfo.specs_pot.barrier_scale;
    cs.use_barrier_rescaling() = potinfo.specs_pot.use_barrier_rescaling;
    cs.use_scales()    = potinfo.specs_pot.use_data_scales;

    cs.finalize_setup();







    cond_conv.functolabs = potinfo.specs_pot.functolabs;
    cond_conv.functolrel = potinfo.specs_pot.functolrel;
    cond_conv.gradtolabs = potinfo.specs_pot.gradtolabs;
    cond_conv.steptolabs = potinfo.specs_pot.steptolabs;
    cond_conv.steptolrel = potinfo.specs_pot.steptolrel;
    cond_conv.nitermin   = potinfo.specs_pot.nitermin;
    cond_conv.nitermax   = potinfo.specs_pot.nitermax;
    cond_conv.niterrestart = potinfo.specs_pot.niterrestart;

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







    // Vector<double> latcalc(ParamPot & parpot, Vector<CompoundStructureFit> & cmpvec);

    std::cout << std::endl;
    std::cout << "******************************************************" << std::endl;
    std::cout << "******************************************************" << std::endl;
    std::cout << "**                                                  **" << std::endl;
    std::cout << "**          Starting potential fitting ...          **" << std::endl;
    std::cout << "**                                                  **" << std::endl;
    std::cout << "******************************************************" << std::endl;
    std::cout << "******************************************************" << std::endl;
    std::cout << std::endl;
    std::cout << "INFO: Number of free fitting prameters: " << param.NXfree() << std::endl;
    std::cout << "INFO: Number of data points           : " << complistfit.NData() << std::endl;
    std::cout << std::endl;



    if (param.NXfree() == 0)
      aborterror("ERROR: There are no fitting parameters used, so fitting cannot be performed. Exiting.");
    

    int seed = potinfo.specs_pot.seed;

    /* ###############################################################################
       Fit function.
       ############################################################################### */
    if (potinfo.specs_pot.fitmet=="CG"){
      // Conjugate Gradients
      std::cout << "Using conjugate gradients method." << std::endl;
      ConjGrad< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="GN"){
      // Gauss-Newton
      std::cout << "Using Gauss-Newton method." << std::endl;
      GaussNewton< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="LM"){
      // Levenberg-Marquardt
      std::cout << "Using Levenberg-Marquardt method." << std::endl;
      LeveMarq< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="DL"){
      // Powell dog-leg
      std::cout << "Using Powell dog-leg method." << std::endl;
      PowellDogLeg< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 potinfo.specs_pot.dogleg_radius,
			 potinfo.specs_pot.dogleg_minradius,
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="SA"){
      // Simulated Annealing
      SimAnn< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 potinfo.specs_pot.simann_delta_rel,
			 seed,
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="MD"){
      // Molecular Dynamics
      std::cout << "Using Molecular Dynamics method." << std::endl;
      MolDynFit< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 potinfo.specs_pot.moldyn_min_dx,
			 potinfo.specs_pot.moldyn_max_dx,
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="PM"){
      // Powell's method
      std::cout << "Using Powell's method." << std::endl;
      Powell< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="SM"){
      // Simplex method
      std::cout << "Using Simplex method." << std::endl;
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
      DiffEvol< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 seed,
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="PS"){
      // Particle Swarm
      PartSwarm< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 seed,
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="BC"){
      // Bee colony
      BeeColony< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 seed,
			 cond_conv, cond_debug, cond_print);
    }
    else if (potinfo.specs_pot.fitmet=="GS"){
      // Gravitational Search
      GravSearch< ChiSqFunc<ParamPot, CompoundStructureFit, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
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
    //cout << "Made it here" << std::endl;
    //cout << "Xopt - cs.Param().X() is " << std::endl;
    //cout << Xopt - cs.Param().X() << std::endl;

    std::cout << std::endl;
    std::cout << "******************************************************" << std::endl;
    std::cout << "******************************************************" << std::endl;
    std::cout << "**                                                  **" << std::endl;
    std::cout << "**          Potential fitting terminated.           **" << std::endl;
    std::cout << "**                                                  **" << std::endl;
    std::cout << "******************************************************" << std::endl;
    std::cout << "******************************************************" << std::endl;
    std::cout << std::endl;


    // Calculate properties using optimized parameters:
    cs(Xopt);
    // Report:
    report_pot_prop( cs.Param(), cs.DataX(), cs.DataY(), cs.ModelDataY() );

    DX = cs.DataX();

  }
  // ***********************************************************************************
  // ***********************************************************************************
  // ***********************************************************************************
  // ***********************************************************************************





  std::cout << "Dumping potential parameters and compound properties to files ..." << std::endl;
  std::ofstream fout1, fout2;
  std::string fname1 = "report-potpar.dat";
  std::string fname2 = "report-compprop.dat";
  fout1.open(fname1.c_str());
  fout2.open(fname2.c_str());

  report_pot( &potinfo, true, true, fout1 );

  report_prop( DX, fout2 );

  fout1.close(); fout1.clear();
  fout2.close(); fout2.clear();




  double ts2 = omp_get_wtime();
  //clock_t clock2 = std::clock();
  
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "Total time (seconds) spent in program: "
       << ts2 - ts1
       << std::endl;
  std::cout << "**********************************************************************" << std::endl;




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

  std::cout << "####################################################################" << std::endl;
  std::cout << "Fitted parameters and their uncertainties:" << std::endl;
  std::cout << "Parameter    Value      Uncertainty    Rel. uncertainty:" << std::endl;
  std::cout << "####################################################################" << std::endl;
  for (i=0; i<NXf; i++){
    re = 0.0;
    if (abs(Xfopt[i])>eps) re = abs(dXfopt[i] / Xfopt[i]);
    std::cout << "x" << i << " " << Xfopt[i] << " " << dXfopt[i] << " " << re << std::endl;
  }
  std::cout << "####################################################################" << std::endl;
#endif
