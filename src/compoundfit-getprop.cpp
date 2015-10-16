
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
#include "utils-matrix-Choleskydecomp.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix-LUdecomp.hpp"
#include "utils-matrix-QRdecomp.hpp"
#include "utils-string.hpp"
#include "utils-vector.hpp"
#include "utils-errors.hpp"

#include "atomsystem.hpp"
//#include "compound.hpp"
#include "compoundfit.hpp"
#include "elem-iacs.hpp"
#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mtwister.hpp"
#include "param-pot.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "propfun.hpp"
#include "specs-fit-prop-pot.hpp"


#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;


using namespace utils;
using namespace constants;
using namespace physconst;
using namespace funcfit;
using boost::format;






void CompoundStructureFit::getprop(ParamPot & param){

  int N1,N2,N3;
  Vector3<double> b_ini(0), b_fin(0), boxlen_bak(0);
  Vector< Vector3<double> > bdir_ini(3, Vector3<double>(0)), bdir_fin(3, Vector3<double>(0));
  int i,j,k,p;
  Vector< Vector3<double> > pos_bak;
  MatrixSq3<double> boxdir_bak;
  MatrixSq3<double> alpha(0);
  Vector<double> X_opt;
  double E0, V0;




  MDSystem mds;

  // Set potential:
  mds.p_potinfo = param.p_potinfo;

  // Set other settings:
  mds.elem = param.p_potinfo->elem;
  mds.iacs = param.p_potinfo->iacs;


  mds.name = name;
    
  //cout << "Getting properties of compound " << name << " ..." << endl;


  //cout << "filename: " << cmpfit.filename << endl;


  // ###################################################################
  // MD relaxation of compound:
  // ###################################################################

  // Initialize common MD settings for all compounds:
  mds.specs_common = param.p_potinfo->specs_prop.mds_specs_common;
  // Initialize MD settings:
  mds.specs        = mds_specs;

  // Other MD settings from different places ...:
  mds.omp_info     = param.p_potinfo->omp_info;
  mds.rcut         = mds.rcut_max = mds.p_potinfo->get_rcut_max( elemnames );
  mds.skint        = mds.specs.skint;

  // These settings are applied to inherited base class members:
  mds.use_def_xyz_fmt = param.p_potinfo->specs_prop.mds_specs_common.use_def_dump_xyz_fmt;
  mds.def_xyz_fmt     = param.p_potinfo->specs_prop.mds_specs_common.def_dump_xyz_fmt;




  //cout << "rcut_max and skint are " << mds.rcut_max << " " << mds.skint << endl;


  // std::cout << "NB: compoundfit-getprop: Tstart " << mds.specs.Tstart << std::endl;


  bool retry;
  int mds_box_iter=0;
  int max_mds_box_iter = 10;
  while (true){

    std::cout << "Creating atom system ... " << std::endl;
    double rm = mds.rcut_max + mds.skint;
    mds.create_from_structure(*this, 2.0*rm); // removes old atoms

    /*
    if (mds_box_iter>0){
      // debug phase: dump box, and exit
      ofstream fdump;
      string dumpfn = "mds-retried-frame-" + mds.name + ".xyz";
      fdump.open(dumpfn.c_str());
      mds.dumpframe(fdump);
      fdump.close();
      //exit(EXIT_SUCCESS);
    }
    */


    b_ini[0] = mds.boxlen[0];
    b_ini[1] = mds.boxlen[1];
    b_ini[2] = mds.boxlen[2];
    bdir_ini[0] = mds.boxdir.col(0);
    bdir_ini[1] = mds.boxdir.col(1);
    bdir_ini[2] = mds.boxdir.col(2);
    for (i=0; i<mds.natoms(); i++){
      mds.type[i] = param.p_potinfo->elem.atomtype(mds.matter[i]);
    }


    std::cout << "Relaxing atom system ... " << std::endl;
    mds.relax();
    /*
    if (mds_box_iter>0){
      exit(EXIT_SUCCESS);
    }
    */


    E0 = mds.Ep_tot / mds.natoms();
    V0 = mds.vol / mds.natoms();

    retry = false;
    if (mds.N[0]<0 || mds.N[1]<0 || mds.N[2]<0) retry=true;
    if (mds.N[0]<0) mds.N[0] = -mds.N[0]+1;
    if (mds.N[1]<0) mds.N[1] = -mds.N[1]+1;
    if (mds.N[2]<0) mds.N[2] = -mds.N[2]+1;

    if (retry) mds_box_iter++;

    if (mds_box_iter > max_mds_box_iter)
      aborterror("ERROR: Too many (" + tostring(max_mds_box_iter)
		 + ") iterations of box construction "
		 + "when trying to relax compound " + name + ". Exiting.");

    if (! retry) break;
  }



  /* ############################################################################
     ############################################################################


     Predicted forces vs read-in forces


     ############################################################################
     ############################################################################
   */
  if (prop_use.frc){
    int nb = basis_elems.size();

    // Need to average over forces when compound is enlargened for mds !!!
    // Use: prop_pred.frc[ sitetype[i] ][k] += mds.frc[i][k] ...

    for (int i=0; i<nb; ++i){
      prop_pred.frc[i][0] = mds.frc[i][0];
      prop_pred.frc[i][1] = mds.frc[i][1];
      prop_pred.frc[i][2] = mds.frc[i][2];
    }
  }





  // ###################################################################
  // Get predicted properties
  // ###################################################################

  b_fin[0] = mds.boxlen[0];
  b_fin[1] = mds.boxlen[1];
  b_fin[2] = mds.boxlen[2];
  bdir_fin[0] = mds.boxdir.col(0);
  bdir_fin[1] = mds.boxdir.col(1);
  bdir_fin[2] = mds.boxdir.col(2);

  // ************************************************************************************
  // Due to ambiguity of how to specify a,b,c, etc in LAT files of read-in compounds
  // we just compare relative changes in box directions and apply them to user
  // specified values. If these specified values disagree with those in the LAT file
  // (due to outdated LAT file, e.g.), then there will be incorrect results reported.
  // ************************************************************************************
  if (prop_use.a) prop_pred.a = b_fin[0]/b_ini[0] * prop_readin.a;
  if (prop_use.b) prop_pred.b = b_fin[1]/b_ini[1] * prop_readin.b;
  if (prop_use.c) prop_pred.c = b_fin[2]/b_ini[2] * prop_readin.c;
  
  if (prop_use.bpa)
    prop_pred.bpa = (b_fin[1]/b_ini[1]) / (b_fin[0]/b_ini[0]) * prop_readin.bpa;
  if (prop_use.cpa)
    prop_pred.cpa = (b_fin[2]/b_ini[2]) / (b_fin[0]/b_ini[0]) * prop_readin.cpa;
  // ************************************************************************************
  // ************************************************************************************

  if (prop_use.r0){
    if (mds.natoms()!=2)
      aborterror("Error: No dimer in cell, cannot calculate r0 distance. Exiting.");
    Vector3<double> tv = mds.pos[0] - mds.pos[1];
    prop_pred.r0 = tv.magn();
  }
  
  if (prop_use.angle_ab)
    prop_pred.angle_ab = acos( (bdir_fin[0]        * bdir_fin[1]) /
				      (bdir_fin[0].magn() * bdir_fin[1].magn())
				      )/(2*PI) * 360.0;
  if (prop_use.angle_ac)
    prop_pred.angle_ac = acos( (bdir_fin[0]        * bdir_fin[2]) /
				      (bdir_fin[0].magn() * bdir_fin[2].magn())
				      )/(2*PI) * 360.0;
  if (prop_use.angle_bc)
    prop_pred.angle_bc = acos( (bdir_fin[1]        * bdir_fin[2]) /
				      (bdir_fin[1].magn() * bdir_fin[2].magn())
				      )/(2*PI) * 360.0;
    
  if (prop_use.Vatom)
    prop_pred.Vatom = mds.vol/mds.natoms();

      
  if (prop_use.Ecoh)
    prop_pred.Ecoh = mds.Ep_tot/mds.natoms();

  if (prop_use.Ecoh_delta)
    prop_pred.Ecoh_delta = mds.Ep_tot/mds.natoms();


  if (prop_use.Emix){
    prop_pred.Emix = mds.Ep_tot;
    //cout << "Total potential energy of compound is " << mds.Ep_tot << endl;

    // Number of different elements present:
    Vector<int> npres(mds.elem.nelem(), 0);
    for (j=0; j<mds.natoms(); ++j){
      int k = param.p_potinfo->elem.name2idx( mds.matter[j] );
      ++(npres[k]);
    }
    for (j=0; j<mds.elem.nelem(); ++j){
      std::string nname = param.p_potinfo->elem.idx2name(j);
      cout << "Found " << npres[j] << " atoms with element index " << j
	   << " with name " << nname
	   << " and MD type " << param.p_potinfo->elem.atomtype( nname )
	   << " each with cohesive energy " << param.p_potinfo->Ecoh_ref[ j ] << endl;
      
      prop_pred.Emix -= npres[j] * param.p_potinfo->Ecoh_ref[ j ];
    }
    prop_pred.Emix /= mds.natoms();
  }

  prop_pred.Fmax = mds.F_max;
  prop_pred.Pmax = mds.P_max;
  prop_pred.displmax = mds.displ_max;



  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Make sure additional calculations of properties
  // start with the relaxed system
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pos_bak    = mds.pos;
  boxdir_bak = mds.boxdir;
  boxlen_bak = mds.boxlen;
  V0         = mds.vol/mds.natoms();
  E0         = mds.Ep_tot/mds.natoms();



  // ############################################################################
  // Bulk modulus and its pressure derivative
  // ############################################################################
  if (prop_use.B || prop_use.Bp){
    std::cout << "Calculating for B and B'(P) ..." << std::endl;
    get_B_Bp(mds, param, pos_bak, boxdir_bak, boxlen_bak, E0, V0);
  }


  // ############################################################################
  // Elastic constants
  // ############################################################################
  j=0;
  for (k=0; k<6; ++k){
    for (p=0; p<6; ++p){
      if (prop_use.C.elem(k,p)) j++;
    }
  }
  if (j>0){
    std::cout << "Calculating for Cij ..." << flush << std::endl;
    get_Cij(mds, param, pos_bak, boxdir_bak, boxlen_bak, E0, V0);
  }

    




}






void CompoundStructureFit::get_B_Bp(MDSystem             & mds,
				    ParamPot             & param,
				    Vector< Vector3<double> > & pos_bak,
				    MatrixSq3<double>           & boxdir_bak,
				    Vector3<double>           & boxlen_bak,
				    double & E0,
				    double & V0
				    ){

  int j,k,p,Nf;
  MatrixSq3<double> alpha(0);
  double fmin, fmax, ef, f, df, g, td;
  Vector<double> xp, yp, dyp;
  Vector<double> Xopt;
  double B0, Bp0, dV;
  double d2E, d3E;
  bool fit_OK;
  double min_V=V0, min_E=E0, guess_B, guess_Bp;
  double eps = std::numeric_limits<double>::epsilon();
  double small = sqrt(eps);
  bool rel_sys;

  rel_sys = param.p_potinfo->specs_prop.BM_rel_sys;
  fmin    = param.p_potinfo->specs_prop.BM_fmin;
  fmax    = param.p_potinfo->specs_prop.BM_fmax;
  Nf      = param.p_potinfo->specs_prop.BM_Nf;
  ef      = param.p_potinfo->specs_prop.BM_ef;
  if (Nf<3) Nf=3;
  df      = (fmax - fmin)/Nf;
  xp.resize(Nf); yp.resize(Nf); dyp.resize(Nf);


  // Loop over particular symmetry changes
  // (only one => no loop)

  // Loop over differential changes
  for (j=0; j<Nf; ++j){
    f = fmin + j*df;
    g = 1 + f;
    xp[j] = V0 * g*g*g;
    for (k=0; k<3; ++k) for (p=0; p<3; ++p) alpha.elem(k,p)=0;
    alpha.elem(0,0) = alpha.elem(1,1) = alpha.elem(2,2) = g; // Note!!!
    mds.transform_cell(alpha);


    // **************************************************************
    // **************************************************************
    double Epa=0;
    if (! rel_sys){
      // +++++++++++++++++++++++++++++++++++++++++++++++++
      // Option A: Do not allow internal relaxation:
      // +++++++++++++++++++++++++++++++++++++++++++++++++

      mds.get_all_neighborcollections();
      Epa = mds.calc_potential_energy() / mds.natoms();
    }
    else {
      // +++++++++++++++++++++++++++++++++++++++++++++++++
      // Option B: Allow internal relaxation:
      // +++++++++++++++++++++++++++++++++++++++++++++++++

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // Update specs for this relaxation:
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MDSettings specs_bak = mds.specs;
      mds.specs.use_Pcontrol  = false;
      mds.specs.quench_always = true;

      // Relax using current system. Only enlargen it if position scaling
      // has made the system too small.
      bool retry;
      while (true){
	mds.relax();

	Epa = mds.Ep_tot / mds.natoms();

	retry = false;
	if (mds.N[0]<0 || mds.N[1]<0 || mds.N[2]<0) retry=true;
	if (mds.N[0]<0) mds.N[0] = -mds.N[0]+1;
	if (mds.N[1]<0) mds.N[1] = -mds.N[1]+1;
	if (mds.N[2]<0) mds.N[2] = -mds.N[2]+1;

	if (! retry) break;

	double rm = mds.rcut_max + mds.skint;
	mds.create_from_structure(*this, 2.0*rm); // removes old atoms
	for (int i=0; i<mds.natoms(); i++){
	  mds.type[i] = param.p_potinfo->elem.atomtype(mds.matter[i]);
	}
      }
      mds.specs = specs_bak;
    }
    // **************************************************************
    // **************************************************************

    yp[j] = Epa;
    dyp[j] = ((Epa<0)? -Epa: Epa) * ef;

    {
      std::ofstream fdump;
      std::string dumpfn = "mds-frame-B-calc-" + mds.name + "-step-" + tostring(j) + ".xyz";
      fdump.open(dumpfn.c_str());
      mds.dumpframe(fdump);
      fdump.close();
    }

    // Reset:
    mds.pos    = pos_bak;
    mds.boxdir = boxdir_bak;
    mds.boxlen = boxlen_bak;


    /*
    cout << "f g V Epa "
	 << f << " "
	 << g << " "
	 << xp[j] << " "
	 << yp[j] << " "
	 << dyp[j] << endl;
    */
    /*
    cout << xp[j] << " "
	 << yp[j] << " "
	 << dyp[j] << endl;
    */


  }


  {
    std::ofstream fdump;
    std::string dumpfn = "data-V-Epa-dEpa-" + name + ".dat";
    fdump.open(dumpfn.c_str());
    for (j=0; j<Nf; ++j)
      fdump << format("%20.10e  %20.10e  %20.10e") % xp[j] % yp[j] % dyp[j] << std::endl;
    fdump.close();
  }




  // Fit:
  Param pp;
  ChiSqFunc<Param, double, double> cs;
  Cond_Conv  cond_conv;
  Cond_Debug cond_debug;
  Cond_Print cond_print;


  // General physically motivated guess:
  guess_B  = 100.0 * GPa_to_eVA3; // unit now: eV/Ang^3
  guess_Bp = 4.0;

  pp.init(4);
  pp.X(0) = min_E;
  pp.X(1) = min_V;
  pp.X(2) = guess_B;
  pp.X(3) = guess_Bp;


  cs.Param() = pp;
  cs.DataX() = xp;
  cs.DataY() = yp;
  cs.DataUncertaintyY() = dyp;
  cs.ModelFuncPointer() = fun_bmeos;
  cs.finalize_setup();

  fit_OK = false;



  cs.barrier_scale() = param.p_potinfo->specs_prop.barrier_scale;
  cs.use_scales()    = param.p_potinfo->specs_prop.use_data_scales;

  cond_conv.functolabs = param.p_potinfo->specs_prop.functolabs;
  cond_conv.functolrel = param.p_potinfo->specs_prop.functolrel;
  cond_conv.gradtolabs = param.p_potinfo->specs_prop.gradtolabs;
  cond_conv.steptolabs = param.p_potinfo->specs_prop.steptolabs;
  cond_conv.steptolrel = param.p_potinfo->specs_prop.steptolrel;
  cond_conv.nitermin   = param.p_potinfo->specs_prop.nitermin;
  cond_conv.niterrestart = param.p_potinfo->specs_prop.niterrestart;

  cond_conv.report_conv = param.p_potinfo->specs_prop.report_conv;
  cond_conv.prefix_report_conv = "propfit B,B' conv: ";

  cond_debug.debug_fit_level0 = param.p_potinfo->specs_prop.debug_fit_level0;
  cond_debug.debug_fit_level1 = param.p_potinfo->specs_prop.debug_fit_level1;
  cond_debug.debug_fit_level2 = param.p_potinfo->specs_prop.debug_fit_level2;
  cond_debug.debug_fit_level3 = param.p_potinfo->specs_prop.debug_fit_level3;
  cond_debug.debug_fit_level4 = param.p_potinfo->specs_prop.debug_fit_level4;
  cond_debug.prefix_debug_fit_level0 = "propfit B,B' debug0: ";
  cond_debug.prefix_debug_fit_level1 = "propfit B,B' debug1: ";
  cond_debug.prefix_debug_fit_level2 = "propfit B,B' debug2: ";
  cond_debug.prefix_debug_fit_level3 = "propfit B,B' debug3: ";
  cond_debug.prefix_debug_fit_level4 = "propfit B,B' debug4: ";

  cond_print.report_iter = param.p_potinfo->specs_prop.report_iter;
  cond_print.report_warn = param.p_potinfo->specs_prop.report_warn;
  cond_print.report_error= param.p_potinfo->specs_prop.report_error;
  cond_print.prefix_report_iter  = "propfit B,B' iter: ";
  cond_print.prefix_report_warn  = "propfit B,B' warn: ";
  cond_print.prefix_report_error = "propfit B,B' error: ";

  if (cond_debug.debug_fit_level0){
    cond_conv.report_conv  = true;
    cond_print.report_iter = true;
    cond_print.report_warn = true;
    cond_print.report_error= true;
    cs.debug();
  }




  int seed = param.p_potinfo->specs_prop.seed;

  if (param.p_potinfo->specs_prop.fitmet=="CG"){
    // Conjugate Gradients
    ConjGrad< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="GN"){
    // Gauss-Newton
    GaussNewton< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="LM"){
    // Levenberg-Marquardt
    LeveMarq< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="DL"){
    // Powell dog-leg
    PowellDogLeg< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       param.p_potinfo->specs_prop.dogleg_radius,
		       param.p_potinfo->specs_prop.dogleg_minradius,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="SA"){
    // Simulated Annealing
    SimAnn< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       param.p_potinfo->specs_prop.simann_delta_rel,
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="MD"){
    // Molecular Dynamics
    //std::cout << "Using Molecular Dynamics method." << std::endl;
    MolDynFit< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       param.p_potinfo->specs_prop.moldyn_min_dx,
		       param.p_potinfo->specs_prop.moldyn_max_dx,
		       cond_conv, cond_debug, cond_print);
  }




  else if (param.p_potinfo->specs_prop.fitmet=="PM"){
    // Powell's method
    Powell< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="SM"){
    // Simplex method
    SimplexFit< ChiSqFunc<Param, double, double> > fm(cs);
    /*
    Vector<double> X_displ( cs.Param().X().size(), 
			    param.p_potinfo->specs_prop.simplex_delta );
    */
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="DE"){
    // Differential evolution
    DiffEvol< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="PS"){
    // Particle Swarm
    PartSwarm< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="BC"){
    // Bee colony
    BeeColony< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="GS"){
    // Gravitational Search
    GravSearch< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else
    aborterror("Error: Unknown fitting method " + param.p_potinfo->specs_prop.fitmet + ". " +
	       "Exiting.");

  /*
  cout << "Xopt: " << Xopt << endl;
  cout << "Xmax: " << cs.Param().Xmin() << endl;
  cout << "Xmax: " << cs.Param().Xmax() << endl;

  cs(Xopt);
  cout << " ************************ B, B': Done with recalc after Xopt obtained." << flush << endl;
  */



  // If fitting went alright, then we are almost done. Else we need to get
  // a numerical estimate of the required derivatives.

  if (fit_OK){
    E0  = Xopt[0];
    V0  = Xopt[1];
    B0  = Xopt[2] * eVA3_to_GPa; // unit now: GPa
    Bp0 = Xopt[3];
  }
  else {
    // **************************************************************
    // Numerical derivative
    // **************************************************************
    E0  = min_E;
    V0  = min_V;
    
    double eps = std::numeric_limits<double>::epsilon();
    int j, nset, nmid;

    nset = 6;
    Vector<double> Ep(nset + 1, 0);
    nmid = nset/2;   // e.g. 3

    df = pow(eps, 1.0/3.0);
    g = 1 + df;
    dV = V0*g*g*g - V0;

    for (j=0; j<=nset; ++j){
      if (j==nmid){
	Ep[j]=E0;
	continue;
      }
      g = 1 + (j - nmid) * df;
      for (k=0; k<3; ++k) for (p=0; p<3; ++p) alpha.elem(k,p)=0;
      alpha.elem(0,0) = alpha.elem(1,1) = alpha.elem(2,2) = g; // Note!!!
      mds.transform_cell(alpha);
      mds.get_all_neighborcollections();
      Ep[j] = mds.calc_potential_energy() / mds.natoms();
    }

    d2E = Ep[nmid+1] - 2*Ep[nmid] + Ep[nmid-1];
    guess_B = B0 = - d2E/(dV*dV) * V0;
    //B0  = d2E/(4.0*dV*dV) * V0 * eVA3_to_GPa;
    d3E = Ep[nmid+2] - 2*Ep[nmid+1] + 2*Ep[nmid-1] - Ep[nmid-2];
    Bp0 = -1.0 - V0*V0/B0 * d3E/(2.0*dV*dV*dV);
    guess_Bp = Bp0;

    B0  = guess_B * eVA3_to_GPa;
    Bp0 = guess_Bp;
    std::cout << "Warning: Failed to obtain bulk modulus and/or its pressure derivative for compound " + name + ". " +
      "Using numerical estimates " << B0 << " GPa and " << Bp0 << std::endl;
    // **************************************************************
    // **************************************************************
  }



  if (prop_use.Ecoh)  prop_pred.Ecoh  = E0;
  if (prop_use.Vatom) prop_pred.Vatom = V0;
  if (prop_use.B)     prop_pred.B     = B0;
  if (prop_use.Bp)    prop_pred.Bp    = Bp0;

  // cout << "Done with fitting B, B'" << flush << endl;

  return;
}






void CompoundStructureFit::get_Cij(MDSystem             & mds,
				   ParamPot             & param,
				   Vector< Vector3<double> > & pos_bak,
				   MatrixSq3<double>           & boxdir_bak,
				   Vector3<double>           & boxlen_bak,
				   double & E0,
				   double & V0
				   ){
		
  bool quit, fit_OK;
  int j,k,p,Nf,isym;
  Vector<double> xp, yp, dyp, Xopt;
  int NC=21;
  Matrix<double> Clincomb(7,7,0);
  double epsi[7];
  MatrixSq3<double> alpha(0);
  Matrix<double> C(6,6,0);
  double fmin, fmax, ef, df, f, td;
  double C11, C12, C13, C14, C15, C16;
  double C22, C23, C24, C25, C26;
  double C33, C34, C35, C36;
  double C44, C45, C46;
  double C55, C56;
  double C66;
  double epsi1, epsi2, epsi3, epsi4, epsi5, epsi6;

  double eps = std::numeric_limits<double>::epsilon(), min_f=0.0, min_E=E0, guess_C;
  double small = sqrt(eps);
  bool rel_sys;

  rel_sys = param.p_potinfo->specs_prop.C_rel_sys;
  fmin    = param.p_potinfo->specs_prop.C_fmin;
  fmax    = param.p_potinfo->specs_prop.C_fmax;
  Nf      = param.p_potinfo->specs_prop.C_Nf;
  ef      = param.p_potinfo->specs_prop.C_ef;
  if (Nf<3) Nf=3;
  df      = (fmax - fmin)/Nf;
  xp.resize(Nf); yp.resize(Nf); dyp.resize(Nf);


  C11 = C12 = C13 = C14 = C15 = C16 = 0.0;
  C22 = C23 = C24 = C25 = C26 = 0.0;
  C33 = C34 = C35 = C36 = 0.0;
  C44 = C45 = C46 = 0.0;
  C55 = C56 = 0.0;
  C66 = 0.0;


  /* i:
    0  1  2  3  4  5
       6  7  8  9  10
          11 12 13 14
	     15 16 17
	        18 19
		   20  */

  
  Matrix<bool> Cuse(7,7,false);
  get_Cuse(Cuse);

  // Diagonal elements needed for many non-diagonal elements.
  // These may not have been specified by user, so we must force
  // usage of them, without messing out the ''officially used''
  // elements.
  Matrix<bool> Cuse_hidden(7,7,false);
  for (int i=1; i<=6; ++i)
    Cuse_hidden.elem(i,i) = true;



  /* ###############################################################################
     Find the elastic constants.
       
     For a cubic lattice we have the elastic constants
     C11, C12, C44
       
     For a hexagonal lattice we have the elastic constants
     C11, C12, C44, C13, C33
     ############################################################################### */

  /*
           eps1  eps6  eps5
    eps =  eps6  eps2  eps4
           eps5  eps4  eps3

            1  0  0
    alpha = 0  1  0 + eps
            0  0  1

    E - E0 = 0.5 * V0 * (
          C11 * eps1^2
    +     C22 * eps2^2
    +     C33 * eps3^2
    + 4 * C44 * eps4^2
    + 4 * C55 * eps5^2
    + 4 * C66 * eps6^2
    + 2 * C12 * eps1 * eps2
    + 2 * C13 * eps1 * eps3
    + 2 * C23 * eps2 * eps3
    + 4 * eps1 * ( C14 * eps4 + C15 * eps5 + C16 * eps6 )
    + 4 * eps2 * ( C24 * eps4 + C25 * eps5 + C26 * eps6 )
    + 4 * eps3 * ( C34 * eps4 + C35 * eps5 + C36 * eps6 )
    + 8 * eps4 * ( C45 * eps5 + C46 * eps6 )
    + 8 * eps5 *   C56 * eps6
    )

   */

  quit = false;


  // Loop over particular symmetry changes
  for (int ik=1; ik<=6; ++ik){
    for (int ip=1; ip<=6; ++ip){

      int choice_sum=0;
      if (Cuse.elem(ik,ip)) choice_sum++;
      if (ik==ip && Cuse_hidden.elem(ik,ip)) choice_sum++;

      if (choice_sum==0) continue;


      // Loop over differential changes
      for (j=0; j<Nf; ++j){
	f = fmin + j*df;

	/* Establish the transformation. */
	for (k=0; k<3; ++k) for (p=0; p<3; ++p) alpha.elem(k,p) = 0.0;
	for (k=0; k<3; ++k) alpha.elem(k,k) = 1.0;
	for (k=0; k<=6; ++k) epsi[k]=0.0;
	epsi1 = epsi2 = epsi3 = epsi4 = epsi5 = epsi6 = 0.0;

	epsi[ik] = f;
	epsi[ip] = f;

	epsi1 = epsi[1];
	epsi2 = epsi[2];
	epsi3 = epsi[3];
	epsi4 = epsi[4];
	epsi5 = epsi[5];
	epsi6 = epsi[6];

    
	xp[j] = f;
	//xp[j] = V0 * g*g*g;
	//for (k=0; k<3; ++k) for (p=0; p<3; ++p) alpha.elem(k,p)=0;
	//alpha.elem(0,0) = alpha.elem(1,1) = alpha.elem(2,2) = g; // Note!!!


	alpha.elem(0,0) += epsi1;
	alpha.elem(0,1) += epsi6;
	alpha.elem(0,2) += epsi5;

	alpha.elem(1,0) += epsi6;
	alpha.elem(1,1) += epsi2;
	alpha.elem(1,2) += epsi4;

	alpha.elem(2,0) += epsi5;
	alpha.elem(2,1) += epsi4;
	alpha.elem(2,2) += epsi3;


	mds.transform_cell(alpha);


	// **************************************************************
	// **************************************************************
	double Epa=0;
	if (! rel_sys){
	  // +++++++++++++++++++++++++++++++++++++++++++++++++
	  // Option A: Do not allow internal relaxation:
	  // +++++++++++++++++++++++++++++++++++++++++++++++++
	
	  mds.get_all_neighborcollections();
	  Epa = mds.calc_potential_energy() / mds.natoms();
	}
	else {
	  // +++++++++++++++++++++++++++++++++++++++++++++++++
	  // Option B: Allow internal relaxation:
	  // +++++++++++++++++++++++++++++++++++++++++++++++++

	  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  // Update specs for this relaxation:
	  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  MDSettings specs_bak = mds.specs;
	  mds.specs.use_Pcontrol  = false;
	  mds.specs.quench_always = true;

	  // Relax using current system. Only enlargen it if position scaling
	  // has made the system too small.
	  bool retry;
	  while (true){
	    mds.relax();

	    Epa = mds.Ep_tot / mds.natoms();

	    retry = false;
	    if (mds.N[0]<0 || mds.N[1]<0 || mds.N[2]<0) retry=true;
	    if (mds.N[0]<0) mds.N[0] = -mds.N[0]+1;
	    if (mds.N[1]<0) mds.N[1] = -mds.N[1]+1;
	    if (mds.N[2]<0) mds.N[2] = -mds.N[2]+1;

	    if (! retry) break;

	    double rm = mds.rcut_max + mds.skint;
	    mds.create_from_structure(*this, 2.0*rm); // removes old atoms
	    for (int i=0; i<mds.natoms(); i++){
	      mds.type[i] = param.p_potinfo->elem.atomtype(mds.matter[i]);
	    }
	  }
	  mds.specs = specs_bak;
	}
	// **************************************************************
	// **************************************************************


	yp[j] = Epa;
	dyp[j] = ((Epa<0)? -Epa: Epa) * ef;


	{
	  std::ofstream fdump;
	  std::string dumpfn = "mds-frame-C-calc-" + mds.name + "-isym-"
	    + tostring(ik) + tostring(ip)
	    + "-step-" + tostring(j) + ".xyz";
	  fdump.open(dumpfn.c_str());
	  mds.dumpframe(fdump);
	  fdump.close();
	}


	// Reset:
	mds.pos    = pos_bak;
	mds.boxdir = boxdir_bak;
	mds.boxlen = boxlen_bak;

      }



      {
	std::ofstream fdump;
	std::string dumpfn = "data-f-Epa-dEpa-" + name + "-" + tostring(ik) + tostring(ip) + ".dat";
	fdump.open(dumpfn.c_str());
	for (j=0; j<Nf; ++j)
	  fdump << format("%20.10e  %20.10e  %20.10e") % xp[j] % yp[j] % dyp[j] << std::endl;
	fdump.close();
      }








      /* Fit the energies 'Evec' as a function of 'fvec' to a 2nd degree polynomial.
	 E(f) = c_0 + c_1 * f^2
      */
    
      // Fit:
      Param pp;
      ChiSqFunc<Param, double, double> cs;
      Cond_Conv  cond_conv;
      Cond_Debug cond_debug;
      Cond_Print cond_print;


      // General physically motivated guess:
      guess_C = 0.5 * V0 * (100.0 * GPa_to_eVA3); // unit now: eV/Ang^3

      pp.init(3);
      pp.X(0) = min_f;
      pp.X(1) = min_E;
      pp.X(2) = guess_C;

   
      cs.Param() = pp;
      cs.DataX() = xp;
      cs.DataY() = yp;
      cs.DataUncertaintyY() = dyp;
      cs.ModelFuncPointer() = fun_poly2;
      cs.finalize_setup();


      fit_OK = false;

      cs.barrier_scale() = param.p_potinfo->specs_prop.barrier_scale;
      cs.use_scales()    = param.p_potinfo->specs_prop.use_data_scales;

      cond_conv.functolabs = param.p_potinfo->specs_prop.functolabs;
      cond_conv.functolrel = param.p_potinfo->specs_prop.functolrel;
      cond_conv.gradtolabs = param.p_potinfo->specs_prop.gradtolabs;
      cond_conv.steptolabs = param.p_potinfo->specs_prop.steptolabs;
      cond_conv.steptolrel = param.p_potinfo->specs_prop.steptolrel;
      cond_conv.nitermin   = param.p_potinfo->specs_prop.nitermin;
      cond_conv.nitermax   = param.p_potinfo->specs_prop.nitermax;
      cond_conv.niterrestart = param.p_potinfo->specs_prop.niterrestart;

      cond_conv.report_conv = param.p_potinfo->specs_prop.report_conv;
      cond_conv.prefix_report_conv = "propfit Cij conv: ";

      cond_debug.debug_fit_level0 = param.p_potinfo->specs_prop.debug_fit_level0;
      cond_debug.debug_fit_level1 = param.p_potinfo->specs_prop.debug_fit_level1;
      cond_debug.debug_fit_level2 = param.p_potinfo->specs_prop.debug_fit_level2;
      cond_debug.debug_fit_level3 = param.p_potinfo->specs_prop.debug_fit_level3;
      cond_debug.debug_fit_level4 = param.p_potinfo->specs_prop.debug_fit_level4;
      cond_debug.prefix_debug_fit_level0 = "propfit Cij debug0: ";
      cond_debug.prefix_debug_fit_level1 = "propfit Cij debug1: ";
      cond_debug.prefix_debug_fit_level2 = "propfit Cij debug2: ";
      cond_debug.prefix_debug_fit_level3 = "propfit Cij debug3: ";
      cond_debug.prefix_debug_fit_level4 = "propfit Cij debug4: ";

      cond_print.report_iter = param.p_potinfo->specs_prop.report_iter;
      cond_print.report_warn = param.p_potinfo->specs_prop.report_warn;
      cond_print.report_error= param.p_potinfo->specs_prop.report_error;
      cond_print.prefix_report_iter  = "propfit Cij iter: ";
      cond_print.prefix_report_warn  = "propfit Cij warn: ";
      cond_print.prefix_report_error = "propfit Cij error: ";


      if (cond_debug.debug_fit_level0){
	cond_conv.report_conv  = true;
	cond_print.report_iter = true;
	cond_print.report_warn = true;
	cond_print.report_error= true;
	cs.debug();
      }


      int seed = param.p_potinfo->specs_prop.seed;

      if (param.p_potinfo->specs_prop.fitmet=="CG"){
	// Conjugate Gradients
	ConjGrad< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="GN"){
	// Gauss-Newton
	GaussNewton< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="LM"){
	// Levenberg-Marquardt
	LeveMarq< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="DL"){
	// Powell dog-leg
	PowellDogLeg< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   param.p_potinfo->specs_prop.dogleg_radius,
			   param.p_potinfo->specs_prop.dogleg_minradius,
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="SA"){
	// Simulated Annealing
	SimAnn< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   param.p_potinfo->specs_prop.simann_delta_rel,
			   seed,
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="MD"){
	// Molecular Dynamics
	//std::cout << "Using Molecular Dynamics method." << std::endl;
	MolDynFit< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   param.p_potinfo->specs_prop.moldyn_min_dx,
			   param.p_potinfo->specs_prop.moldyn_max_dx,
			   cond_conv, cond_debug, cond_print);
      }



      else if (param.p_potinfo->specs_prop.fitmet=="PM"){
	// Powell's method
	Powell< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="SM"){
	// Simplex method
	SimplexFit< ChiSqFunc<Param, double, double> > fm(cs);
	/*
	  Vector<double> X_displ( cs.Param().X().size(), 
	  param.p_potinfo->specs_prop.simplex_delta );
	*/
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   seed,
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="DE"){
	// Differential evolution
	DiffEvol< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   seed,
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="PS"){
	// Particle Swarm
	PartSwarm< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   seed,
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="BC"){
	// Bee colony
	BeeColony< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   seed,
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else if (param.p_potinfo->specs_prop.fitmet=="GS"){
	// Gravitational Search
	GravSearch< ChiSqFunc<Param, double, double> > fm(cs);
	Xopt = fm.minimize(cs.Param().X(),
			   cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			   seed,
			   cond_conv, cond_debug, cond_print);
	fit_OK = fm.status.fit_OK;
      }
      else
	aborterror("Error: Unknown fitting method " + param.p_potinfo->specs_prop.fitmet + ". " +
		   "Exiting.");



      /*
	cout << "Xopt: " << Xopt << endl;
	cout << "Xmax: " << cs.Param().Xmin() << endl;
	cout << "Xmax: " << cs.Param().Xmax() << endl;

	cs(Xopt);
	cout << " ************************ Cij: Done with recalc after Xopt obtained." << flush << endl;
      */





      // If fitting went alright, then we are almost done. Else we need to get
      // a numerical estimate of the required derivatives.

      if (fit_OK){
	E0  = Xopt[1];
	//f0  = Xopt[1];
	Clincomb.elem(ik,ip) = Xopt[2] / V0 * eVA3_to_GPa; // unit now: GPa
      }
      else {
	// **************************************************************
	// Numerical estimate
	// **************************************************************
	int min_idx=-1;
	for (k=0; k<yp.size(); ++k){
	  if (k==0 || (k>0 && yp[k]<td)){
	    td = yp[k]; min_idx = k;
	  }
	}
	Vector3<double> vx(0), vy(0);
	int im=0;
	if ((min_idx-1)>=0 && (min_idx+1)<yp.size()) im = min_idx - 1;
	else if (min_idx==0)                         im = min_idx;
	else if (min_idx==yp.size()-1)               im = min_idx - 2;
	vx[0] = xp[im];
	vx[1] = xp[im+1];
	vx[2] = xp[im+2];
	vy[0] = yp[im];
	vy[1] = yp[im+1];
	vy[2] = yp[im+2];
	get_parabolic_fit_from_triplet(vx, vy, min_E, guess_C, min_f, cond_debug.debug_fit_level1);
	E0 = min_E;
	//f0 = min_f;
	Clincomb.elem(ik,ip) = guess_C / V0 * eVA3_to_GPa; // unit now: GPa
	std::cout << "Warning: Failed to obtain combination of linear constants for compound " + name + ". " +
	  "Using estimate " << guess_C << " GPa." << std::endl;
      }







    }
  }
  // End of loop over isym


  // Get all Cij:
  get_Cresolved(Clincomb, C);


  // Get the used Cij only:
  for (k=0; k<6; ++k){
    for (p=0; p<6; ++p){
      if (prop_use.C.elem(k,p))
	prop_pred.C.elem(k,p) = C.elem(k,p);
    }
  }


}




