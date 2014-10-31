
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
#include "utils-matrix-Choleskydecomp.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix-LUdecomp.hpp"
#include "utils-matrix-QRdecomp.hpp"
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
#include "propfun.hpp"
#include "specs-fit-prop-pot.hpp"

using namespace std;
using namespace utils;
using namespace exiterrors;
using namespace constants;
using namespace physconst;
using namespace funcfit;
using boost::format;



void CompoundStructureFit::getprop(ParamPot & param){

  int N1,N2,N3;
  Vector<double> b_ini(3,0), b_fin(3,0), boxlen_bak(3,0);
  Vector< Vector<double> > bdir_ini(3, Vector<double>(3,0)), bdir_fin(3, Vector<double>(3,0));
  int i,j,k,p;
  Vector< Vector<double> > pos_bak;
  Matrix<double> boxdir_bak;
  Matrix<double> alpha(3,3,0);
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

  // Initialize common MD settings:
  mds.specs_common = param.p_potinfo->specs_prop.mds_specs_common;
    
  mds.omp_info = param.p_potinfo->omp_info;


  // Initialize MD settings:
  mds.specs    = mds_specs;

  // Other settings:
  mds.rcut = mds.rcut_max = mds.p_potinfo->get_rcut_max( elemnames );
  mds.skint = mds.specs.skint;



  double rm = mds.rcut_max + mds.skint;
  N1 = N2 = N3 = 0;
  //cout << "rcut_max and skint are " << mds.rcut_max << " " << mds.skint << endl;



  bool retry;
  while (true){
    // Set atoms:
    cout << "Creating atom system ... " << endl;
    mds.create_from_structure(*this, N1, N2, N3, 2.0*rm); // removes old atoms

    b_ini[0] = mds.boxlen[0];
    b_ini[1] = mds.boxlen[1];
    b_ini[2] = mds.boxlen[2];
    bdir_ini[0] = mds.boxdir.col(0);
    bdir_ini[1] = mds.boxdir.col(1);
    bdir_ini[2] = mds.boxdir.col(2);
    for (i=0; i<mds.natoms(); i++){
      mds.type[i] = param.p_potinfo->elem.atomtype(mds.matter[i]);
    }


      
    // Perform MD relaxation:
    cout << "Relaxing atom system ... " << endl;
    if (mds.specs.is_ref_comp)
      mds.relax(N1, N2, N3, false); // never quick mode on a reference compound
    else
      mds.relax(N1, N2, N3, mds.specs_common.quick_mode);


    E0 = mds.Ep_tot / mds.natoms();
    V0 = mds.V / mds.natoms();
    retry = false;
    if (N1<0 || N2<0 || N3<0) retry=true;
    if (N1<0) N1 = -N1+1;
    if (N2<0) N2 = -N2+1;
    if (N3<0) N3 = -N3+1;
    
    if (! retry) break;
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
    Vector<double> tv = mds.pos[0] - mds.pos[1];
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
    prop_pred.Vatom = mds.V/mds.natoms();
      
  if (prop_use.Ecoh)
    prop_pred.Ecoh = mds.Ep_tot/mds.natoms();

  if (prop_use.Emix){
    prop_pred.Emix = mds.Ep_tot;
    //cout << "Total potential energy of compound is " << mds.Ep_tot << endl;

    Vector<int> ntype(mds.elem.nelem(), 0);
    for (j=0; j<mds.natoms(); ++j)
      ++(ntype[ (mds.type[j]<0) ? -mds.type[j] : mds.type[j] ]);
    for (j=0; j<mds.elem.nelem(); ++j){
      /*
	cout << "Found " << ntype[j] << " atoms of type " << j
	<< " each with cohesive energy " << param.p_potinfo->Ecoh_ref[ j ] << endl;
      */
      prop_pred.Emix -= ntype[j] * param.p_potinfo->Ecoh_ref[ j ];
    }
    prop_pred.Emix /= mds.natoms();
  }


  pos_bak    = mds.pos;
  boxdir_bak = mds.boxdir;
  boxlen_bak = mds.boxlen;
  V0         = mds.V/mds.natoms();
  E0         = mds.Ep_tot/mds.natoms();


  // ############################################################################
  // Bulk modulus and its pressure derivative
  // ############################################################################
  if (prop_use.B || prop_use.Bp){
    cout << "Calculating for B and B'(P) ..." << endl;
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
    cout << "Calculating for Cij ..." << endl;
    get_Cij(mds, param, pos_bak, boxdir_bak, boxlen_bak, E0, V0);
  }

    
  if (prop_use.Fmax)     prop_pred.Fmax = mds.F_max;
  if (prop_use.Pmax)     prop_pred.Pmax = mds.P_max;
  if (prop_use.displmax) prop_pred.displmax = mds.displ_max;




}






void CompoundStructureFit::get_B_Bp(MDSystem             & mds,
				    ParamPot             & param,
				    Vector< Vector<double> > & pos_bak,
				    Matrix<double>           & boxdir_bak,
				    Vector<double>           & boxlen_bak,
				    double & E0,
				    double & V0
				    ){
		
  int j,k,p,Nf;
  Matrix<double> alpha(3,3,0);
  double fmin, fmax, df, f, g;
  double td;
  Vector<double> xp, yp, dyp;
  Vector<double> Xopt;
  double B0, Bp0, dV;
  double d2E, d3E;
  bool fit_OK;
  double small = sqrt( std::numeric_limits<double>::epsilon() );
  double min_V=V0, min_E=E0, guess_B, guess_Bp;

  fmin = param.p_potinfo->specs_prop.BM_fmin;
  fmax = param.p_potinfo->specs_prop.BM_fmax;
  Nf   = param.p_potinfo->specs_prop.BM_Nf;
  if (Nf<3) Nf=3;
  df   = (fmax - fmin)/Nf;
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
    mds.get_all_neighborcollections();
    td = mds.calc_potential_energy() / mds.natoms();
    yp[j] = td;
    dyp[j] = ((td<0)? -td: td) * small;

    cout << "f g V Epa "
	 << f << " "
	 << g << " "
	 << xp[j] << " "
	 << yp[j] << endl;

    // Reset:
    mds.pos    = pos_bak;
    mds.boxdir = boxdir_bak;
    mds.boxlen = boxlen_bak;
  }
  exit(1);



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

  cond_conv.functolabs = param.p_potinfo->specs_prop.functolabs;
  cond_conv.functolrel = param.p_potinfo->specs_prop.functolrel;
  cond_conv.gradtolabs = param.p_potinfo->specs_prop.gradtolabs;
  cond_conv.steptolabs = param.p_potinfo->specs_prop.steptolabs;
  cond_conv.steptolrel = param.p_potinfo->specs_prop.steptolrel;
  cond_conv.nitermin   = param.p_potinfo->specs_prop.nitermin;
  cond_conv.nitermax   = param.p_potinfo->specs_prop.nitermax;
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


  int seed = param.p_potinfo->specs_prop.seed;

  if (param.p_potinfo->specs_prop.fitmet=="CG"){
    // Conjugate Gradients
    ConjGrad< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="PM"){
    // Powell's method
    Powell< ChiSqFunc<Param, double, double> > fm(cs);
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
    cs.barrier_scale() = 0.0;
    DiffEvol< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="PS"){
    // Particle Swarm
    cs.barrier_scale() = 0.0;
    PartSwarm< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="BC"){
    // Bee colony
    cs.barrier_scale() = 0.0;
    BeeColony< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="GS"){
    // Gravitational Search
    cs.barrier_scale() = 0.0;
    GravSearch< ChiSqFunc<Param, double, double> > fm(cs);
    Xopt = fm.minimize(cs.Param().X(),
		       cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
		       seed,
		       cond_conv, cond_debug, cond_print);
    fit_OK = fm.status.fit_OK;
  }
  else if (param.p_potinfo->specs_prop.fitmet=="SA"){
    // Simulated Annealing
    cs.barrier_scale() = 0.0;
    SimAnn< ChiSqFunc<Param, double, double> > fm(cs);

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
    fit_OK = fm.status.fit_OK;
  }
  else
    aborterror("Error: Unknown fitting method " + param.p_potinfo->specs_prop.fitmet + ". " +
	       "Exiting.");



  cs(Xopt);



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
    guess_B = B0 = d2E/(dV*dV) * V0;
    //B0  = d2E/(4.0*dV*dV) * V0 * eVA3_to_GPa;
    d3E = Ep[nmid+2] - 2*Ep[nmid+1] + 2*Ep[nmid-1] - Ep[nmid-2];
    Bp0 = -1.0 - V0*V0/B0 * d3E/(2.0*dV*dV*dV);
    guess_Bp = Bp0;

    B0  = guess_B * eVA3_to_GPa;
    Bp0 = guess_Bp;
    cout << "Warning: Failed to obtain bulk modulus and/or its pressure derivative for compound " + name + ". " +
      "Using numerical estimates " << B0 << " GPa and " << Bp0 << endl;
    // **************************************************************
    // **************************************************************
  }



  if (prop_use.Ecoh)  prop_pred.Ecoh  = E0;
  if (prop_use.Vatom) prop_pred.Vatom = V0;
  if (prop_use.B)     prop_pred.B     = B0;
  if (prop_use.Bp)    prop_pred.Bp    = Bp0;
  
}






void CompoundStructureFit::get_Cij(MDSystem             & mds,
				   ParamPot             & param,
				   Vector< Vector<double> > & pos_bak,
				   Matrix<double>           & boxdir_bak,
				   Vector<double>           & boxlen_bak,
				   double & E0,
				   double & V0
				   ){
		
  bool quit, fit_OK;
  int j,k,p,Nf,NC,isym;
  Vector<double> xp, yp, dyp, Clincomb, Xopt;
  Matrix<double> alpha(3,3,0), C(6,6,0);
  double fmin, fmax, df, f, td, C11, C12, C13, C22, C23, C33, C44, C55, C66;
  double small = sqrt( std::numeric_limits<double>::epsilon() ), min_f=0.0, min_E=E0, guess_C;


  fmin = param.p_potinfo->specs_prop.C_fmin;
  fmax = param.p_potinfo->specs_prop.C_fmax;
  Nf   = param.p_potinfo->specs_prop.C_Nf;
  if (Nf<3) Nf=3;
  df   = (fmax - fmin)/Nf;
  xp.resize(Nf); yp.resize(Nf); dyp.resize(Nf);


  /* ###############################################################################
     Find the elastic constants.
       
     For a cubic lattice we have the elastic constants
     C11, C12, C44
       
     For a hexagonal lattice we have the elastic constants
     C11, C12, C44, C13, C33
     ############################################################################### */

  if      (csystem=="cubic") NC=3;
  else if (csystem=="hexagonal") NC=5;
  else if (csystem=="orthorombic") NC=9;
  else NC=0;

  Clincomb.resize(NC);
  for (isym=0; isym<NC; ++isym) Clincomb[isym]=0;



  quit = false;

  // Loop over particular symmetry changes
  for (isym=0; isym<NC; ++isym){


    // Loop over differential changes
    for (j=0; j<Nf; ++j){
      f = fmin + j*df;

      /* Establish the transformation. */
      for (k=0; k<3; ++k) for (p=0; p<3; ++p) alpha.elem(k,p)=0;

      if (csystem=="cubic"){
	if (isym==0){
	  /* Result: E - E0 = V0 * (C11 + C12) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==1){
	  /* Result: E - E0 = V0 * (C11 - C12) * f^2*/
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0 - f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==2){
	  /* Result: E - E0 = 2 * V0 * C44 * f^2*/
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	  alpha.elem(2,1) = f;
	}
	else quit=true;
      }
      else if (csystem=="hexagonal"){
	if (isym==0){
	  /* Result: E - E0 = V0 * (C11 + C12) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==1){
	  /* Result: E - E0 = V0 * (C11 - C12) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0 - f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==2){
	  /* Result: E - E0 = V0 * C33/2 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==3){
	  /* Result: E - E0 = V0 * 0.5 * (C11 + C33 + 2*C13) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==4){
	  /* Result: E - E0 = 2 * V0 * C44 * f^2*/
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	  alpha.elem(2,1) = f;
	}
	else quit=true;
      }
      else if (csystem=="orthorombic"){
	if (isym==0){
	  /* Result: E - E0 = 0.5 * V0 * C11 * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==1){
	  /* Result: E - E0 = 0.5 * V0 * C33 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==2){
	  /* Result: E - E0 = 0.5 * V0 * C33 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==3){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + C22 + 2 * C12) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==4){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + C33 + 2 * C13) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==5){
	  /* Result: E - E0 = 0.5 * V0 * (C22 + C33 + 2 * C23) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==6){
	  /* Result: E - E0 = 2.0 * V0 * C44 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	  alpha.elem(2,1) = f;
	}
	else if (isym==7){
	  /* Result: E - E0 = 2.0 * V0 * C55 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,2) = f;
	  alpha.elem(2,0) = f;
	}
	else if (isym==8){
	  /* Result: E - E0 = 2.0 * V0 * C66 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,1) = f;
	  alpha.elem(1,0) = f;
	}
	else quit=true;
      }
      else quit=true;

      if (quit) break;

    
      xp[j] = f;
      //xp[j] = V0 * g*g*g;
      //for (k=0; k<3; ++k) for (p=0; p<3; ++p) alpha.elem(k,p)=0;
      //alpha.elem(0,0) = alpha.elem(1,1) = alpha.elem(2,2) = g; // Note!!!

      mds.transform_cell(alpha);
      mds.get_all_neighborcollections();
      td = mds.calc_potential_energy() / mds.natoms();
      yp[j] = td;
      dyp[j] = ((td<0)? -td: td) * small;


      // Reset:
      mds.pos    = pos_bak;
      mds.boxdir = boxdir_bak;
      mds.boxlen = boxlen_bak;
    }

    if (quit) break;



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

    cond_conv.functolabs = param.p_potinfo->specs_prop.functolabs;
    cond_conv.functolrel = param.p_potinfo->specs_prop.functolrel;
    cond_conv.gradtolabs = param.p_potinfo->specs_prop.gradtolabs;
    cond_conv.steptolabs = param.p_potinfo->specs_prop.steptolabs;
    cond_conv.steptolrel = param.p_potinfo->specs_prop.steptolrel;
    cond_conv.nitermin   = param.p_potinfo->specs_prop.nitermin;
    cond_conv.nitermax   = param.p_potinfo->specs_prop.nitermax;
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


    int seed = param.p_potinfo->specs_prop.seed;

    if (param.p_potinfo->specs_prop.fitmet=="CG"){
      // Conjugate Gradients
      ConjGrad< ChiSqFunc<Param, double, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 cond_conv, cond_debug, cond_print);
      fit_OK = fm.status.fit_OK;
    }
    else if (param.p_potinfo->specs_prop.fitmet=="PM"){
      // Powell's method
      Powell< ChiSqFunc<Param, double, double> > fm(cs);
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
      cs.barrier_scale() = 0.0;
      DiffEvol< ChiSqFunc<Param, double, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 seed,
			 cond_conv, cond_debug, cond_print);
      fit_OK = fm.status.fit_OK;
    }
    else if (param.p_potinfo->specs_prop.fitmet=="PS"){
      // Particle Swarm
      cs.barrier_scale() = 0.0;
      PartSwarm< ChiSqFunc<Param, double, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 seed,
			 cond_conv, cond_debug, cond_print);
      fit_OK = fm.status.fit_OK;
    }
    else if (param.p_potinfo->specs_prop.fitmet=="BC"){
      // Bee colony
      cs.barrier_scale() = 0.0;
      BeeColony< ChiSqFunc<Param, double, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 seed,
			 cond_conv, cond_debug, cond_print);
      fit_OK = fm.status.fit_OK;
    }
    else if (param.p_potinfo->specs_prop.fitmet=="GS"){
      // Gravitational Search
      cs.barrier_scale() = 0.0;
      GravSearch< ChiSqFunc<Param, double, double> > fm(cs);
      Xopt = fm.minimize(cs.Param().X(),
			 cs.Param().Xmin(), cs.Param().Xmax(), cs.Param().Xtype(),
			 seed,
			 cond_conv, cond_debug, cond_print);
      fit_OK = fm.status.fit_OK;
    }
    else if (param.p_potinfo->specs_prop.fitmet=="SA"){
      // Simulated Annealing
      cs.barrier_scale() = 0.0;
      SimAnn< ChiSqFunc<Param, double, double> > fm(cs);

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
      fit_OK = fm.status.fit_OK;
    }
    else
      aborterror("Error: Unknown fitting method " + param.p_potinfo->specs_prop.fitmet + ". " +
		 "Exiting.");










    // If fitting went alright, then we are almost done. Else we need to get
    // a numerical estimate of the required derivatives.

    if (fit_OK){
      E0  = Xopt[0];
      //f0  = Xopt[1];
      Clincomb[isym] = Xopt[2] / V0 * eVA3_to_GPa; // unit now: GPa
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
      Vector<double> vx(3,0), vy(3,0);
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
      Clincomb[isym] = guess_C / V0 * eVA3_to_GPa; // unit now: GPa
      cout << "Warning: Failed to obtain combination of linear constants for compound " + name + ". " +
	"Using estimate " << guess_C << " GPa." << endl;
    }

  } // End of loop over isym



  
  if (csystem=="cubic"){
    C44 = Clincomb[2]/2.0;
    C11 = (Clincomb[0] + Clincomb[1])/2.0;
    C12 = Clincomb[0] - C11;
    C.elem(0,0) = C11;
    C.elem(0,1) = C12;
    C.elem(3,3) = C44;
  }
  else if (csystem=="hexagonal"){
    C33 = 2.0 * Clincomb[2];
    C11 = (Clincomb[0] + Clincomb[1])/2.0;
    C12 = Clincomb[0] - C11;
    C13 = (2.0*Clincomb[3] - C11 - C33)/2.0;
    C44 = Clincomb[4]/2.0;
    C.elem(0,0) = C11;
    C.elem(0,1) = C12;
    C.elem(3,3) = C44;
    C.elem(2,2) = C33;
    C.elem(0,2) = C13;
  }
  else if (csystem=="orthorombic"){
    C11 = 2.0 * Clincomb[0];
    C22 = 2.0 * Clincomb[1];
    C33 = 2.0 * Clincomb[2];
    C12 = 0.5 * ( 2.0 * Clincomb[3] - C11 - C22 );
    C13 = 0.5 * ( 2.0 * Clincomb[4] - C11 - C33 );
    C23 = 0.5 * ( 2.0 * Clincomb[5] - C22 - C33 );
    C44 = Clincomb[6]/2.0;
    C55 = Clincomb[7]/2.0;
    C66 = Clincomb[8]/2.0;
    C.elem(0,0) = C11;
    C.elem(1,1) = C22;
    C.elem(2,2) = C33;
    C.elem(0,1) = C12;
    C.elem(0,2) = C13;
    C.elem(1,2) = C23;
    C.elem(3,3) = C44;
    C.elem(4,4) = C55;
    C.elem(5,5) = C66;
  }


  for (k=0; k<6; ++k){
    for (p=0; p<6; ++p){
      if (prop_use.C.elem(k,p))
	prop_pred.C.elem(k,p) = C.elem(k,p);
    }
  }


}



