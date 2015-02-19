



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
#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mdsettings.hpp"
#include "mtwister.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "potinfo.hpp"
#include "specs-fit-prop-pot.hpp"
#include "errors.hpp"

using namespace std;
using namespace utils;
using namespace constants;
using boost::format;
using std::ofstream;





// ##############################################################################
// ##############################################################################
//
// Main MD routine
//
// ##############################################################################
// ##############################################################################


void MDSystem::relax(void){

  int i,j;
  Vector<double> tv(3,0.0), tv1(3,0.0), tv2(3,0.0);
  Vector<double> boxlen_orig(3,0);
  double tmp1, tmp2, tmp3, tmp4, td;
  double eps = std::numeric_limits<double>::epsilon();
  Vector<double> drs(3), drc(3);
  Vector<double> u(3), v(3), w(3);
  double Ek_tot, lambda, Tnew;
  Matrix<double> W(3,3,0);
  double mu[3];
  double drsq, drsq_max, drsq_max2;
  double dt_v_max=-1, dt_a_max=-1, dt_F_max=-1;
  Vector<double> dt_candidate(4, specs.max_dt);
  int nat;
  bad_mds err_bad_mds;


  bool periodic = true;
  if (!pbc[0] || !pbc[1] || !pbc[2]) periodic = false;


  const double liml=1.0e-6, limu=1.0e+6;
  string fmt, fmtf = "%12.6f", fmte = "%12.6e";







  get_pot_force  = true;
  get_pot_energy = true;;

  boxlen_orig = boxlen;


  {
    ofstream fdump;
    string dumpfn = "mds-frame0-" + name + ".xyz";
    fdump.open(dumpfn.c_str());
    dumpframe(fdump);
    fdump.close();
  }

  handle_pbc_of_positions();

  {
    ofstream fdump;
    string dumpfn = "mds-frame1-" + name + ".xyz";
    fdump.open(dumpfn.c_str());
    dumpframe(fdump);
    fdump.close();
  }


  calc_volume();
  cout << "Atomic volume " << vol_atom << endl;

  int type1;
  double mass1;

  
  type1 = elem.atomtype(matter[0]);
  mass1 = elem.mass(matter[0]); //type1


  rand_mtwister mtwister(specs.seed);
  Vector<double> vrms(elem.nelem(), 0.0);
  td = sqrt( 8.817343e-5 * 2 * 0.5*specs.Tstart * 1.602 / 1.661 * 0.01 );
  if (sys_single_elem)
    vrms[0] = td / sqrt( mass1 );
  else {
    for (j=0; j<elem.nelem(); j++)
      vrms[j] = td / sqrt( elem.mass(elem.idx2name(j)) );
  }



  Vector<int> nspecies( elem.nelem(), 0);
  Vector<string> elems_present;


  if (pos.size()==0)
    aborterror("Error: Trying to initilize a MDSystem containing 0 atoms. Exiting.");
  
  nat = natoms();

  vel.resize(nat);
  acc.resize(nat);
  frc.resize(nat);
  virials.resize(nat);
  Ep.resize(nat);
  Ek.resize(nat);
  dpos.resize(nat);
  pos_int_ini.resize(nat);
  pos_int_fin.resize(nat);

  if (specs_common.debug_forces){
    frc_num.resize(nat);
  }




  // Memory allocation:
  for (i=0; i<nat; ++i){
    vel[i] = Vector<double>(3, 0.0);
    acc[i] = Vector<double>(3, 0.0);
    frc[i] = Vector<double>(3, 0.0);
    virials[i] = Matrix<double>(3, 3, 0.0);
    Ek[i] = 0.0;
    Ep[i] = 0.0;
    dpos[i] = Vector<double>(3, 0.0);
    pos_int_ini[i] = Vector<double>(3, 0.0);
    pos_int_fin[i] = Vector<double>(3, 0.0);
    // --------------------------------------------------
    get_coords_cart2skew(pos[i], pos_int_ini[i], -1);
    pos_int_ini[i][0] /= boxlen[0];
    pos_int_ini[i][1] /= boxlen[1];
    pos_int_ini[i][2] /= boxlen[2];
    // --------------------------------------------------
    if (specs_common.debug_forces)
      frc_num[i] = Vector<double>(3,0);

    if (sys_single_elem) type[i] = type1;
    else                 type[i] = elem.atomtype(matter[i]);//name2idx(matter[i]);
    
    tmp1 = mtwister.gauss(); // has internal state which is updated after call
    tmp2 = mtwister.gauss();
    tmp3 = mtwister.gauss();
    int n2i = elem.name2idx(matter[i]);
    vel[i][0] = vrms[ n2i ] * tmp1; // type[i]
    vel[i][1] = vrms[ n2i ] * tmp2;
    vel[i][2] = vrms[ n2i ] * tmp3;

    nspecies[ n2i ]++;
  }
  //cout << "MD system vectors resized to correct sizes." << endl;



  drsq_max = drsq_max2 = 0.0;

  // Original positions should be recorded in a way that is not dependent on the
  // actual box lengths.
  // Internal position in abc system is (ra, rb, rc). Cartesian position is
  //
  // r_xyz = ra * a + rb * b + rc * c


 
  //cout << "Getting internal positions of atoms at start ..." << endl;
  //  Matrix<double> int_pos_matrix(3,3,0), dummy_matrix(3,3,0);



  //cout << "Getting initial velocities ..." << endl;



  for (i=0; i<nspecies.size(); ++i){
    if (nspecies[i]>0) elems_present.push_back( elem.idx2name(i) );
  }
  rcut = rcut_max = p_potinfo->get_rcut_max( elems_present );
  //cout << "rcut_max and skint are " << rcut_max << " " << skint << endl;
 

  double quench_tstart_real = 0.0;
  double quench_Tstart = T;
  bool   quench_active = false;






  /* -------------------------------------------------------------
     Get neighbors.
     ------------------------------------------------------------- */
  //cout << "Getting neighbors of all atoms ..." << endl;
  get_all_neighborcollections( specs_common.report_step );


  T = 0.0;
  P = Px = Py = Pz = 0.0;
  Ep_tot = 0.0;
  Ek_tot = 0.0;
  dt = specs.dt;
  time = 0.0;

  //cout << "made it here 01*" << endl;
  // calc_volume();
  if (!periodic)
    calc_closepacked_volume();
  else {
    calc_volume();
  }

  //cout << "made it here 02*" << endl;
  // calc_P();

  /*
  cout << "natoms: " << nat << endl;
  cout << "boxlen[0]: " << boxlen[0] << endl;
  cout << "boxlen[1]: " << boxlen[1] << endl;
  cout << "boxlen[2]: " << boxlen[2] << endl;
  cout << "type1: " << type1 << endl;
  cout << "mass1: " << mass1 << endl;
  */




  //  printf("nat_fixed = %ld       volume = %20.10f\n", nat_fixed, V);
  //  printf("P_vir(GPa): %20.10f  %20.10f  %20.10f\n", W_x/ V * 160.21773, W_y/ V * 160.21773, W_z/ V * 160.21773);

  get_virials(nat, W);

  /* Calculate Cartesian pressure components: */
  Px = stresstensor_xyz.elem(0,0) = W.elem(0,0);  /* Unit now: GPa. */
  Py = stresstensor_xyz.elem(1,1) = W.elem(1,1);  /* Unit now: GPa. */
  Pz = stresstensor_xyz.elem(2,2) = W.elem(2,2);  /* Unit now: GPa. */
  P = 1.0/3.0 * (Px + Py + Pz);

  //cout << "Stress tensor (Cartesian system): " << stresstensor_xyz << endl;
  //printf("Px Py Pz:  %15.10e  %15.10e  %15.10e\n", Px, Py, Pz);
  //  fflush(stdout);
  
  /* Get pressure components in skewed system: */
  {
    double td0,td1,td2;
    double te0,te1,te2;
    td0 = stresstensor_xyz.elem(0,0);
    td1 = stresstensor_xyz.elem(1,1);
    td2 = stresstensor_xyz.elem(2,2);

    // Default to Cartesian box:
    stresstensor_abc.elem(0,0) = td0;
    stresstensor_abc.elem(1,1) = td1;
    stresstensor_abc.elem(2,2) = td2;

    te0 = td0;
    te1 = td1;
    te2 = td2;
    if (! isCart){
      te0 = 0.0
	+ Bravaismatrix_inv.elem(0,0) * td0
	+ Bravaismatrix_inv.elem(0,1) * td1
	+ Bravaismatrix_inv.elem(0,2) * td2;
      te1 = 0.0
	+ Bravaismatrix_inv.elem(1,0) * td0
	+ Bravaismatrix_inv.elem(1,1) * td1
	+ Bravaismatrix_inv.elem(1,2) * td2;
      te2 = 0.0
	+ Bravaismatrix_inv.elem(2,0) * td0
	+ Bravaismatrix_inv.elem(2,1) * td1
	+ Bravaismatrix_inv.elem(2,2) * td2;
    }
    stresstensor_abc.elem(0,0) = te0;
    stresstensor_abc.elem(1,1) = te1;
    stresstensor_abc.elem(2,2) = te2;
  }


  //cout << "Stress tensor (skew system): " << stresstensor_abc << endl;
  //    printf("  Px Py Pz = %20.10f  %20.10f  %20.10f\n", Px, Py, Pz);
  //printf("  Pa Pb Pc = %20.10f  %20.10f  %20.10f\n", *Pa, *Pb, *Pc);




  //cout << "made it here 03*" << endl;
  //calc_T();
  Ek_tot = 0.0;
  if (sys_single_elem){
#pragma omp parallel for reduction(+:Ek_tot) schedule(static)
    for (i=0; i<nat; ++i){
      Ek_tot += mass1 * ( vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] );
    }
  }
  else {
#pragma omp parallel for reduction(+:Ek_tot) schedule(static)
    for (i=0; i<nat; ++i){
      int n2i = elem.name2idx( matter[i] );
      Ek_tot += elem.mass( n2i ) * ( vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] );
      //Ek_tot += elem.mass(matter[i]) * td;
    }
  }
  Ek_tot *= 0.5;
  T = 2.0 * Ek_tot / (3.0 * nat * 8.817343e-5) * 1.660538782/1.60217653 * 100;
  // amu * Ang^2/fs^2 = amu * 1e-20/1e-30 m^2/s^2 = amu * 1e10 m^2/s^2
  // = 1.660538782e-27 * 1e10 kg m^2/s^2
  // = 1.660538782e-17 J
  // = 1.660538782e-17 1/1.60217653e-19 eV
  // = 1.660538782/1.60217653e * 100 eV
  
  //cout << "made it here 04*" << endl;






  /*
    cout << "boxlen([0]: " << boxlen[0] << endl;
    cout << "boxlen([1]: " << boxlen[1] << endl;
    cout << "boxlen([2]: " << boxlen[2] << endl;
    cout << "boxdir([0]: " << boxdir.col(0) << endl;
    cout << "boxdir([1]: " << boxdir.col(1) << endl;
    cout << "boxdir([2]: " << boxdir.col(2) << endl;
  */




  /* #############################################################
     #############################################################
     #############################################################

     LOOP OVER TIME

     #############################################################
     #############################################################
     ############################################################# */

  int istep = 0;
  //cout << "Looping over time steps ..." << endl;

  drsq_max = drsq_max2 = 0.0;

  while (true){
    //  while (fabs(P) > 1e-3 || Fmax > 1e-3 ){

    /* ----------------------------------------------------------------------
       Velocity-Verlet predictor,
       time step check,
       neighbor list update check
       ---------------------------------------------------------------------- */
    //cout << "made it here 01" << endl;
    //predict();


#pragma omp parallel for schedule(static)
    for (i=0; i<nat; ++i){

      double td1 = dt * vel[i][0] + 0.5 * dt*dt * acc[i][0];
      double td2 = dt * vel[i][1] + 0.5 * dt*dt * acc[i][1];
      double td3 = dt * vel[i][2] + 0.5 * dt*dt * acc[i][2];

      if (specs.fixed_geometry)
	td1 = td2 = td3 = 0.0;

      // cache trashing???
      pos[i][0] += td1;
      pos[i][1] += td2;
      pos[i][2] += td3;

      dpos[i][0] += td1;
      dpos[i][1] += td2;
      dpos[i][2] += td3;

      if (specs.heating_allowed){
	vel[i][0] += 0.5 * dt * acc[i][0];
	vel[i][1] += 0.5 * dt * acc[i][1];
	vel[i][2] += 0.5 * dt * acc[i][2];
      }
      if (specs.quench_always){
	vel[i][0] = 0.0;
	vel[i][1] = 0.0;
	vel[i][2] = 0.0;
      }
	
      frc[i][0] = 0.0;
      frc[i][1] = 0.0;
      frc[i][2] = 0.0;

      acc[i][0] = 0.0;
      acc[i][1] = 0.0;
      acc[i][2] = 0.0;



    }


    for (i=0; i<nat; ++i){
      // time step
      td = vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2];
      if (td>dt_v_max*dt_v_max) dt_v_max=sqrt(td);
      dt_a_max = 0.0; //td = acc[i][0]*acc[i][0] + acc[i][1]*acc[i][1] + acc[i][2]*acc[i][2]; if (td>a_max) a_max=td;
      dt_F_max = 0.0; //td = frc[i][0]*frc[i][0] + frc[i][1]*frc[i][1] + frc[i][2]*frc[i][2]; if (td>F_max) F_max=td;

      // neighborlist stuff
      drsq = dpos[i][0] * dpos[i][0] + dpos[i][1] * dpos[i][1] + dpos[i][2] * dpos[i][2];
      if (drsq > drsq_max){
	drsq_max2 = drsq_max;
	drsq_max  = drsq;
      }
      else {
	if (drsq > drsq_max2){
	  drsq_max2 = drsq;
	}
      }
    }



    /*
      cout << "After predict():" << endl;
      cout << "boxlen([0]: " << boxlen[0] << endl;
      cout << "boxlen([1]: " << boxlen[1] << endl;
      cout << "boxlen([2]: " << boxlen[2] << endl;
      cout << "boxdir([0]: " << boxdir.col(0) << endl;
      cout << "boxdir([1]: " << boxdir.col(1) << endl;
      cout << "boxdir([2]: " << boxdir.col(2) << endl;
    */

    //cout << "made it here 02" << endl;



    handle_pbc_of_positions();




    /*
      cout << "After handle_pbc...():" << endl;
      cout << "boxlen([0]: " << boxlen[0] << endl;
      cout << "boxlen([1]: " << boxlen[1] << endl;
      cout << "boxlen([2]: " << boxlen[2] << endl;
      cout << "boxdir([0]: " << boxdir.col(0) << endl;
      cout << "boxdir([1]: " << boxdir.col(1) << endl;
      cout << "boxdir([2]: " << boxdir.col(2) << endl;
    */

    //cout << "made it here 03" << endl;
    //check_timestep();
    
    if (dt_v_max > 0.0)
      dt_candidate[0] = specs.max_dr / dt_v_max;
    if (dt_F_max > 0.0 && dt_v_max > 0.0)
      dt_candidate[1] = specs.max_dE / (dt_F_max * dt_v_max);
    dt_candidate[2] = 1.1 * dt;
    dt_candidate[3] = specs.max_dt;

    // Find smallest choice for dt:
    for (int k=0; k<4; ++k){
      if (k==0 || (k>0 && dt_candidate[k]<dt)) dt=dt_candidate[k];
    }


    //cout << "made it here 04" << endl;

    
    /* ----------------------------------------------------------------------
       Check for maximal displacement of atoms, and determine if neighbor list
       needs to be rebuilt.
       ---------------------------------------------------------------------- */
    //check_for_neighbors_update();

    if (drsq_max + drsq_max2 > skint){
      // printf("Updating neighbor list ...\n");
      //fflush(stdout);
      /* Update neighbor list. */
      get_all_neighborcollections( specs_common.report_step );

      for (i=0; i<nat; ++i){
	dpos[i][0] = 0.0;
	dpos[i][1] = 0.0;
	dpos[i][2] = 0.0;
      }
      drsq_max = drsq_max2 = 0.0;
    }


    /*
      cout << "After check_for_neighbors...():" << endl;
      cout << "boxlen([0]: " << boxlen[0] << endl;
      cout << "boxlen([1]: " << boxlen[1] << endl;
      cout << "boxlen([2]: " << boxlen[2] << endl;
      cout << "boxdir([0]: " << boxdir.col(0) << endl;
      cout << "boxdir([1]: " << boxdir.col(1) << endl;
      cout << "boxdir([2]: " << boxdir.col(2) << endl;
    */
    //cout << "made it here 05" << endl;


    /* ######################################################################
       Get numerical forces
       ###################################################################### */
    if (specs_common.debug_forces){
      // Analytical forces will be calculated after this step.
      Vector<double> pos_bak(3,0);
      double t, d, Ep1, Ep2;

      cout << "Debugging forces" << endl;

      t = pow(eps, 1.0/3.0);

      // Displace one atom at a time, and get total potential energy.
      // Numerical derivative of potential energy wrt displacement is the numerical force.
      for (i=0; i<nat; ++i){
	cout << " " << i; cout.flush();
	pos_bak[0] = pos[i][0];
	pos_bak[1] = pos[i][1];
	pos_bak[2] = pos[i][2];

	//cout << "debug_forces: Made it here 01, atom " << i << endl;
	for (int k=0; k<3; ++k){
	  //cout << " " << k;
	  d = abs(pos[i][k]);
	  //d = (d < eps) ? t : t*d;
	  d = (d < eps) ? t*eps : t*d;

	  //cout << " backstep ";
	  pos[i][k] = pos_bak[k] - d;
	  Ep1 = calc_potential_energy();
	  //cout << " done ";
	  //cout << " forwardstep ";
	  pos[i][k] = pos_bak[k] + d;
	  Ep2 = calc_potential_energy();
	  //cout << " done ";
	  frc_num[i][k] = - (Ep2 - Ep1)/(2.0*d);
	  //cout << " reset ";
	  pos[i][k] = pos_bak[k]; // reset
	}
      }
      cout << endl;
    }

    //cout << "made it here 05.5" << endl;



    /* ######################################################################
       Get numerical pressure
       ###################################################################### */
    if (specs_common.debug_pressure){
      
      Vector< Vector<double> > pos_bak = pos;
      
      double bx=boxlen[0], by=boxlen[1], bz=boxlen[2];
      double t = pow(eps, 1.0/3.0);
      double d,dx,dy,dz,rx,ry,rz;
      double V0,Ep1,Ep2;

      V0 = bx*by*bz;
      double s = t / (3.0);
      dx = s * bx;
      dy = s * by;
      dz = s * bz;

      // Point 1:
      boxlen[0] += dx;
      boxlen[1] += dy;
      boxlen[2] += dz;
      // Fractional change:
      rx = boxlen[0] / bx;
      ry = boxlen[1] / by;
      rz = boxlen[2] / bz;
      for (i=0; i<nat; ++i){
	pos[i][0] *= rz;
	pos[i][1] *= ry;
	pos[i][2] *= rz;
      }
      Ep1 = calc_potential_energy();

      // Reset:
      boxlen[0] = bx;
      boxlen[1] = by;
      boxlen[2] = bz;
      for (i=0; i<nat; ++i){
	pos[i][0] = pos_bak[i][0];
	pos[i][1] = pos_bak[i][1];
	pos[i][2] = pos_bak[i][2];
      }

      // Point 2:
      boxlen[0] -= dx;
      boxlen[1] -= dy;
      boxlen[2] -= dz;
      // Fractional change:
      rx = boxlen[0] / bx;
      ry = boxlen[1] / by;
      rz = boxlen[2] / bz;
      for (i=0; i<nat; ++i){
	pos[i][0] *= rx;
	pos[i][1] *= ry;
	pos[i][2] *= rz;
      }
      Ep2 = calc_potential_energy();

      // Reset:
      boxlen[0] = bx;
      boxlen[1] = by;
      boxlen[2] = bz;
      for (i=0; i<nat; ++i){
	pos[i][0] = pos_bak[i][0];
	pos[i][1] = pos_bak[i][1];
	pos[i][2] = pos_bak[i][2];
      }
      
      cout << "Numerical pressure is " << - (Ep1 - Ep2)/(2.0*V0*3*s) * eVA3_to_GPa << endl;
    }


    //cout << "made it here 05.5" << endl;

    /* ----------------------------------------------------------------------
       Get forces and energies:
       ---------------------------------------------------------------------- */
    //calc_forces_and_energies();
    get_pot_force  = true;
    get_pot_energy = true;
    get_forces_and_energies_common();

    /*
      cout << "After calc_forces...():" << endl;
      cout << "boxlen([0]: " << boxlen[0] << endl;
      cout << "boxlen([1]: " << boxlen[1] << endl;
      cout << "boxlen([2]: " << boxlen[2] << endl;
      cout << "boxdir([0]: " << boxdir.col(0) << endl;
      cout << "boxdir([1]: " << boxdir.col(1) << endl;
      cout << "boxdir([2]: " << boxdir.col(2) << endl;
    */

    //cout << "made it here 06" << endl;
    // printf("After energy calculation: Box now: box1 box2 box3  %10.5e %10.5e %10.5e\n", box[1], box[2], box[3]); fflush(stdout);


    if (specs_common.debug_forces) break;


    /* ----------------------------------------------------------------------
       Get maximum force:
       ---------------------------------------------------------------------- */
    //if (debug_mds_prop==True) printf("Getting max. force ...\n");
    {
      int s1, s2, s3, s4;
      s1=s2=s3=s4=1;

      F_max = 0;
      for (i=0; i<nat; ++i){
	s1 = (frc[i][0] < 0) ? -1: 1;
	s2 = (frc[i][1] < 0) ? -1: 1;
	s3 = (frc[i][2] < 0) ? -1: 1;
	tmp1 = s1*frc[i][0];
	tmp2 = s2*frc[i][1];
	tmp3 = s3*frc[i][2];
	tmp4 = 0;

	if      (tmp1>=tmp2 && tmp1>=tmp3) tmp4 = tmp1;
	else if (tmp2>=tmp1 && tmp2>=tmp3) tmp4 = tmp2;
	else                               tmp4 = tmp3;

	if (i==0 || (tmp4 >= F_max)) F_max = tmp4;
      }
    }
  
      


    
    
    /* ****************************************************************
       Perform internal relaxation.
       **************************************************************** */
    //if (debug_mds_prop==True) printf("Velocity Verlet corrector phase ...\n");

    /* ----------------------------------------------------------------------
       Velocity-Verlet corrector, and conversion of force to acceleration.
       ---------------------------------------------------------------------- */
    //correct();

    Ek_tot = 0.0;
#pragma omp parallel for reduction(+:Ek_tot) schedule(static)
    for (i=0; i<nat; ++i){
      int n2i = elem.name2idx( matter[i] );
      double td = 1.0 / elem.mass( n2i ) * 1.60217653 / 1.660538782 * 0.01;
      //td = 1.0 / elem.mass(elem.idx2name(type[i])) * 1.60217653 / 1.660538782 * 0.01;

      /* Convert force to acceleration. */
      acc[i][0] = frc[i][0] * td;
      acc[i][1] = frc[i][1] * td;
      acc[i][2] = frc[i][2] * td;

      /* Get the corrected velocity. */
      if (specs.heating_allowed){
	vel[i][0] += 0.5 * dt * acc[i][0];
	vel[i][1] += 0.5 * dt * acc[i][1];
	vel[i][2] += 0.5 * dt * acc[i][2];
      }
      if (specs.quench_always){
	vel[i][0] = 0.0;
	vel[i][1] = 0.0;
	vel[i][2] = 0.0;
      }
      // temperature
      Ek_tot += elem.mass( n2i ) * ( vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] );
      //Ek_tot += elem.mass(matter[i]) * td;
    }
    Ek_tot *= 0.5;
    T = 2.0 * Ek_tot / (3.0 * nat * 8.817343e-5) * 1.660538782/1.60217653 * 100;


    for (i=0; i<nat; ++i){
      // time step check
      td = vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2]; if (td>dt_v_max*dt_v_max) dt_v_max=sqrt(td);
      td = acc[i][0]*acc[i][0] + acc[i][1]*acc[i][1] + acc[i][2]*acc[i][2]; if (td>dt_a_max*dt_a_max) dt_a_max=sqrt(td);
      td = frc[i][0]*frc[i][0] + frc[i][1]*frc[i][1] + frc[i][2]*frc[i][2]; if (td>dt_F_max*dt_F_max) dt_F_max=sqrt(td);
    }

    get_virials(nat, W);


    // amu * Ang^2/fs^2 = amu * 1e-20/1e-30 m^2/s^2 = amu * 1e10 m^2/s^2
    // = 1.660538782e-27 * 1e10 kg m^2/s^2
    // = 1.660538782e-17 J
    // = 1.660538782e-17 1/1.60217653e-19 eV
    // = 1.660538782/1.60217653e * 100 eV



    /*
      cout << "After correct():" << endl;
      cout << "boxlen([0]: " << boxlen[0] << endl;
      cout << "boxlen([1]: " << boxlen[1] << endl;
      cout << "boxlen([2]: " << boxlen[2] << endl;
      cout << "boxdir([0]: " << boxdir.col(0) << endl;
      cout << "boxdir([1]: " << boxdir.col(1) << endl;
      cout << "boxdir([2]: " << boxdir.col(2) << endl;
    */

    //cout << "made it here 07" << endl;
    //check_timestep();
    
    if (dt_v_max > 0.0)
      dt_candidate[0] = specs.max_dr / dt_v_max;
    if (dt_F_max > 0.0 && dt_v_max > 0.0)
      dt_candidate[1] = specs.max_dE / (dt_F_max * dt_v_max);
    dt_candidate[2] = 1.1 * dt;
    dt_candidate[3] = specs.max_dt;

    // Find smallest choice for dt:
    for (int k=0; k<4; ++k){
      if (k==0 || (k>0 && dt_candidate[k]<dt)) dt=dt_candidate[k];
    }


    //cout << "made it here 08" << endl;



    /* ----------------------------------------------------------------------
       Get some properties:
       ---------------------------------------------------------------------- */
    //if (debug_mds_prop==True) printf("Getting pressure and temperature ...\n");
    if (!periodic)
      calc_closepacked_volume();
    else {
      calc_volume();
    }
    /*
      cout << "After calc_volume():" << endl;
      cout << "boxlen([0]: " << boxlen[0] << endl;
      cout << "boxlen([1]: " << boxlen[1] << endl;
      cout << "boxlen([2]: " << boxlen[2] << endl;
      cout << "boxdir([0]: " << boxdir.col(0) << endl;
      cout << "boxdir([1]: " << boxdir.col(1) << endl;
      cout << "boxdir([2]: " << boxdir.col(2) << endl;
      cout << "V = " << V << endl;
      cout << "V prefactor = " << td << endl;
    */

    P = 0.0;

    // calc_P();


    //  printf("nat_fixed = %ld       volume = %20.10f\n", nat_fixed, V);
    //  printf("P_vir(GPa): %20.10f  %20.10f  %20.10f\n", W_x/ V * 160.21773, W_y/ V * 160.21773, W_z/ V * 160.21773);


    /* Calculate Cartesian pressure components: */
    Px = stresstensor_xyz.elem(0,0) = W.elem(0,0);  /* Unit now: GPa. */
    Py = stresstensor_xyz.elem(1,1) = W.elem(1,1);  /* Unit now: GPa. */
    Pz = stresstensor_xyz.elem(2,2) = W.elem(2,2);  /* Unit now: GPa. */
    P = 1.0/3.0 * (Px + Py + Pz);

    if (specs_common.debug_pressure){
      cout << "Virial matrix: " << W << endl;
    }

    /* Get pressure components in skewed system: */
    double td0,td1,td2;
    double te0,te1,te2;

    td0 = stresstensor_xyz.elem(0,0);
    td1 = stresstensor_xyz.elem(1,1);
    td2 = stresstensor_xyz.elem(2,2);

    // Default to Cartesian box:
    stresstensor_abc.elem(0,0) = td0;
    stresstensor_abc.elem(1,1) = td1;
    stresstensor_abc.elem(2,2) = td2;

    te0 = td0;
    te1 = td1;
    te2 = td2;
    if (! isCart){
      te0 = 0.0
	+ Bravaismatrix_inv.elem(0,0) * td0
	+ Bravaismatrix_inv.elem(0,1) * td1
	+ Bravaismatrix_inv.elem(0,2) * td2;
      te1 = 0.0
	+ Bravaismatrix_inv.elem(1,0) * td0
	+ Bravaismatrix_inv.elem(1,1) * td1
	+ Bravaismatrix_inv.elem(1,2) * td2;
      te2 = 0.0
	+ Bravaismatrix_inv.elem(2,0) * td0
	+ Bravaismatrix_inv.elem(2,1) * td1
	+ Bravaismatrix_inv.elem(2,2) * td2;
    }
    stresstensor_abc.elem(0,0) = te0;
    stresstensor_abc.elem(1,1) = te1;
    stresstensor_abc.elem(2,2) = te2;

    //cout << "Stress tensor (skew system): " << stresstensor_abc << endl;
    //    printf("  Px Py Pz = %20.10f  %20.10f  %20.10f\n", Px, Py, Pz);
    //printf("  Pa Pb Pc = %20.10f  %20.10f  %20.10f\n", *Pa, *Pb, *Pc);


    
    /*
      cout << "After calc_P():" << endl;
      cout << "boxlen([0]: " << boxlen[0] << endl;
      cout << "boxlen([1]: " << boxlen[1] << endl;
      cout << "boxlen([2]: " << boxlen[2] << endl;
      cout << "boxdir([0]: " << boxdir.col(0) << endl;
      cout << "boxdir([1]: " << boxdir.col(1) << endl;
      cout << "boxdir([2]: " << boxdir.col(2) << endl;
    */

    // calc_T();


    //cout << "made it here 09" << endl;


    // printf("After first energy-pressure calculation for this cell: Px Py Pz  %10.5e %10.5e %10.5e\n", Px, Py, Pz);
    // fflush(stdout);
    // if (quickmode == True) break;
    // Leave MD loop over time.
    //    printf("After pressure calculation: Box now: box1 box2 box3  %10.5e %10.5e %10.5e\n", box[1], box[2], box[3]); fflush(stdout);



    /* ****************************************************************
       Perform external relaxation.
       **************************************************************** */
	  
    /* ----------------------------------------------------------------------
       Pressure control.
       ---------------------------------------------------------------------- */
    if (periodic){
      if (specs.use_Pcontrol && !specs.fixed_geometry){
	// control_P();
	// Requires that P has been calculated earlier!

	/*
	  cout << "bpc_P0  bpc_scale  bpc_tau: "
	  << specs.bpc_P0 << " " << specs.bpc_scale << " " << specs.bpc_tau << endl;
	  cout << "stresstensor_xyz.elem(0,0) " << stresstensor_xyz.elem(0,0) << endl;
	  cout << "stresstensor_xyz.elem(1,1) " << stresstensor_xyz.elem(1,1) << endl;
	  cout << "stresstensor_xyz.elem(2,2) " << stresstensor_xyz.elem(2,2) << endl;
	*/

	double third = 1.0/3.0;

	/*
	td = specs.bpc_P0 - stresstensor_xyz.elem(0,0);
	cout << "specs.bpc_P0 - stresstensor_xyz.elem(0,0)  " << td << endl;
	td = specs.bpc_P0 - stresstensor_xyz.elem(1,1);
	cout << "specs.bpc_P0 - stresstensor_xyz.elem(1,1)  " << td << endl;
	td = specs.bpc_P0 - stresstensor_xyz.elem(2,2);
	cout << "specs.bpc_P0 - stresstensor_xyz.elem(2,2)  " << td << endl;
	*/

	td = third * (dt / (specs.bpc_scale * specs.bpc_tau ));
	//cout << "dt / (specs.bpc_scale * specs.bpc_tau ) " << td << endl;
	mu[0] = pow( 1.0 - td * (specs.bpc_P0 - stresstensor_xyz.elem(0,0)), 1.0/3.0);
	mu[1] = pow( 1.0 - td * (specs.bpc_P0 - stresstensor_xyz.elem(1,1)), 1.0/3.0);
	mu[2] = pow( 1.0 - td * (specs.bpc_P0 - stresstensor_xyz.elem(2,2)), 1.0/3.0);
	//cout << "mu(0,1,2): " << mu[0] << " " << mu[1] << " " << mu[2] << endl;

	//printf("Inside pressure control: Px Py Pz:  %15.10e  %15.10e  %15.10e\n", Px, Py, Pz);
	if (specs_common.debug_pressure){
	  printf("Inside pressure control: mx my mz:  %15.10e  %15.10e  %15.10e\n", mu[0], mu[1], mu[2]);
	  fflush(stdout);
	}

	//#pragma omp parallel for schedule(static)
	for (i=0; i<nat; ++i){
	  pos[i][0] *= mu[0];
	  pos[i][1] *= mu[1];
	  pos[i][2] *= mu[2];
	}

	/* -----------------------------------------------------------------------
	   Scaling of box directions:
	   ----------------------------------------------------------------------- */
	// The columns of boxdir are the vectors indicating the primitive directions (normalized).
	boxdir.elem(0, 0) *= mu[0]; // scale x axis 
	boxdir.elem(0, 1) *= mu[0]; // scale x axis 
	boxdir.elem(0, 2) *= mu[0]; // scale x axis 

	boxdir.elem(1, 0) *= mu[1]; // scale y axis 
	boxdir.elem(1, 1) *= mu[1]; // scale y axis 
	boxdir.elem(1, 2) *= mu[1]; // scale y axis 

	boxdir.elem(2, 0) *= mu[2]; // scale z axis 
	boxdir.elem(2, 1) *= mu[2]; // scale z axis 
	boxdir.elem(2, 2) *= mu[2]; // scale z axis 

	// Multiply in the box length:
	boxdir.elem(0,0) *= boxlen[0];
	boxdir.elem(1,0) *= boxlen[0];
	boxdir.elem(2,0) *= boxlen[0];

	boxdir.elem(0,1) *= boxlen[1];
	boxdir.elem(1,1) *= boxlen[1];
	boxdir.elem(2,1) *= boxlen[1];

	boxdir.elem(0,2) *= boxlen[2];
	boxdir.elem(1,2) *= boxlen[2];
	boxdir.elem(2,2) *= boxlen[2];

	// Normalize to get the new lengths:
	boxlen[0] = boxdir.col(0).normalize();
	boxlen[1] = boxdir.col(1).normalize();
	boxlen[2] = boxdir.col(2).normalize();


	
	/* Guard against box becoming too large when relaxing: */
	for (int k=0; k<3; ++k){
	  if (isnan(boxlen[k])){
	    cout << "ERROR: Pressure control: Box length in direction " << k
		 << " for compound " + name + " is NaN." << endl;
	    throw err_bad_mds;
	    //aborterror("ERROR: Box in direction " + tostring(i) + " is NaN. Exiting");
	  }
	}



	/*
	  cout << "Inside control_P()" << endl;
	  cout << "boxlen([0]: " << boxlen[0] << endl;
	  cout << "boxlen([1]: " << boxlen[1] << endl;
	  cout << "boxlen([2]: " << boxlen[2] << endl;
	  cout << "boxdir([0]: " << boxdir.col(0) << endl;
	  cout << "boxdir([1]: " << boxdir.col(1) << endl;
	  cout << "boxdir([2]: " << boxdir.col(2) << endl;
	*/    
	//cout << "Updating box geometry" << endl;
	update_box_geometry();

	/*
	  cout << "boxlen([0]: " << boxlen[0] << endl;
	  cout << "boxlen([1]: " << boxlen[1] << endl;
	  cout << "boxlen([2]: " << boxlen[2] << endl;
	  cout << "boxdir([0]: " << boxdir.col(0) << endl;
	  cout << "boxdir([1]: " << boxdir.col(1) << endl;
	  cout << "boxdir([2]: " << boxdir.col(2) << endl;

	  cout << "Doing calc_P()" << endl;
	*/


	if (!periodic)
	  calc_closepacked_volume();
	else {
	  calc_volume();
	}


	get_virials(nat, W, mu[0], mu[1], mu[2]);


	/* Calculate Cartesian pressure components: */
	Px = stresstensor_xyz.elem(0,0) = W.elem(0,0);  /* Unit now: GPa. */
	Py = stresstensor_xyz.elem(1,1) = W.elem(1,1);  /* Unit now: GPa. */
	Pz = stresstensor_xyz.elem(2,2) = W.elem(2,2);  /* Unit now: GPa. */
	P = 1.0/3.0 * (Px + Py + Pz);

	if (specs_common.debug_pressure){
	  cout << "After P control: Virial matrix: " << W << endl;
	}

	//cout << "Stress tensor (Cartesian system): " << stresstensor_xyz << endl;
	//printf("Px Py Pz:  %15.10e  %15.10e  %15.10e\n", Px, Py, Pz);
	//  fflush(stdout);
  
	/* Get pressure components in skewed system: */
	double td0,td1,td2;
	double te0,te1,te2;

	td0 = stresstensor_xyz.elem(0,0);
	td1 = stresstensor_xyz.elem(1,1);
	td2 = stresstensor_xyz.elem(2,2);

	// Default to Cartesian box:
	stresstensor_abc.elem(0,0) = td0;
	stresstensor_abc.elem(1,1) = td1;
	stresstensor_abc.elem(2,2) = td2;

	te0 = td0;
	te1 = td1;
	te2 = td2;
	if (! isCart){
	  te0 = 0.0
	    + Bravaismatrix_inv.elem(0,0) * td0
	    + Bravaismatrix_inv.elem(0,1) * td1
	    + Bravaismatrix_inv.elem(0,2) * td2;
	  te1 = 0.0
	    + Bravaismatrix_inv.elem(1,0) * td0
	    + Bravaismatrix_inv.elem(1,1) * td1
	    + Bravaismatrix_inv.elem(1,2) * td2;
	  te2 = 0.0
	    + Bravaismatrix_inv.elem(2,0) * td0
	    + Bravaismatrix_inv.elem(2,1) * td1
	    + Bravaismatrix_inv.elem(2,2) * td2;
	}
	stresstensor_abc.elem(0,0) = te0;
	stresstensor_abc.elem(1,1) = te1;
	stresstensor_abc.elem(2,2) = te2;
	
	/*
	  cout << "After control_P():" << endl;
	  cout << "boxlen([0]: " << boxlen[0] << endl;
	  cout << "boxlen([1]: " << boxlen[1] << endl;
	  cout << "boxlen([2]: " << boxlen[2] << endl;
	  cout << "boxdir([0]: " << boxdir.col(0) << endl;
	  cout << "boxdir([1]: " << boxdir.col(1) << endl;
	  cout << "boxdir([2]: " << boxdir.col(2) << endl;
	*/
	//cout << "made it here 10" << endl;



	/*
	cout << "After calc_volume():" << endl;
	printf("boxlen([0]: %15.10f\n", boxlen[0]);
	printf("boxlen([1]: %15.10f\n", boxlen[1]);
	printf("boxlen([2]: %15.10f\n", boxlen[2]);
	cout << "boxdir([0]: " << boxdir.col(0) << endl;
	cout << "boxdir([1]: " << boxdir.col(1) << endl;
	cout << "boxdir([2]: " << boxdir.col(2) << endl;
	*/


	/* Guard against box becoming too small when relaxing: */
	if (boxlen[0] <= 2.0*(rcut+skint)) N[0] = - floor(boxlen_orig[0]/boxlen[0] * N[0]);
	if (boxlen[1] <= 2.0*(rcut+skint)) N[1] = - floor(boxlen_orig[1]/boxlen[1] * N[1]);
	if (boxlen[2] <= 2.0*(rcut+skint)) N[2] = - floor(boxlen_orig[2]/boxlen[2] * N[2]);
	
	if (N[0]<0 || N[1]<0 || N[2]<0){
	  cout << "Warning: Some box length(s) is (are) less than 2 * (rcut + skint): N[0] N[1] N[2]: "
	       << N[0] << " "
	       << N[1] << " "
	       << N[2] << endl;
	  cout << "Warning: Returning to box construction and redoing relaxation with bigger box ..." << endl;
	  return;
	}
      }
    
    }





    /* ----------------------------------------------------------------------
       Get maximum pressure:
       ---------------------------------------------------------------------- */
    //if (debug_mds_prop==True) printf("Getting max. pressure ...\n");
    {
      int s1, s2, s3;
      s1=s2=s3=1;
      
      s1 = (stresstensor_xyz.elem(0,0) < 0) ? -1 : 1;
      s2 = (stresstensor_xyz.elem(1,1) < 0) ? -1 : 1;
      s3 = (stresstensor_xyz.elem(2,2) < 0) ? -1 : 1;

      tmp1 = s1*stresstensor_xyz.elem(0,0);
      tmp2 = s2*stresstensor_xyz.elem(1,1);
      tmp3 = s3*stresstensor_xyz.elem(2,2);

      P_max = tmp3;
      if      (tmp1 >= tmp2 && tmp1 >= tmp3) P_max = tmp1;
      else if (tmp2 >= tmp1 && tmp2 >= tmp3) P_max = tmp2;
    }
      
    //cout << "made it here 11" << endl;    
      
    
    /* ----------------------------------------------------------------------
       Temperature control:
       ---------------------------------------------------------------------- */
    if (specs.use_Tcontrol && specs.heating_allowed){
      // control_T();
      // Requires that T has been calculated earlier!
      lambda = 1.0;
      if (! fp_is_small(T))
	lambda = sqrt( 1.0 + (dt / specs.btc_tau) * (specs.btc_T0/T - 1.0) );

      /* Change velocities: */
#pragma omp parallel for schedule(static)
      for (i=0; i<nat; ++i){
	vel[i][0] *= lambda;
	vel[i][1] *= lambda;
	vel[i][2] *= lambda;
      }
      /* Change temperature: */
      Tnew = T * lambda * lambda;
      T = Tnew;
    }


    /* ----------------------------------------------------------------------
       Quench stuff:
       ---------------------------------------------------------------------- */
    if (specs.use_quench && time >= specs.quench_tstart && !quench_active){
      quench_active = true;
      quench_tstart_real = time;
      quench_Tstart = T;
    }
    
    if (quench_active && T>0.0){
      //perform_quench(quench_tstart_real, quench_Tstart);
      if (time >= specs.quench_tstart){
	Tnew = quench_Tstart - specs.quench_rate * (time  - quench_tstart_real);
	if (Tnew < 0.0) Tnew = 0.0;
	/* Change velocities: */
	td = 1.0; if (! fp_is_small(T)) td = Tnew / T;
	td = sqrt(td);
#pragma omp parallel for schedule(static)
	for (i=0; i<nat; ++i){
	  vel[i][0] *= td;
	  vel[i][1] *= td;
	  vel[i][2] *= td;
	}
	T = Tnew;
      }
    }

    if (specs.quench_always)
      T = 0.0;


    /*
    // exponent-form:
    // sign + single_digit_integer_part . fractional_part + e + sign + exponent
    // 1 + 1 + n + 1 + 1 + 3 = 7 + n
    // floating-form:
    // sign + integer_part . fractional_part
    // 1 + m + n

  abs(td)<liml ? fmt=fmte : fmt=fmte;
  abs(td)>limu ? fmt=fmte : fmt=fmtf;

  if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td2;
  else                              fout << format(formatf) % td2;

    */

    if (specs_common.report_step && (istep % specs.ndump == 0)){
      printf("time %15.5e  nat %d  Ecoh %10.5f  T %10.5f  Fmax %12.5e  P %12.5e  "
	     "Px Py Pz  %12.5e %12.5e %12.5e    box1 box2 box3  %12.5e %12.5e %12.5e  Vol_at %12.5e\n",
	     time, nat, Ep_tot/nat, T, F_max, P,
	     stresstensor_xyz.elem(0,0), stresstensor_xyz.elem(1,1),stresstensor_xyz.elem(2,2),
	     boxlen[0], boxlen[1], boxlen[2],
	     vol_atom);
      // fflush(stdout);
    }


    if (specs.fixed_geometry) break;
    if (time >= specs.tend) break;

    time += dt;
    istep++;



    // ///////////////////////////////////////////////////////////////////////////////////
    // Some error checking
    // ///////////////////////////////////////////////////////////////////////////////////
    if (specs.use_error_T_gt && T > specs.error_T_gt){
      cout << "ERROR: Temperature T " << T << " K"
	" for compound " << name << " has grown beyond upper allowed value "
	   << specs.error_T_gt << endl;
      throw err_bad_mds;
    }
    if (specs.use_error_dt_lt && dt < specs.error_dt_lt){
      cout << "ERROR: Time step dt " << dt << " fs"
	" for compound " << name << " has shrunk below lower allowed value "
	   << specs.error_dt_lt << endl;
      throw err_bad_mds;
    }
    if (specs.use_error_boxlen_gt){
      int k=0;
      if (boxlen[0] > specs.error_boxlen_gt) ++k;
      if (boxlen[1] > specs.error_boxlen_gt) ++k;
      if (boxlen[2] > specs.error_boxlen_gt) ++k;
      if (k){
	cout << "Boxlen for direction 1: " << boxlen[0] << endl;
	cout << "Boxlen for direction 2: " << boxlen[1] << endl;
	cout << "Boxlen for direction 3: " << boxlen[2] << endl;
	cout << "ERROR: Box length for one or more directions"
	  " for compound " << name << " has grown beyond upper allowed value "
	     << specs.error_boxlen_gt << endl;
	throw err_bad_mds;
      }
    }

    if (! periodic){
      for (int i=0; i<3; ++i){
	if (isnan(boxlen[i])){
	  cout << "ERROR: Box length in direction " << i
	       << " for compound is NaN" << endl;
	  throw err_bad_mds;
	}
      }
    }
    // ///////////////////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////////////////



  }

  /* #############################################################
     #############################################################
     #############################################################

     END OF LOOP OVER TIME

     #############################################################
     #############################################################
     ############################################################# */


  //exit(1);

  {
    ofstream fdump;
    string dumpfn = "mds-lastframe-" + name + ".xyz";
    fdump.open(dumpfn.c_str());
    dumpframe(fdump);
    fdump.close();
  }

  /*
  double fsumx=0,fsumy=0,fsumz=0;
  double rsumx=0,rsumy=0,rsumz=0;
  for (i=0; i<nat; ++i){
    fsumx += frc[i][0];
    fsumy += frc[i][1];
    fsumz += frc[i][2];
    rsumx += pos[i][0];
    rsumy += pos[i][1];
    rsumz += pos[i][2];
  }
  cout << "Net position (A)   : " << rsumx << " " << rsumy << " " << rsumz << endl;
  cout << "CM position (A)    : " << rsumx/nat << " " << rsumy/nat << " " << rsumz/nat << endl;
  cout << "Net force    (eV/A): " << fsumx << " " << fsumy << " " << fsumz << endl;
  */

  if (specs_common.debug_forces){
    double td1;
    string dumpfile = "debugforces-" + name + ".out";
    ofstream fout;

    fout.open(dumpfile.c_str());
    for (i=0; i<nat; ++i){
      td1 = 0;
      for (int k=0; k<3; ++k){
	td1 += (frc[i][k] - frc_num[i][k]) * (frc[i][k] - frc_num[i][k]);
      }
      td1 = sqrt(td1);

      fout << format(" atom %10d  ") % i << format("  abs err %15.10e  ") % td1
	   << format(" x %15.10f y %15.10f z %15.10f  ") % pos[i][0]     % pos[i][1]     % pos[i][2]
	   << format(" ana %15.10f  %15.10f  %15.10f  ") % frc[i][0]     % frc[i][1]     % frc[i][2]
	   << format(" num %15.10f  %15.10f  %15.10f  ") % frc_num[i][0] % frc_num[i][1] % frc_num[i][2]
	   << endl;
    }
    fout.close();


    return;
  }


  /* -------------------------------------------------------------
     Compare initial positions with relaxed positions.
     Ignore changes in the overall box size.
     ------------------------------------------------------------- */
  for (i=0; i<nat; ++i){

    get_coords_cart2skew(pos[i], pos_int_fin[i], -1);

    pos_int_fin[i][0] /= boxlen[0];
    pos_int_fin[i][1] /= boxlen[1];
    pos_int_fin[i][2] /= boxlen[2];


    tv1[0] = pos_int_fin[i][0] - pos_int_ini[i][0];
    while (pbc[0] && (tv1[0] <  -0.5)) tv1[0] += 1;
    while (pbc[0] && (tv1[0] >=  0.5)) tv1[0] -= 1;
    tv1[0] *= boxlen[0];

    tv1[1] = pos_int_fin[i][1] - pos_int_ini[i][1];
    while (pbc[1] && (tv1[1] <  -0.5)) tv1[1] += 1;
    while (pbc[1] && (tv1[1] >=  0.5)) tv1[1] -= 1;
    tv1[1] *= boxlen[1];

    tv1[2] = pos_int_fin[i][2] - pos_int_ini[i][2];
    while (pbc[2] && (tv1[2] <  -0.5)) tv1[2] += 1;
    while (pbc[2] && (tv1[2] >=  0.5)) tv1[2] -= 1;
    tv1[2] *= boxlen[2];


    get_coords_skew2cart(tv1, tv2, -1);
    td = tv2.magn();
    if (i==0 || (i>0 && (td>displ_max))) displ_max = td;
  }



  update_box_geometry();

    
  return;
}
    


void MDSystem::get_virials(int nat, Matrix<double> & W,
			   double mux, double muy, double muz){

  int i,v1,v2;
  double td, mu[3];

  mu[0]=mux;
  mu[1]=muy;
  mu[2]=muz;

  // Virials: W_ij = sum frc_i * dpos_j
  for (v1=0; v1<3; ++v1){
    for (v2=0; v2<3; ++v2){
      W.elem(v1,v2) = 0.0;
      for (i=0; i<nat; ++i){
	W.elem(v1,v2) += virials[i].elem(v1,v2) * mu[v2];
      }
    }
  }
  
  calc_volume();

  // Volume and unit transformation factors:
  td = (1.0/vol) * eVA3_to_GPa;
  for (v1=0; v1<3; ++v1){
    for (v2=0; v2<3; ++v2){
      W.elem(v1,v2) *= td;
    }
  }
      
  return;
}

