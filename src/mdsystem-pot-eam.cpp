



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

#include "exiterrors.hpp"
#include "constants.hpp"
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

using namespace std;
using namespace utils;
using namespace exiterrors;
using boost::format;


double MDSystem::force_EAM(){
  int iat, jat, j;
  double V2_ij, rho_s_ij, rho_p_ij, rho_d_ij;
  double dV2_ij, drho_s_ij, drho_p_ij, drho_d_ij;
  double V2_i,  F_i,  F_s_i, F_p_i, F_d_i, rho_s_i,  rho_p_i , rho_d_i;
  double dF_s_i, dF_p_i, dF_d_i;
  double drsq, rcutsq;
  double td, dx,dy,dz, eps=numeric_limits<double>::epsilon(), dr_ij, dx_ij, dy_ij, dz_ij;
  int ivec;
  Vector<double> pos1(3,0.0), pos2(3,0.0), drvec(3,0.0);
  Vector<double> dF_drho_d(natoms(), 0), dF_drho_p(natoms(), 0), dF_drho_s(natoms(), 0);
  string s1, s2;
  double Ep_tot_local = 0.0;
  int type1, type2;


  if (sys_single_elem){
    type1 = abs(type[0]);
    ivec  = p_potinfo->basepot_vecidx(type1, type1);
  }



  /* Loop over atom 'i'. */
  for (iat=0; iat<natoms(); ++iat){
    Ep[iat]  = 0.0;
    V2_i     = 0.0;
    F_i      = 0.0;
    rho_s_i  = 0.0;
    rho_p_i  = 0.0;
    rho_d_i  = 0.0;

    pos1 = pos[iat];
    s1   = matter[iat];
    type1 = abs(type[iat]);
      
    /* Loop over atom 'j'. */
    for (j=0; j<neighborcollection[iat].size(); ++j){
      jat = neighborcollection[iat][j];
      if (iat==jat) continue;
	  
      V2_ij   = 0;

      pos2 = pos[jat];
      s2   = matter[jat];
      type2 = abs(type[jat]);

      /* Is the interaction correct ? */
      if (! iac_pure_EAM)
	if (p_potinfo->basepot(type1,type2) != "EAM") continue;

      if (! sys_single_elem)
	ivec = p_potinfo->basepot_vecidx(type1,type2);
      if (ivec<0) continue;


      /* --------------------------------------------------------------
	 EAM
	 -------------------------------------------------------------- */
	  
      /* Get the cutoff radius. */
      rcut = p_potinfo->pot_EAM[ivec].rcut();
      rcutsq = rcut * rcut;


      //	    printf("Position of first  atom: %f %f %f\n", x1, y1, z1);
      //	    printf("Position of second atom: %f %f %f\n", x2, y2, z2);

      get_atom_distance_vec(pos1, pos2, drvec);
      drsq = drvec[0]*drvec[0] + drvec[1]*drvec[1] + drvec[2]*drvec[2];


      if (drsq < eps*eps) continue;
      if (drsq > rcutsq) continue;

      dr_ij = sqrt(drsq);

      //	    if (iat==1) printf("dr = %10.5f\n", dr_ij);
      // printf("Distance between atoms Is %12.6f\n", dr_ij);

      /* ------------------------------------------------------------------
	 Calculate the energy of this pair of atoms.
	 ------------------------------------------------------------------ */

      /* Get the pair energy. */
      V2_ij = splint(p_potinfo->pot_EAM[ivec].r,
		     p_potinfo->pot_EAM[ivec].V2,
		     p_potinfo->pot_EAM[ivec].d2_V2,
		     dr_ij);
      V2_ij *= 0.5;



      /* Sum up partial electron densities. */
      rho_s_ij = 0.0;
      rho_p_ij = 0.0;
      rho_d_ij = 0.0;
      if (p_potinfo->pot_EAM[ivec].Nrho_s >0)
	rho_s_ij = splint(p_potinfo->pot_EAM[ivec].r,
			  p_potinfo->pot_EAM[ivec].rho_s,
			  p_potinfo->pot_EAM[ivec].d2_rho_s,
			  dr_ij);
      if (p_potinfo->pot_EAM[ivec].Nrho_p >0)
	rho_p_ij = splint(p_potinfo->pot_EAM[ivec].r,
			  p_potinfo->pot_EAM[ivec].rho_p,
			  p_potinfo->pot_EAM[ivec].d2_rho_p,
			  dr_ij);
      if (p_potinfo->pot_EAM[ivec].Nrho_d >0)
	rho_d_ij = splint(p_potinfo->pot_EAM[ivec].r,
			  p_potinfo->pot_EAM[ivec].rho_d,
			  p_potinfo->pot_EAM[ivec].d2_rho_d,
			  dr_ij);
	    

      //	printf(" dr_ij=%10.6lg  rho_d_ij=%10.6lg\n", dr_ij, rho_d_ij);
      V2_i      += V2_ij;
      rho_s_i   += rho_s_ij;
      rho_p_i   += rho_p_ij;
      rho_d_i   += rho_d_ij;
	    

	    

	    
      /* Go to next atom 'j'. */
    }
    /* End of loop over neighbors 'j'. */



	
    /* The total electron density for atom i is now OK.       */
    /* Now get the embedding energies, which depend on the total electron density. */
    F_s_i = 0.0;
    F_p_i = 0.0;
    F_d_i = 0.0;
    if (p_potinfo->pot_EAM[ivec].Nrho_s >0)
      F_s_i = splint(p_potinfo->pot_EAM[ivec].F_s_rho_s,
		     p_potinfo->pot_EAM[ivec].F_s,
		     p_potinfo->pot_EAM[ivec].d2_F_s,
		     rho_s_i);
    if (p_potinfo->pot_EAM[ivec].Nrho_p >0)
      F_p_i = splint(p_potinfo->pot_EAM[ivec].F_p_rho_p,
		     p_potinfo->pot_EAM[ivec].F_p,
		     p_potinfo->pot_EAM[ivec].d2_F_p,
		     rho_p_i);
    if (p_potinfo->pot_EAM[ivec].Nrho_d >0)
      F_d_i = splint(p_potinfo->pot_EAM[ivec].F_d_rho_d,
		     p_potinfo->pot_EAM[ivec].F_d,
		     p_potinfo->pot_EAM[ivec].d2_F_d,
		     rho_d_i);
    /* Total embedding energy: */
    F_i     = F_s_i + F_p_i + F_d_i;
    // printf("Total embedding energy for atom %ld is %20.10lg\n", iat, F_i);
	
    /* -------------------------------------------------------------------
       EAM energy:
       ------------------------------------------------------------------- */
    Ep[iat]  = F_i + V2_i;
    Ep_tot_local  += Ep[iat];



    // printf("Atom i: rho_d=%10.6lg\n", rho_d_i);
    // printf("Atom i: 1/2 V2=%10.6lg\n", V2_i);
    // printf("Atom i: F_s=%10.6lg   F_p=%10.6lg   F_d=%10.6lg\n", F_s_i, F_p_i, F_d_i);
      
    // printf( "Energy of atom i is %10.6lg\n", Epot_i);



    if (get_pot_force){

      /* -------------------------------------------------------------------
	 Get the derivative of F with respect to the total electron density rho.
	 ------------------------------------------------------------------- */
      dF_s_i = 0.0;
      dF_p_i = 0.0;
      dF_d_i = 0.0;
      if (p_potinfo->pot_EAM[ivec].Nrho_s > 0)
	dF_s_i = splint_dy(p_potinfo->pot_EAM[ivec].F_s_rho_s,
			   p_potinfo->pot_EAM[ivec].F_s,
			   p_potinfo->pot_EAM[ivec].d2_F_s,
			   rho_s_i);
      if (p_potinfo->pot_EAM[ivec].Nrho_p > 0)
	dF_p_i = splint_dy(p_potinfo->pot_EAM[ivec].F_p_rho_p,
			   p_potinfo->pot_EAM[ivec].F_p,
			   p_potinfo->pot_EAM[ivec].d2_F_p,
			   rho_p_i);
      if (p_potinfo->pot_EAM[ivec].Nrho_d > 0)
	dF_d_i = splint_dy(p_potinfo->pot_EAM[ivec].F_d_rho_d,
			   p_potinfo->pot_EAM[ivec].F_d,
			   p_potinfo->pot_EAM[ivec].d2_F_d,
			   rho_d_i);
      dF_drho_d[iat] += dF_d_i;
      dF_drho_p[iat] += dF_p_i;
      dF_drho_s[iat] += dF_s_i;
    }
      
	
    /* Go to next atom 'iat'. */
  }





    
  if (! get_pot_force) return Ep_tot_local;






  /* -------------------------------------------------------------------
     Add the many-body contributions from the embedding energy to the forces.
     ------------------------------------------------------------------- */
  for (iat=0; iat<natoms(); ++iat){
    pos1 = pos[iat];
    s1   = matter[iat];
    type1 = abs(type[iat]);


    /* Now loop over neighbors, calculating the derivatives of
       the partial electron densities. */
    //	printf("Looping over %ld neighbors of atom %ld\n", neighbors[iat].nn, iat);
      
    for (j=0; j<neighborcollection[iat].size(); ++j){
      jat = neighborcollection[iat][j];
      if (iat==jat) continue;

      pos2 = pos[jat];
      s2   = matter[jat];
      type2 = abs(type[jat]);


      /* Is the interaction correct ? */
      if (! iac_pure_EAM)
	if (p_potinfo->basepot(type1,type2) != "EAM") continue;
      
      if (! sys_single_elem)
	ivec = p_potinfo->basepot_vecidx(type1,type2);
      if (ivec<0) continue;


      /* --------------------------------------------------------------
	 EAM
	 -------------------------------------------------------------- */
	
      /* Get the cutoff radius. */
      rcut = p_potinfo->pot_EAM[ivec].rcut();
      rcutsq = rcut * rcut;


      //	    printf("Position of first  atom: %f %f %f\n", x1, y1, z1);
      //	    printf("Position of second atom: %f %f %f\n", x2, y2, z2);

      get_atom_distance_vec(pos1, pos2, drvec);
      dr_ij = drvec.magn();

      if (dr_ij < eps) continue;
      if (dr_ij > rcut) continue;


      drho_s_ij = 0.0;
      drho_p_ij = 0.0;
      drho_d_ij = 0.0;
      if (p_potinfo->pot_EAM[ivec].Nrho_s >0)
	drho_s_ij = splint_dy(p_potinfo->pot_EAM[ivec].r,
			      p_potinfo->pot_EAM[ivec].rho_s,
			      p_potinfo->pot_EAM[ivec].d2_rho_s,
			      dr_ij);
      if (p_potinfo->pot_EAM[ivec].Nrho_p >0)
	drho_p_ij = splint_dy(p_potinfo->pot_EAM[ivec].r,
			      p_potinfo->pot_EAM[ivec].rho_p,
			      p_potinfo->pot_EAM[ivec].d2_rho_p,
			      dr_ij);
      if (p_potinfo->pot_EAM[ivec].Nrho_d >0)
	drho_d_ij = splint_dy(p_potinfo->pot_EAM[ivec].r,
			      p_potinfo->pot_EAM[ivec].rho_d,
			      p_potinfo->pot_EAM[ivec].d2_rho_d,
			      dr_ij);


      dx_ij = drvec[0];
      dy_ij = drvec[1];
      dz_ij = drvec[2];
	    
      dx  = dx_ij;   dy  = dy_ij;   dz  = dz_ij;
      dx /= dr_ij;   dy /= dr_ij;   dz /= dr_ij;

      /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	 Pair part:
	 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
	
      dV2_ij = splint_dy(p_potinfo->pot_EAM[ivec].r,
			 p_potinfo->pot_EAM[ivec].V2,
			 p_potinfo->pot_EAM[ivec].d2_V2,
			 dr_ij);

      frc[iat][0] += - 2  * 0.5 * dV2_ij * dx ;
      frc[iat][1] += - 2  * 0.5 * dV2_ij * dy ;
      frc[iat][2] += - 2  * 0.5 * dV2_ij * dz ;

      virials[iat].elem(0,0) += -dx_ij * (- 1  * 0.5 * dV2_ij * -dx ) ;
      virials[iat].elem(1,1) += -dy_ij * (- 1  * 0.5 * dV2_ij * -dy ) ;
      virials[iat].elem(2,2) += -dz_ij * (- 1  * 0.5 * dV2_ij * -dz ) ;




      /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	 Many-body part:
	 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

      td =  - (dF_drho_s[iat] * drho_s_ij +
		 dF_drho_p[iat] * drho_p_ij +
		 dF_drho_d[iat] * drho_d_ij);
	    
      frc[iat][0] += td * dx ;
      frc[iat][1] += td * dy ;
      frc[iat][2] += td * dz ;


      td =  - (dF_drho_s[jat] * drho_s_ij +
		 dF_drho_p[jat] * drho_p_ij +
		 dF_drho_d[jat] * drho_d_ij);
	
      frc[iat][0] += td * dx ;
      frc[iat][1] += td * dy ;
      frc[iat][2] += td * dz ;

      virials[iat].elem(0,0) += -dx_ij * td * -dx;
      virials[iat].elem(1,1) += -dy_ij * td * -dy;
      virials[iat].elem(2,2) += -dz_ij * td * -dz;


      /* Go to next atom 'jat'. */
    }
    /* End of loop over neighbors 'j'. */
    /* Go to next atom 'iat'. */
  }


  //printf("Done.\n");


  return Ep_tot_local;
}


