


#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
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
#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mtwister.hpp"
#include "param-pot.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "potinfo.hpp"
//#include "specs-fit-prop-pot.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;


using namespace utils;
using namespace constants;
using boost::format;







double MDSystem::force_ABOP(){
  int i;

  int nat = natoms();
  double Ep_tot_local = 0.0;

  int se_typei = itype[0]; //elem.name2idx( matter[0] );
  int se_ivecij = p_potinfo->basepot_vecidx(se_typei, se_typei);
  int se_ivec_reppot = p_potinfo->reppot_vecidx(se_typei, se_typei);
  bool se_use_reppot = p_potinfo->use_reppot(se_typei, se_typei);


  double se_alphaijk = 0.0;
  double se_omegaijk = 1.0;
  double se_twomuik  = 0.0;
  if (p_potinfo->use_abop_alpha.elem(se_typei, se_typei, se_typei))
    se_alphaijk = p_potinfo->abop_alpha.elem(se_typei, se_typei, se_typei);
  if (p_potinfo->use_abop_omega.elem(se_typei, se_typei, se_typei))
    se_omegaijk = p_potinfo->get_abop_omega(se_typei, se_typei, se_typei);
  if (p_potinfo->use_abop_2mu.elem(se_typei, se_typei))
    se_twomuik = p_potinfo->abop_2mu.elem(se_typei, se_typei);

  


  /*
  std::cout << "sys_single_elem: " << sys_single_elem << std::endl;
  std::cout << "iac_pure_ABOP  : " << iac_pure_ABOP << std::endl;
  exit(0);
  */




  int ij,j,ik,k,p,q,s,is;
  int typei, typej, typek,types;
  double D0ij,r0ij,betaij,Sij;
  double Rij,Dij,Rik,Dik;
  double gammaik,cik,dik,hik,r0ik;
  double fcij,rij,rcutij;
  double fcik,rik,rcutik;
  double ris,rcutis;
  double V1, dV1, fermi, dfermi, bfermi, rfermi;
  double VRij,VAij,bij,Chiij,gijk,c2,d2,cost,hcost,hcost2,alphaijk,twomuik;
  double pij;

  //double eps = std::numeric_limits<double>::epsilon();
  double VRij_r, VAij_r, dVRij_r, dVAij_r;
  double dVRij, dVAij, dgijk;
  double dfcij, dfcik, threebodyfactor;

  Vector3<double> dposij,dposik,dposjk,dposis,dposks;
  Vector3<double> dcost_i, dcost_j, dcost_k;
  Vector3<double> dcost_ij, dcost_ji, dcost_ik, dcost_ki;
  Vector3<double> dgijk_i, dgijk_j, dgijk_k;
  Vector3<double> dgijk_ij, dgijk_ji, dgijk_ik, dgijk_ki;
  Vector3<double> frci, frcj, frck;

  double td, td1,td2;
  int ivecij, ivecik, ivecjk, ivecis, ivecks, ivec_reppot, Nr;
  double Epij;

  double frcij;
  Vector3<double> frc_ij, frc_ik, frc_jk, frc_is, frc_ks;
  Vector3<double> frc_ji, frc_ki, frc_kj, frc_si, frc_sk;

  double F1, F2, dF1, dF2;
  Vector3<double> dF1_i, dF1_j, dF1_k;
  Vector3<double> dF2_i, dF2_j, dF2_k;
  Vector3<double> dF1_ij, dF1_ji, dF1_ik, dF1_ki;
  Vector3<double> dF2_ij, dF2_ji, dF2_ik, dF2_ki;


  bool pair_ij_tersoff;
  bool pair_ik_tersoff;

  bool pair_ij_perriot_cut;
  bool pair_ij_perriot_scr;
  bool pair_ik_perriot_cut;
  bool pair_ik_perriot_scr;
  bool pair_jk_perriot_cut;
  bool pair_jk_perriot_scr;

  double rminij, rmaxij, rminik, rmaxik;



  int nelem = elem.nelem();
  Matrix<double> rcut_all(nelem,nelem);
  Matrix<CutoffScreeningPair> rcs_all(nelem,nelem);
  Matrix<Potential_ABOPPair> abop_params_all(nelem,nelem);
  Matrix<std::string> basepot_all(nelem,nelem);
  Matrix<int> basepot_vecidx_all(nelem,nelem);
  Matrix<int> use_reppot_all(nelem,nelem);
  Matrix<int> reppot_vecidx_all(nelem,nelem);

  for (int i=0; i<nelem; ++i){
    for (int j=0; j<nelem; ++j){
      ivecij = p_potinfo->basepot_vecidx(i,j);
      basepot_vecidx_all.elem(i,j)=ivecij;
      if (ivecij<0) rcut_all.elem(i,j)=0.0;

      if (p_potinfo->basepot(i,j) != "ABOP") continue;
      basepot_all.elem(i,j)="ABOP";

      rcut_all.elem(i,j)=p_potinfo->pot_ABOP[ivecij].rcut();

      if (p_potinfo->pot_ABOP[ivecij].rcs.name=="tersoff"
	  &&
	  p_potinfo->pot_ABOP[ivecij].rcs.mode=="cut"){
	rcs_all.elem(i,j).tersoff = true;
	rcs_all.elem(i,j).R    = p_potinfo->pot_ABOP[ivecij].get_parval("R");
	rcs_all.elem(i,j).D    = p_potinfo->pot_ABOP[ivecij].get_parval("D");
	rcs_all.elem(i,j).rcut = p_potinfo->pot_ABOP[ivecij].rcut();
      }
      else if (p_potinfo->pot_ABOP[ivecij].rcs.name=="perriot"
	       &&
	       p_potinfo->pot_ABOP[ivecij].rcs.mode=="cut"){
	rcs_all.elem(i,j).perriot_cut = true;
	rcs_all.elem(i,j).prmin = p_potinfo->pot_ABOP[ivecij].get_parval("prmin");
	rcs_all.elem(i,j).prmax = p_potinfo->pot_ABOP[ivecij].get_parval("prmax");
	rcs_all.elem(i,j).rcut  = p_potinfo->pot_ABOP[ivecij].rcut();
      }
      else if (p_potinfo->pot_ABOP[ivecij].rcs.name=="perriot"
	       &&
	       p_potinfo->pot_ABOP[ivecij].rcs.mode=="scr"){
	rcs_all.elem(i,j).perriot_scr = true;

	rcs_all.elem(i,j).pn    = p_potinfo->pot_ABOP[ivecij].get_parval("pn");
	rcs_all.elem(i,j).pm    = p_potinfo->pot_ABOP[ivecij].get_parval("pm");
	rcs_all.elem(i,j).rcut  = p_potinfo->pot_ABOP[ivecij].get_parval("prcut");
	rcs_all.elem(i,j).prmin = p_potinfo->pot_ABOP[ivecij].get_parval("prmin");
	rcs_all.elem(i,j).prmax = p_potinfo->pot_ABOP[ivecij].get_parval("prmax");
      }
      else
	rcs_all.elem(i,j).rcut = 0.0;

      abop_params_all.elem(i,j).D0 = p_potinfo->pot_ABOP[ivecij].get_parval("D0");
      abop_params_all.elem(i,j).r0 = p_potinfo->pot_ABOP[ivecij].get_parval("r0");
      abop_params_all.elem(i,j).beta = p_potinfo->pot_ABOP[ivecij].get_parval("beta");
      abop_params_all.elem(i,j).S = p_potinfo->pot_ABOP[ivecij].get_parval("S");
      abop_params_all.elem(i,j).p = p_potinfo->pot_ABOP[ivecij].get_parval("p");
      abop_params_all.elem(i,j).gamma = p_potinfo->pot_ABOP[ivecij].get_parval("gamma");
      abop_params_all.elem(i,j).c = p_potinfo->pot_ABOP[ivecij].get_parval("c");
      abop_params_all.elem(i,j).d = p_potinfo->pot_ABOP[ivecij].get_parval("d");
      abop_params_all.elem(i,j).h = p_potinfo->pot_ABOP[ivecij].get_parval("h");
      abop_params_all.elem(i,j).rfermi = p_potinfo->pot_ABOP[ivecij].get_parval("rfermi");
      abop_params_all.elem(i,j).bfermi = p_potinfo->pot_ABOP[ivecij].get_parval("bfermi");

      use_reppot_all.elem(i,j) = p_potinfo->use_reppot(i,j);
      reppot_vecidx_all.elem(i,j) = p_potinfo->reppot_vecidx(i,j);
    }
  }

  Matrix3<double> abop_omega(nelem,nelem,nelem);
  
  for (int i=0; i<nelem; ++i){
    for (int j=0; j<nelem; ++j){
      for (int k=0; k<nelem; ++k){

	if (p_potinfo->use_abop_alpha.elem(i,j,k)){
	  if (p_potinfo->use_abop_omega.elem(i,j,k)){
	    abop_omega.elem(i,j,k)=p_potinfo->get_abop_omega(i,j,k);
	  }
	}

      }
    }
  }





  double Kij, Kik;
  Vector<double> vecKik(0);
  int ijkp;
  int ijp;

  
  Vector<bool>   hv0_perriot_scr_K_frc_ij;
  Vector<double> hv1_perriot_scr_K_frc_ij, hv2_perriot_scr_K_frc_ij, hv3_perriot_scr_K_frc_ij;
  Vector<double> hv4_perriot_scr_K_frc_ij, hv5_perriot_scr_K_frc_ij;
  Vector< Vector3<double> > hv6_perriot_scr_K_frc_ij, hv7_perriot_scr_K_frc_ij;

  Vector< Vector<bool> >   hv0_perriot_scr_K_frc_ik;
  Vector< Vector<double> > hv1_perriot_scr_K_frc_ik, hv2_perriot_scr_K_frc_ik, hv3_perriot_scr_K_frc_ik;
  Vector< Vector<double> > hv4_perriot_scr_K_frc_ik, hv5_perriot_scr_K_frc_ik;
  Vector< Vector< Vector3<double> > > hv6_perriot_scr_K_frc_ik, hv7_perriot_scr_K_frc_ik;


  int ijh, ijkh;


  // iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  // Loop over atoms i
  // iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

  //#pragma omp parallel for reduction(+:Ep_tot_local)

  // Allocate helper vectors once. We are mostly dealing with equilibrium structures
  // where the number of neighbors of each atom is on average the same. In e.g.
  // collision cascades some atoms may have many neighbors, and some very few.
  // In those cases allocation once like here might waste memory.
  int nijmax=0, ns;
  for (i=0; i<nat; ++i){
    ns = neighborcollection[i].size();
    if (ns>=nijmax) nijmax=ns;
  }
  Vector<bool>   abop_is_neigh_ij(nijmax, false);
  Vector<double> abop_fc_ij(nijmax, 0.0);
  Vector<double> abop_dfc_ij(nijmax, 0.0);
  Vector< Vector3<double> > abop_dpos_ij(nijmax, Vector3<double>(0) );



  for (i=0; i<nat; ++i){
    ijp = 0;
    typei = itype[i]; //elem.name2idx( matter[i] );

    // jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
    // Loop over j
    // jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj

    // Initialize:
    for (ij=0; ij<nijmax; ij++){
      abop_is_neigh_ij[ij]=false;
      abop_fc_ij[ij]=0.0;
      abop_dfc_ij[ij]=0.0;
      abop_dpos_ij[ij][0]=0;
      abop_dpos_ij[ij][1]=0;
      abop_dpos_ij[ij][2]=0;
    }

    ijh = -1;
    for (ij=0; ij<neighborcollection[i].size(); ij++){
      ijh++;
      abop_is_neigh_ij[ijh] = false;

      j = neighborcollection[i][ij];
      typej = itype[j]; //elem.name2idx( matter[j] );

      if (i==j) aborterror("ERROR: Neighbor of atom i is atom i itself!");
      if (! iac_pure_ABOP) if (basepot_all.elem(typei,typej) != "ABOP") continue;
      ivecij = se_ivecij;
      if (! sys_single_elem) ivecij = basepot_vecidx_all.elem(typei, typej);
      if (ivecij<0) continue;

      /* ############################ cutoff/screening ############################ */
      rcutij = rcs_all.elem(typei,typej).rcut;
      get_atom_distance_vec(pos[i], pos[j], dposij);
      rij = dposij.magn();
      if (rij > rcutij) continue;
      /* ################################################################### */

      abop_is_neigh_ij[ijh] = true;
      abop_dpos_ij[ijh]     = dposij;

      /* ############################ cutoff/screening ############################ */
      pair_ij_tersoff     = rcs_all.elem(typei,typej).tersoff;
      pair_ij_perriot_cut = rcs_all.elem(typei,typej).perriot_cut;
      pair_ij_perriot_scr = rcs_all.elem(typei,typej).perriot_scr;
      Rij = rcs_all.elem(typei,typej).R;
      Dij = rcs_all.elem(typei,typej).D;
      rminij = rcs_all.elem(typei,typej).prmin;
      rmaxij = rcs_all.elem(typei,typej).prmax;
      rcutij = rcs_all.elem(typei,typej).rcut;
      /* ################################################################### */
      
      D0ij   = abop_params_all.elem(typei,typej).D0;
      r0ij   = abop_params_all.elem(typei,typej).r0;
      betaij = abop_params_all.elem(typei,typej).beta;
      Sij    = abop_params_all.elem(typei,typej).S;
      /*
      std::cout << "D0ij " << D0ij << " r0ij " << r0ij << " betaij " << betaij << " Sij " << Sij
	   << " gammaij " << gammaij << " cij " << cij << " dij " << dij << " hij " << hij
	   << " Rij " << Rij << " Dij " << Dij << std::endl;
      */
      /* ############################ cutoff/screening ############################ */
      /* For the ij pair we can have either a cutoff function (Tersoff/Perriot) or
	 a screening (Perriot).
       */
      fcij  = 1.0;
      dfcij = 0.0;
      abop_fc_ij[ijh]  = fcij;
      abop_dfc_ij[ijh] = dfcij;

      if (pair_ij_tersoff){
	if (rij < Rij-Dij){
	  fcij = 1.0;  dfcij = 0.0;
	  abop_fc_ij[ijh]  = fcij;
	  abop_dfc_ij[ijh] = dfcij;
	}
	else if (rij > Rij+Dij){
	  fcij = 0.0;  dfcij = 0.0;
	  abop_fc_ij[ijh]  = fcij;
	  abop_dfc_ij[ijh] = dfcij;
	  continue;
	}
	else {
	  td = 0.5*PI*(rij-Rij)/Dij;
	  fcij  = 0.5 * (1.0 - sin( td ));
	  dfcij = -0.5 * 0.5*PI/Dij * cos( td );
	  abop_fc_ij[ijh]  = fcij;
	  abop_dfc_ij[ijh] = dfcij;
	}
      }
      else if (pair_ij_perriot_cut){
	if (rij < rminij){
	  fcij  = 1.0;  dfcij = 0.0;
	  abop_fc_ij[ijh]  = fcij;
	  abop_dfc_ij[ijh] = dfcij;
	}
	else if (rij > rmaxij){
	  fcij  = 0.0;  dfcij = 0.0;
	  abop_fc_ij[ijh]  = fcij;
	  abop_dfc_ij[ijh] = dfcij;
	  continue;
	}
	else {
	  double invdrc = 1.0/(rmaxij - rminij);
	  double a = (rij - rminij) * invdrc;

	  double aa = a*a;
	  double aaa = a*aa;

	  double td1 = 6*aa - 15*a + 10;
	  fcij = 1.0 - aaa * td1;
	  dfcij = ( -3*aa * td1 - aaa * (12*a - 15) ) * invdrc;
	  abop_fc_ij[ijh]  = fcij;
	  abop_dfc_ij[ijh] = dfcij;
	}
      }
      /* ################################################################### */
    } // end of loop over j
    // 'Is neighbor?', dpos, fc/dfc saved into vectors for all neighbors ijh.


    // jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
    // Loop over j
    // jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj

    ijh = -1;
    for (ij=0; ij<neighborcollection[i].size(); ij++){
      ijh++;

      if ( ! abop_is_neigh_ij[ijh]) continue;

      j = neighborcollection[i][ij];
      typej = itype[j]; //elem.name2idx( matter[j] );

      /*
      Vector3<double> tv1;
      get_atom_distance_vec(pos[i], pos[j], tv1);
      tv1 = tv1 - dposij;
      std::cout << "ijh: " << ijh << " Difference in stored and direct computed ij distance vectors: " << tv1.magn() << std::endl;
      */

      dposij = abop_dpos_ij[ijh];
      rij    = dposij.magn();
      fcij   = abop_fc_ij[ijh];
      dfcij  = abop_dfc_ij[ijh];


      /* ############################ cutoff/screening ############################ */
      pair_ij_tersoff     = rcs_all.elem(typei,typej).tersoff;
      pair_ij_perriot_cut = rcs_all.elem(typei,typej).perriot_cut;
      pair_ij_perriot_scr = rcs_all.elem(typei,typej).perriot_scr;
      Rij = rcs_all.elem(typei,typej).R;
      Dij = rcs_all.elem(typei,typej).D;
      rminij = rcs_all.elem(typei,typej).prmin;
      rmaxij = rcs_all.elem(typei,typej).prmax;
      rcutij = rcs_all.elem(typei,typej).rcut;
      /* ################################################################### */
      
      D0ij   = abop_params_all.elem(typei,typej).D0;
      r0ij   = abop_params_all.elem(typei,typej).r0;
      betaij = abop_params_all.elem(typei,typej).beta;
      Sij    = abop_params_all.elem(typei,typej).S;

      //std::cout << "rij " << rij << " rcutij " << rcutij << " fcij " << fcij << std::endl;

      td1 = - betaij*sqrt(2.0*Sij);
      td2 = - betaij*sqrt(2.0/Sij);

      VRij =       D0ij/(Sij-1.0) * exp( td1 * (rij - r0ij) );
      VAij = Sij * D0ij/(Sij-1.0) * exp( td2 * (rij - r0ij) );

      dVRij = td1 * VRij;
      dVAij = td2 * VAij;


      // ################################################################
      // reppot
      // ################################################################

      if ( (  sys_single_elem && se_use_reppot) ||
	   (! sys_single_elem && use_reppot_all.elem(typei, typej)) ){

	ivec_reppot = se_ivec_reppot;
	if (! sys_single_elem)
	  ivec_reppot = reppot_vecidx_all.elem(typei, typej);
	//ivec_reppot = p_potinfo->reppot_vecidx(s1,s2);
	if (ivec_reppot<0) break;

	Nr = p_potinfo->pot_Reppot[ivec_reppot].r_rep.size();
	td = p_potinfo->pot_Reppot[ivec_reppot].r_rep[Nr-1];

	bfermi = abop_params_all.elem(typei,typej).bfermi;
	rfermi = abop_params_all.elem(typei,typej).rfermi;
	//bfermi = p_potinfo->pot_Reppot[ivec_reppot].bfermi;
	//rfermi = p_potinfo->pot_Reppot[ivec_reppot].rfermi;
	  
	td = exp(- bfermi * (rij - rfermi));
	fermi  = 1.0/(1.0 + td);
	dfermi = fermi*fermi * bfermi * td;

	if (rij < td){
	  V1 = splint(p_potinfo->pot_Reppot[ivec_reppot].r_rep,
		      p_potinfo->pot_Reppot[ivec_reppot].V_rep,
		      p_potinfo->pot_Reppot[ivec_reppot].d2_V_rep,
		      rij);
	  
	  dV1 = splint_dy(p_potinfo->pot_Reppot[ivec_reppot].r_rep,
			  p_potinfo->pot_Reppot[ivec_reppot].V_rep,
			  p_potinfo->pot_Reppot[ivec_reppot].d2_V_rep,
			  rij);
	}
	else {
	  // reppot ended already, assume 0.0
	  V1  = 0.0;
	  dV1 = 0.0;
	}


	  
	//              call splinereppot(ra(ij),V,df,bf(typei,typej),&
	//                   & rf(typei,typej),typei,typej,fermi,dfermi)
	//              !if (ra(ij) < 2.2) print '(A,4F10.3)','rep2',ra(ij),V,df,fermi
	//
	//              df=df+dfermi*(cpair+cmany)
	  
	
	VRij_r = (1.0 - fermi) * V1 + fermi * VRij;
	VAij_r = (1.0 - fermi) * V1 + fermi * VAij;
	  
	dVRij_r = dV1 * (1-fermi) - V1 * dfermi  +  dVRij * fermi + VRij * dfermi;
	dVAij_r = dV1 * (1-fermi) - V1 * dfermi  +  dVAij * fermi + VAij * dfermi;
	  
	VRij = VRij_r;
	VAij = VAij_r;
	  
	dVRij = dVRij_r;
	dVAij = dVAij_r;
      }

    


      // ################################################################
      // Screening factor Kij:
      // ################################################################


    // Kij: prefactor = - 0.5 * (VRij - bij * VAij) ;
    // Kik: prefactor = threebodyfactor * fcik * gijk * F1 * F2 ;



      // ###################### Perriot screening factor K ######################





      //      vecKij.resize( neighborcollection[i].size() );
      Kij = 1.0;
      if (pair_ij_perriot_scr){
	int nnn = neighborcollection[i].size();
	if (hv0_perriot_scr_K_frc_ij.size()<=nnn) hv0_perriot_scr_K_frc_ij.resize(nnn);
	if (hv1_perriot_scr_K_frc_ij.size()<=nnn) hv1_perriot_scr_K_frc_ij.resize(nnn);
	if (hv2_perriot_scr_K_frc_ij.size()<=nnn) hv2_perriot_scr_K_frc_ij.resize(nnn);
	if (hv3_perriot_scr_K_frc_ij.size()<=nnn) hv3_perriot_scr_K_frc_ij.resize(nnn);
	if (hv4_perriot_scr_K_frc_ij.size()<=nnn) hv4_perriot_scr_K_frc_ij.resize(nnn);
	if (hv5_perriot_scr_K_frc_ij.size()<=nnn) hv5_perriot_scr_K_frc_ij.resize(nnn);
	if (hv6_perriot_scr_K_frc_ij.size()<=nnn) hv6_perriot_scr_K_frc_ij.resize(nnn);
	if (hv7_perriot_scr_K_frc_ij.size()<=nnn) hv7_perriot_scr_K_frc_ij.resize(nnn);
	force_ABOP_perriot_K(i, j, Kij, dposij, rcut_all, rcs_all, basepot_all, basepot_vecidx_all,
			     hv0_perriot_scr_K_frc_ij,
			     hv1_perriot_scr_K_frc_ij, hv2_perriot_scr_K_frc_ij, hv3_perriot_scr_K_frc_ij,
			     hv4_perriot_scr_K_frc_ij, hv5_perriot_scr_K_frc_ij,
			     hv6_perriot_scr_K_frc_ij, hv7_perriot_scr_K_frc_ij);
      }

      // ###################### Perriot screening factor K ######################
      


      // ################################################################
      // Bond order factor bij:
      // ################################################################

      // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      // Loop over k
      // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

      ijkp = 0;
      ijkh = -1;
      Chiij = 0.0;
      for (ik=0; ik<neighborcollection[i].size(); ik++){
	ijkh++;
	
	/*
	  {Neighbors k} + atom j = {Neighbors j} !!!
	*/
	if ( ! abop_is_neigh_ij[ijkh]) continue;

	k = neighborcollection[i][ik];
	if (k==j) continue;
	

	typek = itype[k]; //elem.name2idx( matter[k] );
	//s3 = matter[k]; //p_potinfo->elem.idx2name(typek);

	/* ############################ cutoff/screening ############################ */
	pair_ik_tersoff     = rcs_all.elem(typei,typek).tersoff;
	pair_ik_perriot_cut = rcs_all.elem(typei,typek).perriot_cut;
	pair_ik_perriot_scr = rcs_all.elem(typei,typek).perriot_scr;
	Rik = rcs_all.elem(typei,typek).R;
	Dik = rcs_all.elem(typei,typek).D;
	rminik = rcs_all.elem(typei,typek).prmin;
	rmaxik = rcs_all.elem(typei,typek).prmax;
	rcutik = rcs_all.elem(typei,typek).rcut;
	/* ################################################################### */

	dposik = abop_dpos_ij[ijkh];
	/*
	Vector3<double> tv1;
	get_atom_distance_vec(pos[i], pos[k], tv1);
	tv1 = tv1 - dposik;
	std::cout << "ijkh: " << ijkh << " Difference in stored and direct computed ik distance vectors: " << tv1.magn() << std::endl;
	*/

	rik = dposik.magn();
	if (rik > rcutik) continue;

	/*
	get_atom_distance_vec(pos[i], pos[k], dposik);
	rik = dposik.magn();
	if (rik > rcutik) continue;
	*/
	/* ################################################################### */

	

	r0ik = abop_params_all.elem(typei,typek).r0;

	gammaik = abop_params_all.elem(typei,typek).gamma;
	cik = abop_params_all.elem(typei,typek).c;
	dik = abop_params_all.elem(typei,typek).d;
	hik = abop_params_all.elem(typei,typek).h;
	/*
	std::cout << "D0ik " << D0ik << " r0ik " << r0ik << " betaik " << betaik << " Sik " << Sik
	     << " gammaik " << gammaik << " cik " << cik << " dik " << dik << " hik " << hik
	     << " Rik " << Rik << " Dik " << Dik << std::endl;
	*/

	/* ############################ cutoff/screening ############################ */
	fcik = abop_fc_ij[ijkh];
	dfcik = abop_dfc_ij[ijkh];
	
	/*
	fcik  = 1.0;
	dfcik = 0.0;
	*/
	
	int nnn = neighborcollection[i].size();
	if (vecKik.size() <= nnn) vecKik.resize(nnn);
	vecKik[ijkp] = 1.0;

	if (pair_ik_perriot_scr){
	  if (hv0_perriot_scr_K_frc_ik.size() <= nnn) hv0_perriot_scr_K_frc_ik.resize(nnn);
	  if (hv1_perriot_scr_K_frc_ik.size() <= nnn) hv1_perriot_scr_K_frc_ik.resize(nnn);
	  if (hv2_perriot_scr_K_frc_ik.size() <= nnn) hv2_perriot_scr_K_frc_ik.resize(nnn);
	  if (hv3_perriot_scr_K_frc_ik.size() <= nnn) hv3_perriot_scr_K_frc_ik.resize(nnn);
	  if (hv4_perriot_scr_K_frc_ik.size() <= nnn) hv4_perriot_scr_K_frc_ik.resize(nnn);
	  if (hv5_perriot_scr_K_frc_ik.size() <= nnn) hv5_perriot_scr_K_frc_ik.resize(nnn);
	  if (hv6_perriot_scr_K_frc_ik.size() <= nnn) hv6_perriot_scr_K_frc_ik.resize(nnn);
	  if (hv7_perriot_scr_K_frc_ik.size() <= nnn) hv7_perriot_scr_K_frc_ik.resize(nnn);

	  if (hv0_perriot_scr_K_frc_ik[ijkp].size() <= nnn) hv0_perriot_scr_K_frc_ik[ijkp].resize(nnn);
	  if (hv1_perriot_scr_K_frc_ik[ijkp].size() <= nnn) hv1_perriot_scr_K_frc_ik[ijkp].resize(nnn);
	  if (hv2_perriot_scr_K_frc_ik[ijkp].size() <= nnn) hv2_perriot_scr_K_frc_ik[ijkp].resize(nnn);
	  if (hv3_perriot_scr_K_frc_ik[ijkp].size() <= nnn) hv3_perriot_scr_K_frc_ik[ijkp].resize(nnn);
	  if (hv4_perriot_scr_K_frc_ik[ijkp].size() <= nnn) hv4_perriot_scr_K_frc_ik[ijkp].resize(nnn);
	  if (hv5_perriot_scr_K_frc_ik[ijkp].size() <= nnn) hv5_perriot_scr_K_frc_ik[ijkp].resize(nnn);
	  if (hv6_perriot_scr_K_frc_ik[ijkp].size() <= nnn) hv6_perriot_scr_K_frc_ik[ijkp].resize(nnn);
	  if (hv7_perriot_scr_K_frc_ik[ijkp].size() <= nnn) hv7_perriot_scr_K_frc_ik[ijkp].resize(nnn);
	}


	/*	
	if (pair_ik_tersoff){
	  if (rik < Rik-Dik){
	    fcik  = 1.0;  dfcik = 0.0;
	  }
	  else if (rik > rcutik){
	    fcik  = 0.0;  dfcik = 0.0;
	    continue;
	  }
	  else {
	    td = 0.5*PI*(rik-Rik)/Dik;

	    fcik  = 0.5 * (1.0 - sin( td ));
	    dfcik = -0.5 * 0.5*PI/Dik * cos( td );
	  }
	}
	else if (pair_ik_perriot_cut || pair_ik_perriot_scr){
	  if (rik < rminik){
	    fcik  = 1.0;  dfcik = 0.0;
	  }
	  else if (rik > rmaxik){
	    fcik  = 0.0;  dfcik = 0.0;
	    continue;
	  }
	  else {
	    double invdrc = 1.0/(rmaxik - rminik);
	    double a = (rik - rminik) * invdrc;

	    double aa = a*a;
	    double aaa = a*aa;

	    double td1 = 6*aa - 15*a + 10;
	    fcik = 1.0 - aaa * td1;
	    dfcik = ( -3*aa * td1 - aaa * (12*a - 15) ) * invdrc;
	  }
	*/


	if (pair_ik_perriot_cut || pair_ik_perriot_scr){
	  // ###################### Perriot screening factor K ######################
	  vecKik[ijkp] = 1.0;
	  if (pair_ik_perriot_scr)
	    force_ABOP_perriot_K(i, k, vecKik[ijkp], dposik, rcut_all, rcs_all, basepot_all, basepot_vecidx_all,
				 hv0_perriot_scr_K_frc_ik[ijkp],
				 hv1_perriot_scr_K_frc_ik[ijkp], hv2_perriot_scr_K_frc_ik[ijkp], hv3_perriot_scr_K_frc_ik[ijkp],
				 hv4_perriot_scr_K_frc_ik[ijkp], hv5_perriot_scr_K_frc_ik[ijkp],
				 hv6_perriot_scr_K_frc_ik[ijkp], hv7_perriot_scr_K_frc_ik[ijkp]);

	  // ###################### Perriot screening factor K ######################

	  
	}
	/* ################################################################### */	



	c2 = cik*cik;
	d2 = dik*dik;
	cost = (dposij * dposik) / (rij*rik);
	hcost  = hik + cost;
	hcost2 = hcost * hcost;

	gijk = gammaik * (1.0 + c2/d2 - c2 / (d2 + hcost2) );

	// *****************************************************************
	alphaijk = se_alphaijk;
	twomuik  = se_twomuik;
	F1 = 1.0;
	F2 = se_omegaijk;

	// alpha, omega
	if (p_potinfo->use_abop_alpha.elem(typei, typej, typek)){
	  if (! sys_single_elem)
	    alphaijk = p_potinfo->abop_alpha.elem(typei, typej, typek);
	  F1 = exp( alphaijk * (rij - rik) );
	  // .........................................................
	  if (p_potinfo->use_abop_omega.elem(typei, typej, typek)){
	    if (! sys_single_elem)
	      F2 = abop_omega.elem(typei, typej, typek);
	  }
	  else {
	    td = rij - r0ij - (rik - r0ik);
	    F2 = exp( alphaijk * td );
	  }
	}
	// 2mu
	else if (p_potinfo->use_abop_2mu.elem(typei, typek)){
	  if (! sys_single_elem)
	    twomuik = p_potinfo->abop_2mu.elem(typei, typek);
	  F1 = exp( twomuik * (rij - rik) );
	  F2 = 1.0;
	}
	// *****************************************************************
	
	
	Chiij += fcik * vecKik[ijkp] * gijk * F1 * F2;

	ijkp++;
      } // end of loop over neighbors k



      /*

      if (p_potinfo->use_abop_gamma.elem(typei, typej, typek))
	gamma_ijk = p_potinfo->abop_gamma.elem(typei, typej, typek);
      if (p_potinfo->use_3body_sym && p_potinfo->use_abop_gamma.elem(typek, typej, typei))
	gamma_ijk = p_potinfo->abop_gamma.elem(typek, typej, typei);

      */






      // bij = 1.0/sqrt(1.0 + Chiij);
      pij = abop_params_all.elem(typei,typej).p;
      bij = 1.0 / pow(1.0 + Chiij, pij);
      //std::cout << "bij " << bij << std::endl;


      // ################################################################
      // Energy:
      // ################################################################
      Epij = 0.0;
      if (pair_ij_tersoff || pair_ij_perriot_cut) // cutoff
	Epij = 0.5 * fcij * (VRij - bij * VAij);
      else if (pair_ij_perriot_scr) // screening
	Epij = 0.5 * Kij * (VRij - bij * VAij);
      
      Ep[i] += Epij;
      Ep_tot_local += Epij;




      if (get_pot_force == false) continue;






      // ################################################################
      // Forces which are not derivatives of bij, Kij:
      // ################################################################

      /* ############################ cutoff ############################ */
      for (p=0; p<3; ++p) frc_ij[p]=0.0;
      if (pair_ij_tersoff || pair_ij_perriot_cut){
	for (p=0; p<3; ++p){
	  frc_ij[p] = - 0.5 * ( dfcij * VRij + fcij * dVRij
				- bij * dfcij * VAij - bij * fcij * dVAij) * dposij[p]/rij;
	}
      }
      else if (pair_ij_perriot_scr){
	for (p=0; p<3; ++p){
	  frc_ij[p] = - 0.5 * Kij * ( dVRij - bij * dVAij) * dposij[p]/rij;
	}
      }
          
      for (int v1=0; v1<3; ++v1){
	for (int v2=0; v2<3; ++v2){
	  virials[i].elem(v1,v2) += 0.5 *   frc_ij[v1]  *   dposij[v2];
	  virials[j].elem(v1,v2) += 0.5 * (-frc_ij[v1]) * (-dposij[v2]);
	  ////virials[j].elem(v1,v2) += frcj[v1] * pos[j][v2];
	}
      }

      for (p=0; p<3; ++p){
	frci[p] =  frc_ij[p];
	frcj[p] = -frc_ij[p];
      }

      for (p=0; p<3; ++p){
	frc[i][p] += frci[p];
	frc[j][p] += frcj[p];
      }







      // ################################################################
      // Derivatives coming from the screening factor Kij:
      // ################################################################

      // ###################### Perriot screening factor K force ######################
      if (pair_ij_perriot_scr){
	double pref = - 0.5 * (VRij - bij * VAij);
	force_ABOP_perriot_K_frc(i, j, Kij, dposij, pref, rcut_all, rcs_all, basepot_all, basepot_vecidx_all,
				 hv0_perriot_scr_K_frc_ij,
				 hv1_perriot_scr_K_frc_ij, hv2_perriot_scr_K_frc_ij, hv3_perriot_scr_K_frc_ij,
				 hv4_perriot_scr_K_frc_ij, hv5_perriot_scr_K_frc_ij,
				 hv6_perriot_scr_K_frc_ij, hv7_perriot_scr_K_frc_ij);
      }
      // ###################### Perriot screening factor K force ######################



      // ################################################################
      // Forces which are derivatives of bij:
      // ################################################################

      /* ############################ cutoff ############################ */
      threebodyfactor=0.0;
      if (pair_ij_tersoff || pair_ij_perriot_cut){
	// threebodyfactor = - 0.25 * fcij * VAij * bij*bij*bij;
	threebodyfactor = - 0.5 * pij * fcij * VAij * pow(1.0 + Chiij, -pij-1.0);
      }
      else if (pair_ij_perriot_scr){
	threebodyfactor = - 0.5 * pij * Kij * VAij * pow(1.0 + Chiij, -pij-1.0);
      }



      // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      // Loop over k
      // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      ijkp = 0;
      ijkh = -1;
      for (ik=0; ik<neighborcollection[i].size(); ik++){
	ijkh++;

	/*
	  {Neighbors k} + atom j = {Neighbors j} !!!
	*/
	if ( ! abop_is_neigh_ij[ijkh]) continue;

	k = neighborcollection[i][ik];
	if (k==j) continue;

	typek = itype[k]; //elem.name2idx( matter[k] );

	/* ############################ cutoff/screening ############################ */
	pair_ik_tersoff     = rcs_all.elem(typei,typek).tersoff;
	pair_ik_perriot_cut = rcs_all.elem(typei,typek).perriot_cut;
	pair_ik_perriot_scr = rcs_all.elem(typei,typek).perriot_scr;
	Rik = rcs_all.elem(typei,typek).R;
	Dik = rcs_all.elem(typei,typek).D;
	rminik = rcs_all.elem(typei,typek).prmin;
	rmaxik = rcs_all.elem(typei,typek).prmax;
	rcutik = rcs_all.elem(typei,typek).rcut;
	/* ################################################################### */
	dposik = abop_dpos_ij[ijkh];
	rik = dposik.magn();
	/* ################################################################### */


	/* ############################ cutoff/screening ############################ */
	fcik = abop_fc_ij[ijkh];
	dfcik = abop_dfc_ij[ijkh];

	/*
	fcik  = 1.0;
	dfcik = 0.0;
	if (pair_ik_tersoff){
	  if (rik < Rik-Dik){
	    fcik  = 1.0;  dfcik = 0.0;
	  }
	  else if (rik > rcutik){
	    fcik  = 0.0;  dfcik = 0.0;
	    continue;
	  }
	  else {
	    td = 0.5*PI*(rik-Rik)/Dik;

	    fcik  = 0.5 * (1.0 - sin( td ));
	    dfcik = -0.5 * 0.5*PI/Dik * cos( td );
	  }
	}
	else if (pair_ik_perriot_cut || pair_ik_perriot_scr){
	  if (rik < rminik){
	    fcik  = 1.0;  dfcik = 0.0;
	  }
	  else if (rik > rmaxik){
	    fcik  = 0.0;  dfcik = 0.0;
	    continue;
	  }
	  else {
	    double invdrc = 1.0/(rmaxik - rminik);
	    double a = (rik - rminik) * invdrc;
	    
	    double aa = a*a;
	    double aaa = a*aa;
	    
	    double td1 = 6*aa - 15*a + 10;
	    fcik = 1.0 - aaa * td1;
	    dfcik = ( -3*aa * td1 - aaa * (12*a - 15) ) * invdrc;
	  }
	}
	*/



	r0ik = abop_params_all.elem(typei,typek).r0;

	gammaik = abop_params_all.elem(typei,typek).gamma;
	cik = abop_params_all.elem(typei,typek).c;
	dik = abop_params_all.elem(typei,typek).d;
	hik = abop_params_all.elem(typei,typek).h;
	
	c2 = cik*cik;
	d2 = dik*dik;
	cost = (dposij * dposik) / (rij*rik);
	hcost  = hik + cost;
	hcost2 = hcost * hcost;

	gijk = gammaik * (1.0 + c2/d2 - c2 / (d2 + hcost2) );



	for (p=0; p<3; ++p){
	  /*
	  dcost_i[p] = - cost/(rij*rij) *   dposij[p]
	    - cost/(rik*rik) * dposik[p]
	    + 1.0/(rij*rik) * (dposij[p] + dposik[p]);
	  
	  dcost_j[p] = - cost/(rij*rij) * (-dposij[p])
	    + 0.0
	    + 1.0/(rij*rik) * (-dposik[p]);

	  dcost_k[p] = - cost/(rik*rik) * (-dposik[p])
	    + 0.0
	    + 1.0/(rij*rik) * (-dposij[p]);
	  */

	  dcost_ij[p] = - cost/(rij*rij) * dposij[p]
	    + 1.0/(rij*rik) * dposik[p];
	  dcost_ik[p] = - cost/(rik*rik) * dposik[p]
	    + 1.0/(rij*rik) * dposij[p];
	}

	dgijk = gammaik * c2 /((d2 + hcost2)*(d2 + hcost2)) * 2 * hcost;
	for (p=0; p<3; ++p){
	  /*
	  dgijk_i[p] = dgijk * dcost_i[p];
	  dgijk_j[p] = dgijk * dcost_j[p];
	  dgijk_k[p] = dgijk * dcost_k[p];
	  */
	  dgijk_ij[p] = dgijk * dcost_ij[p];
	  dgijk_ik[p] = dgijk * dcost_ik[p];
	}






	dF1  = 0.0;
	for (p=0; p<3; ++p){
	  /*
	  dF1_i[p] = 0.0;
	  dF1_j[p] = 0.0;
	  dF1_k[p] = 0.0;
	  */
	  
	  dF1_ij[p]=0.0;
	  dF1_ik[p]=0.0;
	}
	dF2 = 0.0;
	for (p=0; p<3; ++p){
	  /*
	  dF2_i[p] = 0.0;
	  dF2_j[p] = 0.0;
	  dF2_k[p] = 0.0;
	  */

	  dF2_ij[p]=0.0;
	  dF2_ik[p]=0.0;
	}


	// *****************************************************************
	alphaijk = se_alphaijk;
	twomuik  = se_twomuik;
	F1  = 1.0;
	dF1 = 0.0;
	F2  = se_omegaijk;
	dF2 = 0.0;

	// alpha, omega
	if (p_potinfo->use_abop_alpha.elem(typei, typej, typek)){
	  if (! sys_single_elem)
	    alphaijk = p_potinfo->abop_alpha.elem(typei, typej, typek);
	  F1 = exp( alphaijk * (rij - rik) );
	  dF1 = alphaijk * F1;
	  for (p=0; p<3; ++p){
	    /*
	    dF1_i[p] = dF1 * ( dposij[p]/rij - dposik[p]/rik );
	    dF1_j[p] = dF1 * (-dposij[p]/rij);
	    dF1_k[p] = dF1 * ( dposik[p]/rik);
	    */
	    dF1_ij[p] = dF1 * dposij[p]/rij;
	    dF1_ik[p] = dF1 * -dposik[p]/rik;
	  }
	  // .........................................................
	  if (p_potinfo->use_abop_omega.elem(typei, typej, typek)){
	    if (! sys_single_elem){
	      F2  = abop_omega.elem(typei, typej, typek);
	      dF2 = 0.0;
	    }
	  }
	  else {
	    td = rij - r0ij - (rik - r0ik);
	    F2  = exp( alphaijk * td );
	    dF2 = alphaijk * F2;
	    for (p=0; p<3; ++p){
	      /*
	      dF2_i[p] = dF2 * ( dposij[p]/rij - dposik[p]/rik );
	      dF2_j[p] = dF2 * (-dposij[p]/rij);
	      dF2_k[p] = dF2 * ( dposik[p]/rik);
	      */
	      dF2_ij[p] = dF2 * dposij[p]/rij;
	      dF2_ik[p] = dF2 * (- dposik[p]/rik );
	    }
	  }
	}
	// 2mu
	else if (p_potinfo->use_abop_2mu.elem(typei, typek)){
	  if (! sys_single_elem)
	    twomuik = p_potinfo->abop_2mu.elem(typei, typek);
	  F1  = exp( twomuik * (rij - rik) );
	  dF1 = twomuik * F1;
	  for (p=0; p<3; ++p){
	    /*
	    dF1_i[p] = dF1 * ( dposij[p]/rij - dposik[p]/rik );
	    dF1_j[p] = dF1 * (-dposij[p]/rij);
	    dF1_k[p] = dF1 * ( dposik[p]/rik);
	    */
	    dF1_ij[p] = dF1 * dposij[p]/rij;
	    dF1_ik[p] = dF1 * (- dposik[p]/rik );
	  }
	  F2  = 1.0;
	  dF2 = 0.0;
	}
	// *****************************************************************


	//std::cout << "alpha omega  " << p_potinfo->abop_alpha.elem(typei, typej, typek) << " " << omegaijk << std::endl;







	// ################################################################
	// Derivatives coming from the screening factor Kik:
	// ################################################################

	// ###################### Perriot screening factor K force ######################
	if (pair_ik_perriot_scr){
	  double pref = threebodyfactor * fcik * gijk * F1 * F2 ;
	  force_ABOP_perriot_K_frc(i, k, vecKik[ijkp], dposik, pref, rcut_all, rcs_all, basepot_all, basepot_vecidx_all,
				   hv0_perriot_scr_K_frc_ik[ijkp],
				   hv1_perriot_scr_K_frc_ik[ijkp], hv2_perriot_scr_K_frc_ik[ijkp], hv3_perriot_scr_K_frc_ik[ijkp],
				   hv4_perriot_scr_K_frc_ik[ijkp], hv5_perriot_scr_K_frc_ik[ijkp],
				   hv6_perriot_scr_K_frc_ik[ijkp], hv7_perriot_scr_K_frc_ik[ijkp]);
	}
	else
	  vecKik[ijkp] = 1.0;

	Kik = vecKik[ijkp];
	// ###################### Perriot screening factor K force ######################



	// Common:
	for (p=0; p<3; ++p){
	  /*
	    frci[p] = 0.0
	    + threebodyfactor * dfcik * dposik[p]/rik * gijk * F1 * F2
	    + threebodyfactor * fcik * dgijk_i[p] * F1 * F2
	    + threebodyfactor * fcik * gijk * dF1_i[p] * F2
	    + threebodyfactor * fcik * gijk * F1 * dF2_i[p];
	  
	    frcj[p] = 0.0
	    + 0.0
	    + threebodyfactor * fcik * dgijk_j[p] * F1 * F2
	    + threebodyfactor * fcik * gijk * dF1_j[p] * F2
	    + threebodyfactor * fcik * gijk * F1 * dF2_j[p];
	  
	    frck[p] = 0.0
	    + threebodyfactor * dfcik * (-dposik[p]/rik) * gijk * F1 * F2
	    + threebodyfactor * fcik * dgijk_k[p] * F1 * F2
	    + threebodyfactor * fcik * gijk * dF1_k[p] * F2
	    + threebodyfactor * fcik * gijk * F1 * dF2_k[p];
	  */


	  frc_ij[p] = 0.0
	    + threebodyfactor * fcik * Kik * dgijk_ij[p] * F1 * F2
	    + threebodyfactor * fcik * Kik * gijk * dF1_ij[p] * F2
	    + threebodyfactor * fcik * Kik * gijk * F1 * dF2_ij[p];
	  frc_ik[p] = 0.0
	    + threebodyfactor * dfcik * dposik[p]/rik * Kik * gijk * F1 * F2
	    + threebodyfactor * fcik * Kik * dgijk_ik[p] * F1 * F2
	    + threebodyfactor * fcik * Kik * gijk * dF1_ik[p] * F2
	    + threebodyfactor * fcik * Kik * gijk * F1 * dF2_ik[p];
	}





	
	  
	for (int v1=0; v1<3; ++v1){
	  for (int v2=0; v2<3; ++v2){
	    /*
	    // ij contribution:
	    virials[i].elem(v1,v2) += 0.5*
	    ( 
	    threebodyfactor * fcik * dgijk * ( dposik[v1]/(rij*rik) - cost*dposij[v1]/(rij*rij))
	    * F1 * F2
	    +
	    threebodyfactor * fcik * gijk * dF1 * dposij[v1]/rij * F2
	    +
	    threebodyfactor * fcik * gijk * F1 * dF2 * dposij[v1]/rij
	    )
	    * dposij[v2];
	    virials[j].elem(v1,v2) += - 0.5*
	    ( 
	    threebodyfactor * fcik * dgijk * ( dposik[v1]/(rij*rik) - cost*dposij[v1]/(rij*rij))
	    * F1 * F2
	    +
	    threebodyfactor * fcik * gijk * dF1 * dposij[v1]/rij * F2
	    +
	    threebodyfactor * fcik * gijk * F1 * dF2 * dposij[v1]/rij
	    )
	    * (-dposij[v2]);

	    // ik contribution:
	    virials[i].elem(v1,v2) += 0.5*
	    ( 
	    threebodyfactor * dfcik * dposik[v1]/rik * gijk * F1 * F2
	    +
	    threebodyfactor * fcik * dgijk * ( dposij[v1]/(rij*rik) - cost*dposik[v1]/(rik*rik))
	    * F1 * F2
	    +
	    threebodyfactor * fcik * gijk * (-dF1 * dposik[v1]/rik) * F2
	    +
	    threebodyfactor * fcik * gijk * F1 * (-dF2 * dposik[v1]/rik)
	    )
	    * dposik[v2];
	    virials[k].elem(v1,v2) += - 0.5*
	    ( 
	    threebodyfactor * dfcik * dposik[v1]/rik * gijk * F1 * F2
	    +
	    threebodyfactor * fcik * dgijk * ( dposij[v1]/(rij*rik) - cost*dposik[v1]/(rik*rik))
	    * F1 * F2
	    +
	    threebodyfactor * fcik * gijk * (-dF1 * dposik[v1]/rik) * F2
	    +
	    threebodyfactor * fcik * gijk * F1 * (-dF2 * dposik[v1]/rik)
	    )
	    * (-dposik[v2]);
	    */

	    virials[i].elem(v1,v2) += 0.5 *   frc_ij[v1]  *   dposij[v2];
	    virials[j].elem(v1,v2) += 0.5 * (-frc_ij[v1]) * (-dposij[v2]);
	    virials[i].elem(v1,v2) += 0.5 *   frc_ik[v1]  *   dposik[v2];
	    virials[k].elem(v1,v2) += 0.5 * (-frc_ik[v1]) * (-dposik[v2]);
	    
	  }
	}

	

	for (p=0; p<3; ++p){
	  /*
	  frc[i][p] += frci[p];
	  frc[j][p] += frcj[p];
	  frc[k][p] += frck[p];
	  */

	  frc[i][p] +=  frc_ij[p] + frc_ik[p];
	  frc[j][p] += -frc_ij[p];
	  frc[k][p] += -frc_ik[p];
	}


	/*
	for (p=0; p<3; ++p){
	  frc[i][p] -= threebodyfactor * dfcik * gijk * expij * expik * omegaijk * dposik[p];
	  frc[k][p] -= threebodyfactor * dfcik * gijk * expij * expik * omegaijk * (-dposik[p]);

	  frc[i][p] -= threebodyfactor * fcik * dgijk_i[p] * expij * expik * omegaijk;
	  frc[j][p] -= threebodyfactor * fcik * dgijk_j[p] * expij * expik * omegaijk;
	  frc[k][p] -= threebodyfactor * fcik * dgijk_k[p] * expij * expik * omegaijk;

	  frc[i][p] -= threebodyfactor * fcik * gijk * dexpij_i[p] * expik * omegaijk;
	  frc[j][p] -= threebodyfactor * fcik * gijk * dexpij_j[p] * expik * omegaijk;

	  frc[i][p] -= threebodyfactor * fcik * gijk * expij * dexpik_i[p] * omegaijk;
	  frc[k][p] -= threebodyfactor * fcik * gijk * expij * dexpik_k[p] * omegaijk;
	}
	*/

	ijkp++;
      } // end of loop over neighbors k

      ijp++;
    } // end of loop over neighbors j
    //std::cout << "Atom i energy " << Ep[i] << std::endl;
  } // end of loop over atoms
  



  return Ep_tot_local;
}





void MDSystem::force_ABOP_perriot_K(int i,
				    int j,
				    double & Kij,
				    Vector3<double> & dposij,
				    Matrix<double> & rcut_all,
				    Matrix<CutoffScreeningPair> & rcs_all,
				    Matrix<std::string> & basepot_all,
				    Matrix<int> & basepot_vecidx_all,
				    Vector<bool>   & hv0_perriot_scr_K_frc,
				    Vector<double> & hv1_perriot_scr_K_frc,
				    Vector<double> & hv2_perriot_scr_K_frc,
				    Vector<double> & hv3_perriot_scr_K_frc,
				    Vector<double> & hv4_perriot_scr_K_frc,
				    Vector<double> & hv5_perriot_scr_K_frc,
				    Vector< Vector3<double> > & hv6_perriot_scr_K_frc,
				    Vector< Vector3<double> > & hv7_perriot_scr_K_frc
				    ){
  
  int k, typei, typej, typek;
  int se_ivecij, ivecij, ivecjk, ivecik;
  bool pair_ik_perriot_scr, pair_jk_perriot_scr;

  typei = itype[i]; //elem.name2idx( matter[i] );
  typej = itype[j]; //elem.name2idx( matter[j] );
  se_ivecij = basepot_vecidx_all.elem(typei, typei);

  Vector3<double> dposik, dposjk;
  double rik, rjk, rcutik, rcutjk, mik, mjk, Xik, Xjk, Tijk_n_sum, Tijk, nij;
  double rij, td, td1, td2;
  double dXik, dXjk;

  Kij = 1.0;

  nij = rcs_all.elem(typei,typej).pn;
  rij = dposij.magn();

  int ijkp=0;


  if (hv0_perriot_scr_K_frc.size() <= neighborcollection[i].size())
    hv0_perriot_scr_K_frc.resize( neighborcollection[i].size() );

  if (hv1_perriot_scr_K_frc.size() <= neighborcollection[i].size())
    hv1_perriot_scr_K_frc.resize( neighborcollection[i].size() );
  if (hv2_perriot_scr_K_frc.size() <= neighborcollection[i].size())
    hv2_perriot_scr_K_frc.resize( neighborcollection[i].size() );
  if (hv3_perriot_scr_K_frc.size() <= neighborcollection[i].size())
    hv3_perriot_scr_K_frc.resize( neighborcollection[i].size() );
  if (hv4_perriot_scr_K_frc.size() <= neighborcollection[i].size())
    hv4_perriot_scr_K_frc.resize( neighborcollection[i].size() );
  if (hv5_perriot_scr_K_frc.size() <= neighborcollection[i].size())
    hv5_perriot_scr_K_frc.resize( neighborcollection[i].size() );

  if (hv6_perriot_scr_K_frc.size() <= neighborcollection[i].size())
    hv6_perriot_scr_K_frc.resize( neighborcollection[i].size() );
  if (hv7_perriot_scr_K_frc.size() <= neighborcollection[i].size())
    hv7_perriot_scr_K_frc.resize( neighborcollection[i].size() );



  // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
  // Loop over k
  // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
  Tijk_n_sum=0.0;
  for (int ik=0; ik<neighborcollection[i].size(); ik++){
    k = neighborcollection[i][ik];

    hv0_perriot_scr_K_frc[ik] = false;

    if (k==j) continue;

    typek = itype[k]; //elem.name2idx( matter[k] );
    if (! iac_pure_ABOP){
      if (basepot_all.elem(typei, typek) != "ABOP" || basepot_all.elem(typej, typek) != "ABOP")
	continue;
    }

    ivecik = se_ivecij;
    ivecjk = se_ivecij;
    if (! sys_single_elem){
      ivecik = basepot_vecidx_all.elem(typei, typek);
      ivecjk = basepot_vecidx_all.elem(typej, typek);
    }
    if (ivecik<0 || ivecjk<0) continue;
	

    /* ############################ cutoff/screening ############################ */
    pair_ik_perriot_scr = rcs_all.elem(typei,typek).perriot_scr;
    pair_jk_perriot_scr = rcs_all.elem(typej,typek).perriot_scr;
    mik    = rcs_all.elem(typei,typek).pm;
    rcutik = rcs_all.elem(typei,typek).rcut;
    mjk    = rcs_all.elem(typej,typek).pm;
    rcutjk = rcs_all.elem(typej,typek).rcut;
    if (! pair_ik_perriot_scr) continue;
    if (! pair_jk_perriot_scr) continue;
    /* ################################################################### */
    get_atom_distance_vec(pos[i], pos[k], dposik);
    rik = dposik.magn();
    if (rik > rcutik) continue;
    get_atom_distance_vec(pos[j], pos[k], dposjk);
    rjk = dposjk.magn();
    if (rjk > rcutjk) continue;
    /* ################################################################### */
    if (fp_is_small(rik - rcutik) || fp_is_small(rjk - rcutjk))
      continue; // If true, then Xik or Xjk or both are very large.
    /* ################################################################### */


    Xik = rik/(1.0 - pow(rik/rcutik, mik));
    Xjk = rjk/(1.0 - pow(rjk/rcutjk, mjk));

    // *****************************************************
    if (Xik + Xjk > 3.0*rij) continue;
    // *****************************************************

    hv0_perriot_scr_K_frc[ik] = true;


    hv6_perriot_scr_K_frc[ijkp] = dposik;
    hv7_perriot_scr_K_frc[ijkp] = dposjk;

    Tijk = -0.5 + rij/(Xik + Xjk - rij);
    Tijk_n_sum += pow(Tijk, nij);

    td1 = Xik/rik;
    td2 = Xjk/rjk;
    dXik = td1 + mik * td1*td1 * pow(rik/rcutik, mik);
    dXjk = td2 + mjk * td2*td2 * pow(rjk/rcutjk, mjk);

    td1 = 1.0/(Xik + Xjk - rij);
    td2 = - rij * td1 * td1;

    hv1_perriot_scr_K_frc[ijkp] = td1;
    hv2_perriot_scr_K_frc[ijkp] = td2;

    hv4_perriot_scr_K_frc[ijkp] = dXik;
    hv5_perriot_scr_K_frc[ijkp] = dXjk;

    /*
    for (int p=0; p<3; ++p){
      dTijk_ij[p] = td1 *   dposij[p]/rij  + td2 * ( (-1)*  dposij[p]/rij  );
      dTijk_ik[p] = td2 * ( dXik *   dposik[p]/rik );

      dTijk_ji[p] = - dTijk_ij[p];
      dTijk_jk[p] = td2 * ( dXjk *   dposjk[p]/rjk );

      dTijk_ki[p] = - dTijk_ik[p];
      dTijk_kj[p] = - dTijk_jk[p];
    }
    */
    //td = Kij * nij * pow(Tijk, nij-1.0);
    td = nij * pow(Tijk, nij-1.0);
    hv3_perriot_scr_K_frc[ijkp] = td;

    ijkp++;
  }
  Kij = exp( - Tijk_n_sum);

}





/* Calculate and store forces on atoms i,j,k and the virials. The variable
   'pref' contains a constant prefactor due to other contributions to
   the forces.
*/
void MDSystem::force_ABOP_perriot_K_frc(int i,
					int j,
					double & Kij,
					Vector3<double> & dposij,
					double & pref,
					Matrix<double> & rcut_all,
					Matrix<CutoffScreeningPair> & rcs_all,
					Matrix<std::string> & basepot_all,
					Matrix<int> & basepot_vecidx_all,
					Vector<bool>   & hv0_perriot_scr_K_frc,
					Vector<double> & hv1_perriot_scr_K_frc,
					Vector<double> & hv2_perriot_scr_K_frc,
					Vector<double> & hv3_perriot_scr_K_frc,
					Vector<double> & hv4_perriot_scr_K_frc,
					Vector<double> & hv5_perriot_scr_K_frc,
					Vector< Vector3<double> > & hv6_perriot_scr_K_frc,
					Vector< Vector3<double> > & hv7_perriot_scr_K_frc
					){

  // Get Kij:
  // force_ABOP_perriot_K(i, j, Kij, dposij, rcut_all, rcs_all, basepot_all, basepot_vecidx_all);

  int k, typei, typej, typek;
  int se_ivecij, ivecij, ivecjk, ivecik;
  bool pair_ik_perriot_scr, pair_jk_perriot_scr;
  Vector3<double> dposik, dposjk;
  double rik, rjk, rij;
  Vector3<double> dTijk_ij, dTijk_ik, dTijk_jk;
  double dXik, dXjk, td, td0, td1, td2, td3;
  int ijkp=0;

  typei = itype[i]; //elem.name2idx( matter[i] );
  typej = itype[j]; //elem.name2idx( matter[j] );
  se_ivecij = basepot_vecidx_all.elem(typei, typei);

  rij = dposij.magn();

  for (int ik=0; ik<neighborcollection[i].size(); ik++){

    if ( ! hv0_perriot_scr_K_frc[ik] ) continue;

    k = neighborcollection[i][ik];
    /*
    get_atom_distance_vec(pos[i], pos[k], dposik);
    rik = dposik.magn();
    get_atom_distance_vec(pos[j], pos[k], dposjk);
    rjk = dposjk.magn();
    */
    dposik = hv6_perriot_scr_K_frc[ijkp];
    dposjk = hv7_perriot_scr_K_frc[ijkp];
    rik = dposik.magn();
    rjk = dposjk.magn();
    

    td1 = hv1_perriot_scr_K_frc[ijkp];
    td2 = hv2_perriot_scr_K_frc[ijkp];
    dXik = hv4_perriot_scr_K_frc[ijkp];
    dXjk = hv5_perriot_scr_K_frc[ijkp];

    for (int p=0; p<3; ++p){
      dTijk_ij[p] = td1 *   dposij[p]/rij  + td2 * ( (-1)*  dposij[p]/rij  );
      dTijk_ik[p] = td2 * ( dXik *   dposik[p]/rik );

      // dTijk_ji[p] = - dTijk_ij[p];
      // //dTijk_ji[p] = td1 * (-dposij[p]/rij) + td2 * ( (-1)*(-dposij[p]/rij) );
      dTijk_jk[p] = td2 * ( dXjk *   dposjk[p]/rjk );

      // dTijk_ki[p] = - dTijk_ik[p];
      // dTijk_kj[p] = - dTijk_jk[p];
      //dTijk_ki[p] = td2 * ( dXik * (-dposik[p]/rik));
      //dTijk_kj[p] = td2 * ( dXjk * (-dposjk[p]/rjk));
    }

    td  = Kij * hv3_perriot_scr_K_frc[ijkp];


#if 0
    for (int p=0; p<3; ++p){
      dKij_ij[p] = - td * (dTijk_ij[p]);
      dKij_ik[p] = - td * (dTijk_ik[p]);
      
      dKij_ji[p] = - td * (dTijk_ji[p]);
      dKij_jk[p] = - td * (dTijk_jk[p]);
      
      dKij_ki[p] = - td * (dTijk_ki[p]);
      dKij_kj[p] = - td * (dTijk_kj[p]);
    }

    // Kij: prefactor = - 0.5 * (VRij - bij * VAij) ;
    // Kik: prefactor = threebodyfactor * fcik * gijk * F1 * F2 ;

    // Use previously calc. results to get total pairwise forces:
    for (int p=0; p<3; ++p){
      frc_ij[p] = dKij_ij[p] * pref;
      frc_ik[p] = dKij_ik[p] * pref;
      frc_jk[p] = dKij_jk[p] * pref;
    }
#endif


    // Add to total atomic forces for all atoms participating:
    for (int p=0; p<3; ++p){
      frc[i][p] +=   - ( dTijk_ij[p] + dTijk_ik[p]) * td * pref;
      frc[j][p] +=   - (-dTijk_ij[p] + dTijk_jk[p]) * td * pref;
      frc[k][p] +=   - (-dTijk_ik[p] - dTijk_jk[p]) * td * pref;
      /*
      frc[i][p] +=   frc_ij[p] + frc_ik[p];
      frc[j][p] += - frc_ij[p] + frc_jk[p];
      frc[k][p] +=  -frc_ik[p] - frc_jk[p];
      */
    }

    // Add to total atomic virials for all atoms participating:
    td0 = - 0.5 * td * pref;
    for (int v1=0; v1<3; ++v1){
      td1 = td0 * dTijk_ij[v1];
      td2 = td0 * dTijk_ik[v1];
      td3 = td0 * dTijk_jk[v1];

      for (int v2=0; v2<3; ++v2){
	// ij
	virials[i].elem(v1,v2) +=  td1 *   dposij[v2];
	virials[j].elem(v1,v2) += -td1 * (-dposij[v2]);
	// ik
	virials[i].elem(v1,v2) +=  td2 *   dposik[v2];
	virials[k].elem(v1,v2) += -td2 * (-dposik[v2]);
	// jk
	virials[j].elem(v1,v2) +=  td3 *   dposjk[v2];
	virials[k].elem(v1,v2) += -td3 * (-dposjk[v2]);


#if 0
	// ij
	virials[i].elem(v1,v2) += 0.5 *   frc_ij[v1]  *   dposij[v2];
	virials[j].elem(v1,v2) += 0.5 * (-frc_ij[v1]) * (-dposij[v2]);
	// ik
	virials[i].elem(v1,v2) += 0.5 *   frc_ik[v1]  *   dposik[v2];
	virials[k].elem(v1,v2) += 0.5 * (-frc_ik[v1]) * (-dposik[v2]);
	// jk
	virials[j].elem(v1,v2) += 0.5 *   frc_jk[v1]  *   dposjk[v2];
	virials[k].elem(v1,v2) += 0.5 * (-frc_jk[v1]) * (-dposjk[v2]);
#endif	
	      /*
	      // ij contributions to i
	      virials[i].elem(v1,v2) += 0.5*
		( - 0.5 * dKij_ij[v1] * (VRij - bij * VAij) ) * dposij[v2];
	      // ij contributions to j
	      virials[j].elem(v1,v2) += 0.5*
		( - 0.5 * dKij_ji[v1] * (VRij - bij * VAij) ) * (-dposij[v2]);

	      // ik contributions to i
	      virials[i].elem(v1,v2) += 0.5*
		( - 0.5 * dKij_ik[v1] * (VRij - bij * VAij) ) * dposik[v2];
	      // ik contributions to k
	      virials[k].elem(v1,v2) += 0.5*
		( - 0.5 * dKij_ki[v1] * (VRij - bij * VAij) ) * (-dposik[v2]);

	      // jk contributions to j
	      virials[j].elem(v1,v2) += 0.5*
		( - 0.5 * dKij_jk[v1] * (VRij - bij * VAij) ) * dposjk[v2];
	      // jk contributions to k
	      virials[k].elem(v1,v2) += 0.5*
		( - 0.5 * dKij_kj[v1] * (VRij - bij * VAij) ) * (-dposjk[v2]);
	      */
      }
    }
   
    ijkp++;
  }



}

