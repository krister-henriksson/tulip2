


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

  int se_typei = elem.name2idx( matter[0] );
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

  bool pair_ij_perriot;
  bool pair_ik_perriot;
  bool pair_jk_perriot;
  bool pair_is_perriot;
  bool pair_ks_perriot;

  double Kij, Tijk_n_sum, Tijk, Xik, Xjk, dXik, dXjk, mik, mjk, nij, rjk, rcutjk;
  double Kik, Tiks_n_sum, Tiks, Xis, Xks, dXis, dXks, mis, mks, nik, rks, rcutks;
  Vector3<double> dTijk_ij,dTijk_ji, dTijk_ik,dTijk_ki, dTijk_jk,dTijk_kj;
  Vector3<double> dTiks_ik,dTiks_ki, dTiks_is,dTiks_si, dTiks_ks,dTiks_sk;
  double rminik, rmaxik;

  Vector3<double> dKij_ij,dKij_ji, dKij_ik,dKij_ki, dKij_jk,dKij_kj;
  Vector3<double> dKik_ik,dKik_ki, dKik_is,dKik_si, dKik_ks,dKik_sk;














  // iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  // Loop over atoms i
  // iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

  //#pragma omp parallel for reduction(+:Ep_tot_local)
  for (i=0; i<nat; ++i){

    typei = elem.name2idx( matter[i] );

    // jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
    // Loop over j
    // jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj

    for (ij=0; ij<neighborcollection[i].size(); ij++){
      j = neighborcollection[i][ij];
      if (i==j) continue;

      typej = elem.name2idx( matter[j] );
     
      if (! iac_pure_ABOP) if (p_potinfo->basepot(typei,typej) != "ABOP") continue;

      ivecij = se_ivecij;
      if (! sys_single_elem) ivecij = p_potinfo->basepot_vecidx(typei, typej);
      if (ivecij<0) continue;
      
      get_atom_distance_vec(pos[i], pos[j], dposij);
      rij = dposij.magn();

      /* ############################ cutoff ############################ */
      pair_ij_tersoff = false;
      pair_ij_perriot = false;
      if (p_potinfo->pot_ABOP[ivecij].rcut_fun=="tersoff"){
	pair_ij_tersoff = true;
	Rij = p_potinfo->pot_ABOP[ivecij].parname2val("R");
	Dij = p_potinfo->pot_ABOP[ivecij].parname2val("D");
	rcutij = p_potinfo->pot_ABOP[ivecij].rcut();
      }
      else if (p_potinfo->pot_ABOP[ivecij].rcut_fun=="perriot"){
	pair_ij_perriot = true;
	nij = p_potinfo->pot_ABOP[ivecij].parname2val("pn");
	rcutij = p_potinfo->pot_ABOP[ivecij].parname2val("prcut");
      }
      if (rij > rcutij) continue;
      /* ################################################################### */


      D0ij   = p_potinfo->pot_ABOP[ivecij].parname2val("D0");
      r0ij   = p_potinfo->pot_ABOP[ivecij].parname2val("r0");
      betaij = p_potinfo->pot_ABOP[ivecij].parname2val("beta");
      Sij    = p_potinfo->pot_ABOP[ivecij].parname2val("S");
      /*
      std::cout << "D0ij " << D0ij << " r0ij " << r0ij << " betaij " << betaij << " Sij " << Sij
	   << " gammaij " << gammaij << " cij " << cij << " dij " << dij << " hij " << hij
	   << " Rij " << Rij << " Dij " << Dij << std::endl;
      */

      /* ############################ cutoff ############################ */
      if (pair_ij_tersoff){

	if (rij < Rij-Dij){
	  fcij = 1.0;  dfcij = 0.0;
	}
	else if (rij > Rij+Dij){
	  fcij = 0.0;  dfcij = 0.0;
	  continue;
	}
	else {
	  td = 0.5*PI*(rij-Rij)/Dij;

	  fcij  = 0.5 * (1.0 - sin( td ));
	  dfcij = -0.5 * 0.5*PI/Dij * cos( td );
	}

      }
      /* ################################################################### */

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
	   (! sys_single_elem && p_potinfo->use_reppot(typei, typej)) ){

	ivec_reppot = se_ivec_reppot;
	if (! sys_single_elem)
	  ivec_reppot = p_potinfo->reppot_vecidx(typei, typej);
	//ivec_reppot = p_potinfo->reppot_vecidx(s1,s2);
	if (ivec_reppot<0) break;

	Nr = p_potinfo->pot_Reppot[ivec_reppot].r_rep.size();
	td = p_potinfo->pot_Reppot[ivec_reppot].r_rep[Nr-1];

	bfermi = p_potinfo->pot_ABOP[ivecij].parname2val("bfermi");
	rfermi = p_potinfo->pot_ABOP[ivecij].parname2val("rfermi");
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


      if (pair_ij_perriot){
	force_ABOP_perriot_K(i, j, Kij, dposij );
#if 0
	// kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
	// Loop over k
	// kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
	Tijk_n_sum=0.0;
	for (ik=0; ik<neighborcollection[i].size(); ik++){
	  k = neighborcollection[i][ik];
	  if (k==i) continue;
	  if (k==j) continue;
	  typek = elem.name2idx( matter[k] );
	  if (! iac_pure_ABOP) if (p_potinfo->basepot(typei, typek) != "ABOP") continue;
	  if (! iac_pure_ABOP) if (p_potinfo->basepot(typej, typek) != "ABOP") continue;

	  ivecik = se_ivecij;
	  if (! sys_single_elem) ivecik = p_potinfo->basepot_vecidx(typei, typek);
	  if (ivecik<0) continue;

	  ivecjk = se_ivecij;
	  if (! sys_single_elem) ivecjk = p_potinfo->basepot_vecidx(typej, typek);
	  if (ivecjk<0) continue;
	
	  get_atom_distance_vec(pos[i], pos[k], dposik);
	  rik = dposik.magn();
	  get_atom_distance_vec(pos[j], pos[k], dposjk);
	  rjk = dposjk.magn();

	  pair_ik_perriot = false;
	  if (p_potinfo->pot_ABOP[ivecik].rcut_fun=="perriot"){
	    pair_ik_perriot = true;
	    rcutik = p_potinfo->pot_ABOP[ivecik].parname2val("prcut");
	    mik = p_potinfo->pot_ABOP[ivecik].parname2val("pm");
	  }
	  pair_jk_perriot = false;
	  if (p_potinfo->pot_ABOP[ivecjk].rcut_fun=="perriot"){
	    pair_jk_perriot = true;
	    rcutjk = p_potinfo->pot_ABOP[ivecjk].parname2val("prcut");
	    mjk = p_potinfo->pot_ABOP[ivecjk].parname2val("pm");
	  }

	  if (pair_ik_perriot==false || pair_jk_perriot==false)
	    continue;

	  if (rik > rcutik) continue;
	  if (rjk > rcutjk) continue;

	  Xik = rik/(1.0 - pow(rik/rcutik, mik));
	  Xjk = rjk/(1.0 - pow(rjk/rcutjk, mjk));

	  // *****************************************************
	  if (Xik + Xjk > 3.0*rij) continue;
	  // *****************************************************

	  Tijk = 0.0;
	  if (Xik + Xjk < 3.0*rij){
	    Tijk = -0.5 + rij/(Xik + Xjk - rij);
	  }
	  td = 0.0;
	  if (Tijk>0.0) td = pow(Tijk, nij);
	  Tijk_n_sum += td;
	}
	Kij = exp( - Tijk_n_sum);
#endif
      }





      // ################################################################
      // Bond order factor bij:
      // ################################################################

      // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      // Loop over k
      // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

      Chiij = 0.0;
      for (ik=0; ik<neighborcollection[i].size(); ik++){
	k = neighborcollection[i][ik];

	if (k==i) continue;
	if (k==j) continue;
	//std::cout << "ijk " << i << j << k << std::endl;

	typek = elem.name2idx( matter[k] );
	//s3 = matter[k]; //p_potinfo->elem.idx2name(typek);

	if (! iac_pure_ABOP) if (p_potinfo->basepot(typei, typek) != "ABOP") continue;
	//if (p_potinfo->basepot(s1,s3) != "ABOP") continue;

	ivecik = se_ivecij;
	if (! sys_single_elem) ivecik = p_potinfo->basepot_vecidx(typei, typek);
	//ivecik = p_potinfo->basepot_vecidx(s1,s3);
	if (ivecik<0) continue;
	
	get_atom_distance_vec(pos[i], pos[k], dposik);
	rik = dposik.magn();

	r0ik = p_potinfo->pot_ABOP[ivecik].parname2val("r0");

	/* ############################ cutoff ############################ */
	pair_ik_tersoff = false;
	pair_ik_perriot = false;
	if (p_potinfo->pot_ABOP[ivecik].rcut_fun=="tersoff"){
	  pair_ik_tersoff = true;
	  Rik = p_potinfo->pot_ABOP[ivecik].parname2val("R");
	  Dik = p_potinfo->pot_ABOP[ivecik].parname2val("D");
	  rcutik = p_potinfo->pot_ABOP[ivecik].rcut();
	}
	if (p_potinfo->pot_ABOP[ivecik].rcut_fun=="perriot"){
	  pair_ik_perriot = true;
	  rcutik = p_potinfo->pot_ABOP[ivecik].parname2val("prcut");
	  rminik = p_potinfo->pot_ABOP[ivecik].parname2val("prmin");
	  rmaxik = p_potinfo->pot_ABOP[ivecik].parname2val("prmax");
	}
	if (rik > rcutik) continue;

	gammaik = p_potinfo->pot_ABOP[ivecik].parname2val("gamma");
	cik = p_potinfo->pot_ABOP[ivecik].parname2val("c");
	dik = p_potinfo->pot_ABOP[ivecik].parname2val("d");
	hik = p_potinfo->pot_ABOP[ivecik].parname2val("h");
	/*
	std::cout << "D0ik " << D0ik << " r0ik " << r0ik << " betaik " << betaik << " Sik " << Sik
	     << " gammaik " << gammaik << " cik " << cik << " dik " << dik << " hik " << hik
	     << " Rik " << Rik << " Dik " << Dik << std::endl;
	*/

	/* ############################ cutoff ############################ */
	fcik  = 1.0;
	dfcik = 0.0;
	Kik   = 1.0;
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
	else if (pair_ik_perriot){

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

	    double td1 = 6*aa - 15*a - 10;
	    fcik = 1.0 - aaa * td1;
	    dfcik = ( -3*aa * td1 - aaa * (12*a - 15) ) * invdrc;
	  }

	  force_ABOP_perriot_K(i, k, Kik, dposik);

#if 0
	  // Screening factor Kik:
	  // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
	  // Loop over s
	  // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
	  Tiks_n_sum=0.0;
	  for (is=0; is<neighborcollection[i].size(); is++){
	    s = neighborcollection[i][is];
	    if (s==i) continue;
	    if (s==k) continue;
	    types = elem.name2idx( matter[s] );
	    if (! iac_pure_ABOP) if (p_potinfo->basepot(typei, types) != "ABOP") continue;
	    if (! iac_pure_ABOP) if (p_potinfo->basepot(typek, types) != "ABOP") continue;

	    ivecis = se_ivecij;
	    if (! sys_single_elem) ivecis = p_potinfo->basepot_vecidx(typei, types);
	    if (ivecis<0) continue;

	    ivecks = se_ivecij;
	    if (! sys_single_elem) ivecks = p_potinfo->basepot_vecidx(typek, types);
	    if (ivecks<0) continue;
	    
	    get_atom_distance_vec(pos[i], pos[s], dposis);
	    ris = dposis.magn();
	    get_atom_distance_vec(pos[k], pos[s], dposks);
	    rks = dposks.magn();

	    pair_is_perriot = false;
	    if (p_potinfo->pot_ABOP[ivecis].rcut_fun=="perriot"){
	      pair_is_perriot = true;
	      rcutis = p_potinfo->pot_ABOP[ivecis].parname2val("prcut");
	      mis = p_potinfo->pot_ABOP[ivecis].parname2val("pm");
	    }
	    pair_ks_perriot = false;
	    if (p_potinfo->pot_ABOP[ivecks].rcut_fun=="perriot"){
	      pair_ks_perriot = true;
	      rcutks = p_potinfo->pot_ABOP[ivecks].parname2val("prcut");
	      mks = p_potinfo->pot_ABOP[ivecks].parname2val("pm");
	    }

	    if (pair_is_perriot==false || pair_ks_perriot==false)
	      continue;

	    if (ris > rcutis) continue;
	    if (rks > rcutks) continue;

	    Xis = ris/(1.0 - pow(ris/rcutis, mis));
	    Xks = rks/(1.0 - pow(rks/rcutks, mks));

	    // *****************************************************
	    if (Xis + Xks > 3.0*rik) continue;
	    // *****************************************************
	    
	    Tiks = 0.0;
	    if (Xis + Xks < 3.0*rik){
	      Tiks = -0.5 + rik/(Xis + Xks - rik);
	    }

	    td = 0.0;
	    if (Tiks>0.0) td = pow(Tiks, nik);
	    Tiks_n_sum += td;
	  }
	  Kik = exp( - Tiks_n_sum);
#endif
	}



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
	      F2 = p_potinfo->get_abop_omega(typei, typej, typek);
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


	Chiij += fcik * Kik * gijk * F1 * F2;



      } // end of loop over neighbors k









      // bij = 1.0/sqrt(1.0 + Chiij);
      pij = p_potinfo->pot_ABOP[ivecij].parname2val("p");
      bij = 1.0 / pow(1.0 + Chiij, pij);
      //std::cout << "bij " << bij << std::endl;


      // ################################################################
      // Energy:
      // ################################################################
      Epij = 0.0;
      if (pair_ij_tersoff)
	Epij = 0.5 * fcij * (VRij - bij * VAij);
      else if (pair_ij_perriot)
	Epij = 0.5 * Kij * (VRij - bij * VAij);

      Ep[i] += Epij;
      Ep_tot_local += Epij;




      if (get_pot_force == false) continue;






      // ################################################################
      // Forces which are not derivatives of bij, Kij:
      // ################################################################

      /* ############################ cutoff ############################ */

      if (pair_ij_tersoff){
	for (p=0; p<3; ++p){
	  frc_ij[p] = - 0.5 * ( dfcij * VRij + fcij * dVRij
				- bij * dfcij * VAij - bij * fcij * dVAij) * dposij[p]/rij;
	}
      }
      else if (pair_ij_perriot){
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
      if (pair_ij_perriot){
	double pref = - 0.5 * (VRij - bij * VAij);
	force_ABOP_perriot_K_frc(i, j, Kij, dposij, pref);

#if 0
	// kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
	// Loop over k
	// kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
	Tijk_n_sum=0.0;
	for (ik=0; ik<neighborcollection[i].size(); ik++){
	  k = neighborcollection[i][ik];
	  if (k==i) continue;
	  if (k==j) continue;
	  typek = elem.name2idx( matter[k] );
	  if (! iac_pure_ABOP) if (p_potinfo->basepot(typei, typek) != "ABOP") continue;
	  if (! iac_pure_ABOP) if (p_potinfo->basepot(typej, typek) != "ABOP") continue;

	  ivecik = se_ivecij;
	  if (! sys_single_elem) ivecik = p_potinfo->basepot_vecidx(typei, typek);
	  if (ivecik<0) continue;

	  ivecjk = se_ivecij;
	  if (! sys_single_elem) ivecjk = p_potinfo->basepot_vecidx(typej, typek);
	  if (ivecjk<0) continue;
	
	  get_atom_distance_vec(pos[i], pos[k], dposik);
	  rik = dposik.magn();
	  get_atom_distance_vec(pos[j], pos[k], dposjk);
	  rjk = dposjk.magn();

	  pair_ik_perriot = false;
	  if (p_potinfo->pot_ABOP[ivecik].rcut_fun=="perriot"){
	    pair_ik_perriot = true;
	    rcutik = p_potinfo->pot_ABOP[ivecik].parname2val("prcut");
	    mik = p_potinfo->pot_ABOP[ivecik].parname2val("pm");
	  }
	  pair_jk_perriot = false;
	  if (p_potinfo->pot_ABOP[ivecjk].rcut_fun=="perriot"){
	    pair_jk_perriot = true;
	    rcutjk = p_potinfo->pot_ABOP[ivecjk].parname2val("prcut");
	    mjk = p_potinfo->pot_ABOP[ivecjk].parname2val("pm");
	  }

	  if (pair_ik_perriot==false || pair_jk_perriot==false)
	    continue;

	  if (rik > rcutik) continue;
	  if (rjk > rcutjk) continue;

	  Xik = rik/(1.0 - pow(rik/rcutik, mik));
	  Xjk = rjk/(1.0 - pow(rjk/rcutjk, mjk));

	  // *****************************************************
	  if (Xik + Xjk > 3.0*rij) continue;
	  // *****************************************************

	  Tijk = 0.0;
	  if (Xik + Xjk < 3.0*rij){
	    Tijk = -0.5 + rij/(Xik + Xjk - rij);
	  }

	  dXik = Xik/rik + mik * (Xik/rik)*(Xik/rik) * pow(rik/rcutik, mik);
	  dXjk = Xjk/rjk + mjk * (Xjk/rjk)*(Xjk/rjk) * pow(rjk/rcutjk, mjk);

	  td1 = 1.0/(Xik + Xjk - rij);
	  td2 = - rij * td1 * td1;

	  for (p=0; p<3; ++p){
	    dTijk_ij[p] = td1 *   dposij[p]/rij  + td2 * ( (-1)*  dposij[p]/rij  );
	    dTijk_ik[p] = td2 * ( dXik *   dposik[p]/rik );

	    dTijk_ji[p] = - dTijk_ij[p];
	    //dTijk_ji[p] = td1 * (-dposij[p]/rij) + td2 * ( (-1)*(-dposij[p]/rij) );
	    dTijk_jk[p] = td2 * ( dXjk *   dposjk[p]/rjk );

	    dTijk_ki[p] = - dTijk_ik[p];
	    dTijk_kj[p] = - dTijk_jk[p];
	    //dTijk_ki[p] = td2 * ( dXik * (-dposik[p]/rik));
	    //dTijk_kj[p] = td2 * ( dXjk * (-dposjk[p]/rjk));
	  }

	  td = nij * pow(Tijk, nij-1.0);
	  for (p=0; p<3; ++p){
	    dKij_ij[p] = - Kij * td * (dTijk_ij[p]);
	    dKij_ik[p] = - Kij * td * (dTijk_ik[p]);

	    dKij_ji[p] = - Kij * td * (dTijk_ji[p]);
	    dKij_jk[p] = - Kij * td * (dTijk_jk[p]);

	    dKij_ki[p] = - Kij * td * (dTijk_ki[p]);
	    dKij_kj[p] = - Kij * td * (dTijk_kj[p]);
	  }





	  // Use previously calc. results to get total pairwise forces:
	  for (p=0; p<3; ++p){
	    frc_ij[p] = - 0.5 * dKij_ij[p] * (VRij - bij * VAij) ;
	    frc_ik[p] = - 0.5 * dKij_ik[p] * (VRij - bij * VAij) ;
	    frc_jk[p] = - 0.5 * dKij_jk[p] * (VRij - bij * VAij) ;
	  }





	  // Add to total atomic forces for all atoms participating:
	  for (p=0; p<3; ++p){
	    frc[i][p] +=   frc_ij[p] + frc_ik[p];
	    frc[j][p] += - frc_ij[p] + frc_jk[p];
	    frc[k][p] +=  -frc_ik[p] - frc_jk[p];
	    /*
	    frc[i][p] += - 0.5 * (dKij_ij[p] + dKij_ik[p]) * (VRij - bij * VAij) ;
	    frc[j][p] += - 0.5 * (dKij_ji[p] + dKij_jk[p]) * (VRij - bij * VAij) ;
	    frc[k][p] += - 0.5 * (dKij_ki[p] + dKij_kj[p]) * (VRij - bij * VAij) ;
	    */
	  }

	  // Add to total atomic virials for all atoms participating:
	  for (int v1=0; v1<3; ++v1){
	    for (int v2=0; v2<3; ++v2){
	      // ij
	      virials[i].elem(v1,v2) += 0.5 *   frc_ij[v1]  *   dposij[v2];
	      virials[j].elem(v1,v2) += 0.5 * (-frc_ij[v1]) * (-dposij[v2]);
	      // ik
	      virials[i].elem(v1,v2) += 0.5 *   frc_ik[v1]  *   dposik[v2];
	      virials[k].elem(v1,v2) += 0.5 * (-frc_ik[v1]) * (-dposik[v2]);
	      // jk
	      virials[j].elem(v1,v2) += 0.5 *   frc_jk[v1]  *   dposjk[v2];
	      virials[k].elem(v1,v2) += 0.5 * (-frc_jk[v1]) * (-dposjk[v2]);

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
	  

	}
#endif

      }
      




      // ################################################################
      // Forces which are derivatives of bij:
      // ################################################################

      /* ############################ cutoff ############################ */
      if (pair_ij_tersoff){
	// threebodyfactor = - 0.25 * fcij * VAij * bij*bij*bij;
	threebodyfactor = - 0.5 * pij * fcij * VAij * pow(1.0 + Chiij, -pij-1.0);
      }
      else if (pair_ij_perriot){
	threebodyfactor = - 0.5 * pij * Kij * VAij * pow(1.0 + Chiij, -pij-1.0);
      }



      // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      // Loop over k
      // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

      for (ik=0; ik<neighborcollection[i].size(); ik++){
	k = neighborcollection[i][ik];

	if (k==i) continue;
	if (k==j) continue;

	typek = elem.name2idx( matter[k] );
	//s3 = matter[k]; //p_potinfo->elem.idx2name(typek);

	if (! iac_pure_ABOP)
	  if (p_potinfo->basepot(typei, typek) != "ABOP") continue;
	//if (p_potinfo->basepot(s1,s3) != "ABOP") continue;

	ivecik = se_ivecij;
	if (! sys_single_elem)
	  ivecik = p_potinfo->basepot_vecidx(typei, typek);
	//ivecik = p_potinfo->basepot_vecidx(s1,s3);
	if (ivecik<0) continue;
	
	get_atom_distance_vec(pos[i], pos[k], dposik);
	rik = dposik.magn();

	r0ik = p_potinfo->pot_ABOP[ivecik].parname2val("r0");

	gammaik = p_potinfo->pot_ABOP[ivecik].parname2val("gamma");
	cik = p_potinfo->pot_ABOP[ivecik].parname2val("c");
	dik = p_potinfo->pot_ABOP[ivecik].parname2val("d");
	hik = p_potinfo->pot_ABOP[ivecik].parname2val("h");


	/* ############################ cutoff ############################ */
	pair_ik_tersoff = false;
	pair_ik_perriot = false;
	if (p_potinfo->pot_ABOP[ivecik].rcut_fun=="tersoff"){
	  pair_ik_tersoff = true;
	  Rik = p_potinfo->pot_ABOP[ivecik].parname2val("R");
	  Dik = p_potinfo->pot_ABOP[ivecik].parname2val("D");
	  rcutik = p_potinfo->pot_ABOP[ivecik].rcut();
	}
	if (p_potinfo->pot_ABOP[ivecik].rcut_fun=="perriot"){
	  pair_ik_perriot = true;
	  rcutik = p_potinfo->pot_ABOP[ivecik].parname2val("prcut");
	  rminik  = p_potinfo->pot_ABOP[ivecik].parname2val("prmin");
	  rmaxik  = p_potinfo->pot_ABOP[ivecik].parname2val("prmax");
	}
	if (rik > rcutik) continue;

	/* ############################ cutoff ############################ */
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
	else if (pair_ik_perriot){

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

	    double td1 = 6*aa - 15*a - 10;
	    fcik = 1.0 - aaa * td1;
	    dfcik = ( -3*aa * td1 - aaa * (12*a - 15) ) * invdrc;
	  }
	}



	Kik = 1.0;

#if 0
	if (pair_ik_perriot){

	  // Screening factor Kik:
	  Tiks_n_sum=0.0;
	  for (is=0; is<neighborcollection[i].size(); is++){
	    s = neighborcollection[i][is];
	    if (s==i) continue;
	    if (s==k) continue;
	    types = elem.name2idx( matter[s] );
	    if (! iac_pure_ABOP) if (p_potinfo->basepot(typei, types) != "ABOP") continue;
	    if (! iac_pure_ABOP) if (p_potinfo->basepot(typek, types) != "ABOP") continue;

	    ivecis = se_ivecij;
	    if (! sys_single_elem) ivecis = p_potinfo->basepot_vecidx(typei, types);
	    if (ivecis<0) continue;

	    ivecks = se_ivecij;
	    if (! sys_single_elem) ivecks = p_potinfo->basepot_vecidx(typek, types);
	    if (ivecks<0) continue;
	    
	    get_atom_distance_vec(pos[i], pos[s], dposis);
	    ris = dposis.magn();
	    get_atom_distance_vec(pos[k], pos[s], dposks);
	    rks = dposks.magn();

	    pair_is_perriot = false;
	    if (p_potinfo->pot_ABOP[ivecis].rcut_fun=="perriot"){
	      pair_is_perriot = true;
	      rcutis = p_potinfo->pot_ABOP[ivecis].parname2val("prcut");
	      mis = p_potinfo->pot_ABOP[ivecis].parname2val("pm");
	    }
	    pair_ks_perriot = false;
	    if (p_potinfo->pot_ABOP[ivecks].rcut_fun=="perriot"){
	      pair_ks_perriot = true;
	      rcutks = p_potinfo->pot_ABOP[ivecks].parname2val("prcut");
	      mks = p_potinfo->pot_ABOP[ivecks].parname2val("pm");
	    }

	    if (pair_is_perriot==false || pair_ks_perriot==false)
	      continue;

	    if (ris > rcutis) continue;
	    if (rks > rcutks) continue;

	    Xis = ris/(1.0 - pow(ris/rcutis, mis));
	    Xks = rks/(1.0 - pow(rks/rcutks, mks));

	    // *****************************************************
	    if (Xis + Xks > 3.0*rik) continue;
	    // *****************************************************
	    
	    Tiks = 0.0;
	    if (Xis + Xks < 3.0*rik){
	      Tiks = -0.5 + rik/(Xis + Xks - rik);
	    }

	    td = 0.0;
	    if (Tiks>0.0) td = pow(Tiks, nik);
	    Tiks_n_sum += td;
	  }
	  Kik = exp( - Tiks_n_sum);
#endif
	}




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
	      F2  = p_potinfo->get_abop_omega(typei, typej, typek);
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
	if (pair_ik_perriot){
	  double pref = threebodyfactor * fcik * gijk * F1 * F2 ;
	  force_ABOP_perriot_K_frc(i, k, Kik, dposik, pref);
	  
	  
#if 0
	  // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
	  // Loop over s
	  // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
	  Tiks_n_sum=0.0;
	  for (is=0; is<neighborcollection[i].size(); is++){
	    s = neighborcollection[i][is];
	    if (s==i) continue;
	    if (s==k) continue;
	    types = elem.name2idx( matter[s] );
	    if (! iac_pure_ABOP) if (p_potinfo->basepot(typei, types) != "ABOP") continue;
	    if (! iac_pure_ABOP) if (p_potinfo->basepot(typek, types) != "ABOP") continue;

	    ivecis = se_ivecij;
	    if (! sys_single_elem) ivecis = p_potinfo->basepot_vecidx(typei, types);
	    if (ivecis<0) continue;

	    ivecks = se_ivecij;
	    if (! sys_single_elem) ivecks = p_potinfo->basepot_vecidx(typek, types);
	    if (ivecks<0) continue;
	
	    get_atom_distance_vec(pos[i], pos[s], dposis);
	    ris = dposis.magn();
	    get_atom_distance_vec(pos[k], pos[s], dposks);
	    rks = dposks.magn();

	    pair_is_perriot = false;
	    if (p_potinfo->pot_ABOP[ivecis].rcut_fun=="perriot"){
	      pair_is_perriot = true;
	      rcutis = p_potinfo->pot_ABOP[ivecis].parname2val("prcut");
	      mis = p_potinfo->pot_ABOP[ivecis].parname2val("pm");
	    }
	    pair_ks_perriot = false;
	    if (p_potinfo->pot_ABOP[ivecks].rcut_fun=="perriot"){
	      pair_ks_perriot = true;
	      rcutks = p_potinfo->pot_ABOP[ivecks].parname2val("prcut");
	      mks = p_potinfo->pot_ABOP[ivecks].parname2val("pm");
	    }

	    if (pair_is_perriot==false || pair_ks_perriot==false)
	      continue;
	    
	    if (ris > rcutis) continue;
	    if (rks > rcutks) continue;

	    Xis = ris/(1.0 - pow(ris/rcutis, mis));
	    Xks = rks/(1.0 - pow(rks/rcutks, mks));

	    // *****************************************************
	    if (Xis + Xks > 3.0*rik) continue;
	    // *****************************************************

	    Tiks = 0.0;
	    if (Xis + Xks < 3.0*rik){
	      Tiks = -0.5 + rik/(Xis + Xks - rik);
	    }




	    
	    dXis = Xis/ris + mis * (Xis/ris)*(Xis/ris) * pow(ris/rcutis, mis);
	    dXks = Xks/rks + mks * (Xks/rks)*(Xks/rks) * pow(rks/rcutks, mks);

	    td1 = 1.0/(Xis + Xks - rik);
	    td2 = - rik * td1 * td1;

	    for (p=0; p<3; ++p){
	      dTiks_ik[p] = td1 *   dposik[p]/rik  + td2 * ( (-1)*  dposik[p]/rik  );
	      dTiks_is[p] = td2 * ( dXis *   dposis[p]/ris );
	      
	      dTiks_ki[p] = - dTiks_ik[p];
	      //dTijk_ji[p] = td1 * (-dposij[p]/rij) + td2 * ( (-1)*(-dposij[p]/rij) );
	      dTiks_ks[p] = td2 * ( dXks *   dposks[p]/rks );

	      dTiks_si[p] = - dTiks_is[p];
	      dTiks_sk[p] = - dTiks_ks[p];
	      //dTijk_ki[p] = td2 * ( dXik * (-dposik[p]/rik));
	      //dTijk_kj[p] = td2 * ( dXjk * (-dposjk[p]/rjk));
	    }

	    td = nik * pow(Tiks, nik-1.0);
	    for (p=0; p<3; ++p){
	      dKik_ik[p] = - Kik * td * (dTiks_ik[p]);
	      dKik_is[p] = - Kik * td * (dTiks_is[p]);
	      
	      dKik_ki[p] = - Kik * td * (dTiks_ki[p]);
	      dKik_ks[p] = - Kik * td * (dTiks_ks[p]);
	      
	      dKik_si[p] = - Kik * td * (dTiks_si[p]);
	      dKik_sk[p] = - Kik * td * (dTiks_ks[p]);
	    }
	  

	    // Establish pairwise forces:
	    for (p=0; p<3; ++p){
	      frc_ik[p] = threebodyfactor * dKik_ik[p] * fcik * gijk * F1 * F2 ;
	      frc_is[p] = threebodyfactor * dKik_is[p] * fcik * gijk * F1 * F2 ;
	      frc_ks[p] = threebodyfactor * dKik_ks[p] * fcik * gijk * F1 * F2 ;
	    }
	    // Add pairwise contributions to virials:
	    for (int v1=0; v1<3; ++v1){
	      for (int v2=0; v2<3; ++v2){
		virials[i].elem(v1,v2) += 0.5 *   frc_ik[v1]  *   dposik[v2];
		virials[k].elem(v1,v2) += 0.5 * (-frc_ik[v1]) * (-dposik[v2]);
		
		virials[i].elem(v1,v2) += 0.5 *   frc_is[v1]  *   dposis[v2];
		virials[s].elem(v1,v2) += 0.5 * (-frc_is[v1]) * (-dposis[v2]);
		
		virials[k].elem(v1,v2) += 0.5 *   frc_ks[v1]  *   dposks[v2];
		virials[s].elem(v1,v2) += 0.5 * (-frc_ks[v1]) * (-dposks[v2]);
	      }
	    }

	    // Add pairwise forces to total foces:
	    for (p=0; p<3; ++p){
	      frc[i][p] +=   frc_ik[p] + frc_is[p];
	      frc[k][p] += - frc_ik[p] + frc_ks[p];
	      frc[s][p] +=  -frc_is[p] - frc_ks[p];
	    }


	  }
#endif

	}




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
      } // end of loop over neighbors k
      
    } // end of loop over neighbors j
    //std::cout << "Atom i energy " << Ep[i] << std::endl;
  } // end of loop over atoms
  



  return Ep_tot_local;
}





void MDSystem::force_ABOP_perriot_K(int i,
				    int j,
				    double & Kij,
				    Vector3<double> dposij
				    ){
  
  int k, typei, typej, typek;
  int se_ivecij, ivecij, ivecjk, ivecik;
  bool pair_ik_perriot, pair_jk_perriot;

  typei = elem.name2idx( matter[i] );
  typej = elem.name2idx( matter[j] );
  typek = elem.name2idx( matter[k] );
  se_ivecij = p_potinfo->basepot_vecidx(typei, typei);

  Vector3<double> dposik, dposjk;
  double rik, rjk, rcutik, rcutjk, mik, mjk, Xik, Xjk, Tijk_n_sum, Tijk, nij;
  double rij, td, td1, td2;



  nij = p_potinfo->pot_ABOP[ivecij].parname2val("pn");
  rij = dposij.magn();



  // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
  // Loop over k
  // kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
  Tijk_n_sum=0.0;
  for (int ik=0; ik<neighborcollection[i].size(); ik++){
    k = neighborcollection[i][ik];
    if (k==i) continue;
    if (k==j) continue;

    typek = elem.name2idx( matter[k] );
    if (! iac_pure_ABOP) if (p_potinfo->basepot(typei, typek) != "ABOP") continue;
    if (! iac_pure_ABOP) if (p_potinfo->basepot(typej, typek) != "ABOP") continue;

    ivecik = se_ivecij;
    if (! sys_single_elem) ivecik = p_potinfo->basepot_vecidx(typei, typek);
    if (ivecik<0) continue;

    ivecjk = se_ivecij;
    if (! sys_single_elem) ivecjk = p_potinfo->basepot_vecidx(typej, typek);
    if (ivecjk<0) continue;
	
    get_atom_distance_vec(pos[i], pos[k], dposik);
    rik = dposik.magn();
    get_atom_distance_vec(pos[j], pos[k], dposjk);
    rjk = dposjk.magn();
    
    pair_ik_perriot = false;
    if (p_potinfo->pot_ABOP[ivecik].rcut_fun=="perriot"){
      pair_ik_perriot = true;
      rcutik = p_potinfo->pot_ABOP[ivecik].parname2val("prcut");
      mik = p_potinfo->pot_ABOP[ivecik].parname2val("pm");
    }

    pair_jk_perriot = false;
    if (p_potinfo->pot_ABOP[ivecjk].rcut_fun=="perriot"){
      pair_jk_perriot = true;
      rcutjk = p_potinfo->pot_ABOP[ivecjk].parname2val("prcut");
      mjk = p_potinfo->pot_ABOP[ivecjk].parname2val("pm");
    }

    if (pair_ik_perriot==false || pair_jk_perriot==false)
      continue;

    if (rik > rcutik) continue;
    if (rjk > rcutjk) continue;

    Xik = rik/(1.0 - pow(rik/rcutik, mik));
    Xjk = rjk/(1.0 - pow(rjk/rcutjk, mjk));

    // *****************************************************
    if (Xik + Xjk > 3.0*rij) continue;
    // *****************************************************

    Tijk = 0.0;
    if (Xik + Xjk < 3.0*rij){
      Tijk = -0.5 + rij/(Xik + Xjk - rij);
    }
    td = 0.0;
    if (Tijk>0.0) td = pow(Tijk, nij);
    Tijk_n_sum += td;
  }
  Kij = exp( - Tijk_n_sum);

}




void MDSystem::force_ABOP_perriot_K_frc(int i,
					int j,
					double & Kij,
					Vector3<double> dposij,
					double pref
					){
  
  force_ABOP_perriot_K(i, j, Kij, dposij);


  Vector3<double> dTijk_ij, dTijk_ik, dTijk_jk, dTijk_ji, dTijk_ki, dTijk_kj;
  Vector3<double> dKij_ij, dKij_ik, dKij_jk, dKij_ji, dKij_ki, dKij_kj;
  Vector3<double> frc_ij, frc_ik, frc_jk;
  
  
  Tijk_n_sum=0.0;
  for (int ik=0; ik<neighborcollection[i].size(); ik++){
    k = neighborcollection[i][ik];
    if (k==i) continue;
    if (k==j) continue;

    typek = elem.name2idx( matter[k] );
    if (! iac_pure_ABOP) if (p_potinfo->basepot(typei, typek) != "ABOP") continue;
    if (! iac_pure_ABOP) if (p_potinfo->basepot(typej, typek) != "ABOP") continue;

    ivecik = se_ivecij;
    if (! sys_single_elem) ivecik = p_potinfo->basepot_vecidx(typei, typek);
    if (ivecik<0) continue;

    ivecjk = se_ivecij;
    if (! sys_single_elem) ivecjk = p_potinfo->basepot_vecidx(typej, typek);
    if (ivecjk<0) continue;
	
    get_atom_distance_vec(pos[i], pos[k], dposik);
    rik = dposik.magn();
    get_atom_distance_vec(pos[j], pos[k], dposjk);
    rjk = dposjk.magn();
    
    pair_ik_perriot = false;
    if (p_potinfo->pot_ABOP[ivecik].rcut_fun=="perriot"){
      pair_ik_perriot = true;
      rcutik = p_potinfo->pot_ABOP[ivecik].parname2val("prcut");
      mik = p_potinfo->pot_ABOP[ivecik].parname2val("pm");
    }

    pair_jk_perriot = false;
    if (p_potinfo->pot_ABOP[ivecjk].rcut_fun=="perriot"){
      pair_jk_perriot = true;
      rcutjk = p_potinfo->pot_ABOP[ivecjk].parname2val("prcut");
      mjk = p_potinfo->pot_ABOP[ivecjk].parname2val("pm");
    }

    if (pair_ik_perriot==false || pair_jk_perriot==false)
      continue;

    if (rik > rcutik) continue;
    if (rjk > rcutjk) continue;

    Xik = rik/(1.0 - pow(rik/rcutik, mik));
    Xjk = rjk/(1.0 - pow(rjk/rcutjk, mjk));

    // *****************************************************
    if (Xik + Xjk > 3.0*rij) continue;
    // *****************************************************

    Tijk = 0.0;
    if (Xik + Xjk < 3.0*rij){
      Tijk = -0.5 + rij/(Xik + Xjk - rij);
    }
    td = 0.0;
    if (Tijk>0.0) td = pow(Tijk, nij);
    Tijk_n_sum += td;
  

    dXik = Xik/rik + mik * (Xik/rik)*(Xik/rik) * pow(rik/rcutik, mik);
    dXjk = Xjk/rjk + mjk * (Xjk/rjk)*(Xjk/rjk) * pow(rjk/rcutjk, mjk);

    td1 = 1.0/(Xik + Xjk - rij);
    td2 = - rij * td1 * td1;




    for (int p=0; p<3; ++p){
      dTijk_ij[p] = td1 *   dposij[p]/rij  + td2 * ( (-1)*  dposij[p]/rij  );
      dTijk_ik[p] = td2 * ( dXik *   dposik[p]/rik );

      dTijk_ji[p] = - dTijk_ij[p];
      //dTijk_ji[p] = td1 * (-dposij[p]/rij) + td2 * ( (-1)*(-dposij[p]/rij) );
      dTijk_jk[p] = td2 * ( dXjk *   dposjk[p]/rjk );

      dTijk_ki[p] = - dTijk_ik[p];
      dTijk_kj[p] = - dTijk_jk[p];
      //dTijk_ki[p] = td2 * ( dXik * (-dposik[p]/rik));
      //dTijk_kj[p] = td2 * ( dXjk * (-dposjk[p]/rjk));
    }

    td = nij * pow(Tijk, nij-1.0);
    for (int p=0; p<3; ++p){
      dKij_ij[p] = - Kij * td * (dTijk_ij[p]);
      dKij_ik[p] = - Kij * td * (dTijk_ik[p]);
      
      dKij_ji[p] = - Kij * td * (dTijk_ji[p]);
      dKij_jk[p] = - Kij * td * (dTijk_jk[p]);
      
      dKij_ki[p] = - Kij * td * (dTijk_ki[p]);
      dKij_kj[p] = - Kij * td * (dTijk_kj[p]);
    }



    // Kij: prefactor = - 0.5 * (VRij - bij * VAij) ;
    // Kik: prefactor = threebodyfactor * fcik * gijk * F1 * F2 ;


    // Use previously calc. results to get total pairwise forces:
    for (int p=0; p<3; ++p){
      frc_ij[p] = - 0.5 * dKij_ij[p] * (VRij - bij * VAij) ;
      frc_ik[p] = - 0.5 * dKij_ik[p] * (VRij - bij * VAij) ;
      frc_jk[p] = - 0.5 * dKij_jk[p] * (VRij - bij * VAij) ;
    }

    // Add to total atomic forces for all atoms participating:
    for (int p=0; p<3; ++p){
      frc[i][p] +=   frc_ij[p] + frc_ik[p];
      frc[j][p] += - frc_ij[p] + frc_jk[p];
      frc[k][p] +=  -frc_ik[p] - frc_jk[p];
      /*
	frc[i][p] += - 0.5 * (dKij_ij[p] + dKij_ik[p]) * (VRij - bij * VAij) ;
	frc[j][p] += - 0.5 * (dKij_ji[p] + dKij_jk[p]) * (VRij - bij * VAij) ;
	frc[k][p] += - 0.5 * (dKij_ki[p] + dKij_kj[p]) * (VRij - bij * VAij) ;
      */
    }

    // Add to total atomic virials for all atoms participating:
    for (int v1=0; v1<3; ++v1){
      for (int v2=0; v2<3; ++v2){
	// ij
	virials[i].elem(v1,v2) += 0.5 *   frc_ij[v1]  *   dposij[v2];
	virials[j].elem(v1,v2) += 0.5 * (-frc_ij[v1]) * (-dposij[v2]);
	// ik
	virials[i].elem(v1,v2) += 0.5 *   frc_ik[v1]  *   dposik[v2];
	virials[k].elem(v1,v2) += 0.5 * (-frc_ik[v1]) * (-dposik[v2]);
	// jk
	virials[j].elem(v1,v2) += 0.5 *   frc_jk[v1]  *   dposjk[v2];
	virials[k].elem(v1,v2) += 0.5 * (-frc_jk[v1]) * (-dposjk[v2]);
	
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
    
  }



}

