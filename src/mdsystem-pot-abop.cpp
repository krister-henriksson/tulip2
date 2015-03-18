


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
  int se_ivecik = se_ivecij;
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




  int ij,j,ik,k,p,q;
  int typei, typej, typek;
  double D0ij,r0ij,betaij,Sij,Rij,Dij,fcij,rij,rcutij;
  double gammaik,cik,dik,hik,r0ik,Rik,Dik,fcik,rik,rcutik;
  double V1, dV1, fermi, dfermi, bfermi, rfermi;
  double VRij,VAij,bij,Chiij,gijk,c2,d2,cost,hcost,hcost2,alphaijk,twomuik;

  //double eps = std::numeric_limits<double>::epsilon();
  double VRij_r, VAij_r, dVRij_r, dVAij_r;
  double dVRij, dVAij, dgijk;
  double dfcij, dfcik, threebodyfactor;

  Vector3<double> dposij,dposik;
  Vector3<double> dcost_i, dcost_j, dcost_k;
  Vector3<double> dgijk_i, dgijk_j, dgijk_k;
  Vector3<double> frci, frcj, frck;

  double td, td1,td2;
  int ivecij, ivecik, ivec_reppot, Nr;
  double Epij;


  double frc_ij, frc_ik, frc_jk;

  double F1, F2, dF1, dF2;
  Vector3<double> dF1_i, dF1_j, dF1_k;
  Vector3<double> dF2_i, dF2_j, dF2_k;




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
      
      if (! iac_pure_ABOP)
	if (p_potinfo->basepot(typei,typej) != "ABOP") continue;
      //if (p_potinfo->basepot(s1,s2) != "ABOP") continue;

      ivecij = se_ivecij;
      if (! sys_single_elem)
	ivecij = p_potinfo->basepot_vecidx(typei, typej);
      //ivecij = p_potinfo->basepot_vecidx(s1,s2);
      if (ivecij<0) continue;
      




      get_atom_distance_vec(pos[i], pos[j], dposij);
      rij = dposij.magn();

      Rij = p_potinfo->pot_ABOP[ivecij].parval[8];
      Dij = p_potinfo->pot_ABOP[ivecij].parval[9];
      rcutij = Rij+Dij;
      if (rij > rcutij) continue;
      
      D0ij   = p_potinfo->pot_ABOP[ivecij].parval[0];
      r0ij   = p_potinfo->pot_ABOP[ivecij].parval[1];
      betaij = p_potinfo->pot_ABOP[ivecij].parval[2];
      Sij    = p_potinfo->pot_ABOP[ivecij].parval[3];
      /*
      std::cout << "D0ij " << D0ij << " r0ij " << r0ij << " betaij " << betaij << " Sij " << Sij
	   << " gammaij " << gammaij << " cij " << cij << " dij " << dij << " hij " << hij
	   << " Rij " << Rij << " Dij " << Dij << std::endl;
      */

      if (rij < Rij-Dij){
	fcij = 1.0;
	dfcij = 0.0;
      }
      else if (rij > rcutij){
	fcij=0;
	dfcij = 0.0;
	continue;
      }
      else {
	td = 0.5*PI*(rij-Rij)/Dij;

	fcij  = 0.5 * (1.0 - sin( td ));
	dfcij = -0.5 * 0.5*PI/Dij * cos( td );
      }
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

	bfermi = p_potinfo->pot_Reppot[ivec_reppot].bfermi;
	rfermi = p_potinfo->pot_Reppot[ivec_reppot].rfermi;
	  
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

	if (! iac_pure_ABOP)
	  if (p_potinfo->basepot(typei, typek) != "ABOP") continue;
	//if (p_potinfo->basepot(s1,s3) != "ABOP") continue;

	ivecik = se_ivecik;
	if (! sys_single_elem)
	  ivecik = p_potinfo->basepot_vecidx(typei, typek);
	//ivecik = p_potinfo->basepot_vecidx(s1,s3);
	if (ivecik<0) continue;
	
	get_atom_distance_vec(pos[i], pos[k], dposik);
	rik = dposik.magn();

	r0ik = p_potinfo->pot_ABOP[ivecik].parval[1];

	Rik = p_potinfo->pot_ABOP[ivecik].parval[8];
	Dik = p_potinfo->pot_ABOP[ivecik].parval[9];
	rcutik = Rik+Dik;

	if (rik > rcutik) continue;

	gammaik = p_potinfo->pot_ABOP[ivecik].parval[4];
	cik = p_potinfo->pot_ABOP[ivecik].parval[5];
	dik = p_potinfo->pot_ABOP[ivecik].parval[6];
	hik = p_potinfo->pot_ABOP[ivecik].parval[7];
	/*
	std::cout << "D0ik " << D0ik << " r0ik " << r0ik << " betaik " << betaik << " Sik " << Sik
	     << " gammaik " << gammaik << " cik " << cik << " dik " << dik << " hik " << hik
	     << " Rik " << Rik << " Dik " << Dik << std::endl;
	*/

	if (rik < Rik-Dik){
	  fcik  = 1.0;
	  dfcik = 0.0;
	}
	else if (rik > rcutik){
	  fcik  = 0;
	  dfcik = 0.0;
	  continue;
	}
	else {
	  td = 0.5*PI*(rik-Rik)/Dik;

	  fcik  = 0.5 * (1.0 - sin( td ));
	  dfcik = -0.5 * 0.5*PI/Dik * cos( td );
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
	    double td = rij - r0ij - (rik - r0ik);
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

	
	Chiij += fcik * gijk * F1 * F2;
      } // end of loop over neighbors k
      bij = 1.0/sqrt(1.0 + Chiij);
      //std::cout << "bij " << bij << std::endl;


      // ################################################################
      // Energy:
      // ################################################################
      Epij = 0.5 * fcij * (VRij - bij * VAij);
      Ep[i] += Epij;
      Ep_tot_local += Epij;




      // ################################################################
      // Two-body forces:
      // ################################################################
      frc_ij = - 0.5 *
	  ( dfcij * VRij
	    + fcij * dVRij
	    - bij * dfcij * VAij
	    - bij * fcij * dVAij);

      for (p=0; p<3; ++p){
	frci[p] = frc_ij *   dposij[p]/rij;
	frcj[p] = frc_ij * (-dposij[p]/rij);
      }


#if 0
      for (p=0; p<3; ++p){
	frci[p] = - 0.5 *
	  ( dfcij * VRij * dposij[p]/rij
	    + fcij * dVRij * dposij[p]/rij
	    - bij * dfcij * VAij * dposij[p]/rij
	    - bij * fcij * dVAij * dposij[p]/rij
	    );

	frcj[p] = - 0.5 *
	  ( dfcij * VRij * (-dposij[p]/rij)
	    + fcij * dVRij * (-dposij[p]/rij)
	    - bij * dfcij * VAij * (-dposij[p]/rij)
	    - bij * fcij * dVAij * (-dposij[p]/rij)
	    );
      }
#endif

      for (int v1=0; v1<3; ++v1){
	for (int v2=0; v2<3; ++v2){
	  virials[i].elem(v1,v2) += (frc_ij * dposij[v1]/rij) * dposij[v2];
	  //virials[j].elem(v1,v2) += frcj[v1] * pos[j][v2];
	}
      }


      for (p=0; p<3; ++p){
	frc[i][p] += frci[p];
	frc[j][p] += frcj[p];
      }




      // ################################################################
      // Three-body forces:
      // ################################################################
      
      threebodyfactor = - 0.25 * fcij * VAij * bij*bij*bij;


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

	ivecik = se_ivecik;
	if (! sys_single_elem)
	  ivecik = p_potinfo->basepot_vecidx(typei, typek);
	//ivecik = p_potinfo->basepot_vecidx(s1,s3);
	if (ivecik<0) continue;
	
	get_atom_distance_vec(pos[i], pos[k], dposik);
	rik = dposik.magn();

	r0ik = p_potinfo->pot_ABOP[ivecik].parval[1];

	Rik = p_potinfo->pot_ABOP[ivecik].parval[8];
	Dik = p_potinfo->pot_ABOP[ivecik].parval[9];
	rcutik = Rik+Dik;
	if (rik > rcutik) continue;
      
	gammaik = p_potinfo->pot_ABOP[ivecik].parval[4];
	cik = p_potinfo->pot_ABOP[ivecik].parval[5];
	dik = p_potinfo->pot_ABOP[ivecik].parval[6];
	hik = p_potinfo->pot_ABOP[ivecik].parval[7];

	if (rik < Rik-Dik){
	  fcik  = 1.0;
	  dfcik = 0.0;
	}
	else if (rik > rcutik){
	  fcik  = 0;
	  dfcik = 0.0;
	  continue;
	}
	else {
	  td = 0.5*PI*(rik-Rik)/Dik;

	  fcik  = 0.5 * (1.0 - sin( td ));
	  dfcik = -0.5 * 0.5*PI/Dik * cos( td );
	}
	
	c2 = cik*cik;
	d2 = dik*dik;
	cost = (dposij * dposik) / (rij*rik);
	hcost  = hik + cost;
	hcost2 = hcost * hcost;

	gijk = gammaik * (1.0 + c2/d2 - c2 / (d2 + hcost2) );



	for (p=0; p<3; ++p){
	  dcost_i[p] = - cost/(rij*rij) *   dposij[p]
	    - cost/(rik*rik) * dposik[p]
	    + 1.0/(rij*rik) * (dposij[p] + dposik[p]);
	  
	  dcost_j[p] = - cost/(rij*rij) * (-dposij[p])
	    + 0.0
	    + 1.0/(rij*rik) * (-dposik[p]);

	  dcost_k[p] = - cost/(rik*rik) * (-dposik[p])
	    + 0.0
	    + 1.0/(rij*rik) * (-dposij[p]);
	}

	dgijk = gammaik * c2 /((d2 + hcost2)*(d2 + hcost2)) * 2 * hcost;
	for (p=0; p<3; ++p){
	  dgijk_i[p] = dgijk * dcost_i[p];
	  dgijk_j[p] = dgijk * dcost_j[p];
	  dgijk_k[p] = dgijk * dcost_k[p];
	}






	dF1  = 0.0;
	for (p=0; p<3; ++p){
	  dF1_i[p] = 0.0;
	  dF1_j[p] = 0.0;
	  dF1_k[p] = 0.0;
	}
	dF2 = 0.0;
	for (p=0; p<3; ++p){
	  dF2_i[p] = 0.0;
	  dF2_j[p] = 0.0;
	  dF2_k[p] = 0.0;
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
	    dF1_i[p] = dF1 * ( dposij[p]/rij - dposik[p]/rik );
	    dF1_j[p] = dF1 * (-dposij[p]/rij);
	    dF1_k[p] = dF1 * ( dposik[p]/rik);
	  }
	  // .........................................................
	  if (p_potinfo->use_abop_omega.elem(typei, typej, typek)){
	    if (! sys_single_elem){
	      F2  = p_potinfo->get_abop_omega(typei, typej, typek);
	      dF2 = 0.0;
	    }
	  }
	  else {
	    double td = rij - r0ij - (rik - r0ik);
	    F2  = exp( alphaijk * td );
	    dF2 = alphaijk * F2;
	    for (p=0; p<3; ++p){
	      dF2_i[p] = dF2 * ( dposij[p]/rij - dposik[p]/rik );
	      dF2_j[p] = dF2 * (-dposij[p]/rij);
	      dF2_k[p] = dF2 * ( dposik[p]/rik);
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
	    dF1_i[p] = dF1 * ( dposij[p]/rij - dposik[p]/rik );
	    dF1_j[p] = dF1 * (-dposij[p]/rij);
	    dF1_k[p] = dF1 * ( dposik[p]/rik);
	  }
	  F2  = 1.0;
	  dF2 = 0.0;
	}
	// *****************************************************************


	//std::cout << "alpha omega  " << p_potinfo->abop_alpha.elem(typei, typej, typek) << " " << omegaijk << std::endl;

	for (p=0; p<3; ++p){
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
	}




	
	for (int v1=0; v1<3; ++v1){
	  for (int v2=0; v2<3; ++v2){
	    // ij contribution:
	    virials[i].elem(v1,v2) += 
	      ( 
	       threebodyfactor * fcik * dgijk * ( dposik[v1]/(rij*rik) - cost*dposij[v1]/(rij*rij))
	       * F1 * F2
	       +
	       threebodyfactor * fcik * gijk * dF1 * dposij[v1]/rij * F2
	       +
	       threebodyfactor * fcik * gijk * F1 * dF2 * dposij[v1]/rij
	        )
	      * dposij[v2];
	    // ik contribution:
	    virials[i].elem(v1,v2) += 
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
	  }
	}
	

	for (p=0; p<3; ++p){
	  frc[i][p] += frci[p];
	  frc[j][p] += frcj[p];
	  frc[k][p] += frck[p];
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



