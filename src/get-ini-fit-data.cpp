
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
#include "funcfit-simplexfit.hpp"
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
#include "compound.hpp"
#include "elem-iacs.hpp"
#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mtwister.hpp"
#include "param-pot.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "potinfo.hpp"
#include "propfun.hpp"
#include "specs-fit-prop-pot.hpp"
#include "get-comp-prop.hpp"
#include "errors.hpp"

#include "get-ini-fit-data.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;






using utils::Vector;
using namespace constants;
using namespace physconst;
using boost::format;





void get_ini_fit_data(ParamPot & param,
		      Vector<CompoundStructureFit> & DX
		      ){

  CompoundStructureFit cmpfit;
  int iDX,sizeDX,i,j,k,p,q;


  // Make sure potentials have been updated from the fitting parameters:
  param.update_pot();

  sizeDX = DX.size();

  if (sizeDX==0){
    std::cout << "No compounds to inspect!" << std::endl; return;
  }


  Vector<std::string> elnames = DX[0].elemnames;
  int nelem = DX[0].nelem;
  std::string name1, name2;

  k=0;
  for (i=0; i<nelem; ++i){
    for (j=i; j<nelem; ++j){
      if (param.p_potinfo->is_fittable(elnames[i], elnames[j])){
	k++;
	name1 = elnames[i];
	name2 = elnames[j];
      }
    }
  }
  if (k>1)
    aborterror("Error: Initial calculations for potential fitting only works"
	       " for a single pair of elements! Exiting.");

  if (param.p_potinfo->basepot(name1,name2) != "ABOP")
    aborterror("Error: Initial calculations for potential fitting only works"
	       " for ABOP type potentials for now! Exiting.");


  int typei = param.p_potinfo->elem.name2idx( name1 );
  int typej = param.p_potinfo->elem.name2idx( name2 );

  std::cout << "********************************************************************" << std::endl;

  std::cout << "Using elements " << name1 << "-" << name2
	    << " to make some initial fitting calculations ..." << std::endl;

  int ivecij = param.p_potinfo->basepot_vecidx(typei, typej);
  double rcutij = param.p_potinfo->pot_ABOP[ivecij].rcut();
  std::cout << "rcutij: " << rcutij << std::endl;

  int ivecii = param.p_potinfo->basepot_vecidx(typei, typei);
  double rcutii = param.p_potinfo->pot_ABOP[ivecii].rcut();
  std::cout << "rcutii: " << rcutii << std::endl;

  int ivecjj = param.p_potinfo->basepot_vecidx(typej, typej);
  double rcutjj = param.p_potinfo->pot_ABOP[ivecjj].rcut();
  std::cout << "rcutjj: " << rcutjj << std::endl;

  std::string ts1, ts2;


  std::cout << "Made it here 2 001" << std::endl;

  Vector<double>        Ep_list(sizeDX);
  Vector<int>           nat1(sizeDX);
  Vector<int>           nat2(sizeDX);
  Vector<int>           nbonds_list(sizeDX);
  Vector<int>           nbonds_list12(sizeDX);
  Vector<int>           nbonds_list21(sizeDX);
  Vector<BondData>      bond_list(sizeDX);
  Vector<BondData>      bond_list12(sizeDX);
  Vector<BondData>      bond_list21(sizeDX);
  Vector<BondAngleData> bondangle_list(sizeDX);



  std::cout << "Made it here 2 002" << std::endl;


  // Loop through compounds:
  for (iDX=0; iDX<sizeDX; ++iDX){
    // *******************************************************
    // Get potential energy
    // *******************************************************
    cmpfit = DX[iDX];
    
    MDSystem mds;

    // Set potential:
    mds.p_potinfo = param.p_potinfo;

    // Set other settings:
    mds.elem = param.p_potinfo->elem;
    mds.iacs = param.p_potinfo->iacs;
    mds.name = cmpfit.name;

    // Initialize common MD settings for all compounds:
    mds.specs_common = param.p_potinfo->specs_prop.mds_specs_common;
    // Initialize MD settings:
    mds.specs        = cmpfit.mds_specs;

    // Other MD settings from different places ...:
    mds.omp_info = param.p_potinfo->omp_info;
    mds.rcut_max = param.p_potinfo->get_rcut_max( cmpfit.elemnames );
    mds.rcut     = mds.rcut_max;
    mds.skint    = mds.specs.skint;
    
    // These settings are applied to inherited base class members:
    mds.use_def_xyz_fmt = mds.specs_common.use_def_dump_xyz_fmt;
    mds.def_xyz_fmt     = mds.specs_common.def_dump_xyz_fmt;
    
    // Create mds:
    double rm = mds.rcut + mds.skint;
    mds.create_from_structure(cmpfit, 2.0*rm); // removes old atoms
    // Single time step relax:
    mds.specs.dt = 1e-7;
    mds.specs.tend = mds.specs.tstart = 0.0;
    mds.relax();


    std::cout << "Made it here 2 003" << std::endl;

    // std::cout << "1 mds.rcut " << mds.rcut << " and mds.rcut_max " << mds.rcut_max << std::endl;

    // *******************************************************
    // Get bond distances and neighbor counts
    // *******************************************************
    // Bond distances d_i and their number n_i
    // Vbond will be calculated from Ec and number of bonds

    // std::cout << "2 mds.rcut " << mds.rcut << " and mds.rcut_max " << mds.rcut_max << std::endl;

    mds.get_bond_list( bond_list[iDX],
		       bond_list12[iDX],
		       bond_list21[iDX],
		       name1, name2,
		       nat1[iDX], nat2[iDX],
		       nbonds_list[iDX],
		       nbonds_list12[iDX],
		       nbonds_list21[iDX],
		       rcutij );
    //mds.get_bond_list( bond_dist_list[iDX], bond_num_list[iDX], name1[i], name2[j], n12 );

    std::cout << "Made it here 2 004" << std::endl;

    // std::cout << "3 mds.rcut " << mds.rcut << " and mds.rcut_max " << mds.rcut_max << std::endl;
    mds.get_bond_angle_list( bondangle_list[iDX], name1, name2, rcutii, rcutjj, rcutij );

    std::cout << "Made it here 2 005" << std::endl;


    // Remove potential energy contributions from known interactions:
    int nat = mds.natoms();
    double Ep = mds.calc_potential_energy();
    std::cout << "Energy from system: " << Ep << std::endl;
    if (cmpfit.prop_use.Ecoh)
      Ep_list[iDX] = nat * cmpfit.prop_readin.Ecoh - Ep;
    else if (cmpfit.prop_use.Emix){
      // Energy of formation: cmpfit.prop_readin.Emix * nat
      // Ef = Etot - n1*Ecoh1 - n2*Ecoh2 - ...
      // => Etot = Ef + sum_i ni * Ecohi
      Ep_list[iDX] = nat * cmpfit.prop_readin.Emix
	+ nat1[iDX] * param.p_potinfo->Ecoh_ref[ typei ]
	+ nat2[iDX] * param.p_potinfo->Ecoh_ref[ typej ]
	- Ep;
    }
    else if (cmpfit.prop_use.Eform){
      // Energy of formation: cmpfit.prop_readin.Emix * nat
      // Ef = Etot - n1*Ecoh1 - n2*Ecoh2 - ...
      // => Etot = Ef + sum_i ni * Ecohi
      Ep_list[iDX] = cmpfit.prop_readin.Eform
	+ nat1[iDX] * param.p_potinfo->Ecoh_ref[ typei ]
	+ nat2[iDX] * param.p_potinfo->Ecoh_ref[ typej ]
	- Ep;
    }
    else
      aborterror("Make initial fit: Error: Ecoh or Emix or Eform not used for compound "
		 + cmpfit.name + "! Cannot obtain readin potential energy! Exiting.");






  }





  std::cout << "----------------------------------------------------------" << std::endl;

  Vector<std::string> tr(8);
  tr[0] = name1 + "-" + name1 + "-" + name1;
  tr[1] = name1 + "-" + name1 + "-" + name2;
  tr[2] = name1 + "-" + name2 + "-" + name1;
  tr[3] = name1 + "-" + name2 + "-" + name2;
  tr[4] = name2 + "-" + name1 + "-" + name1;
  tr[5] = name2 + "-" + name1 + "-" + name2;
  tr[6] = name2 + "-" + name2 + "-" + name1;
  tr[7] = name2 + "-" + name2 + "-" + name2;


  std::cout << "-------------------------" << std::endl;
  std::cout << "BOND ANGLES (cos(theta)):" << std::endl;
  std::cout << "-------------------------" << std::endl;

  for (int it=0; it<8; ++it){
    if (name1 == name2 && it>=1) break;

    std::cout << "Species triplet " << tr[it] << std::endl;
    for (iDX=0; iDX<sizeDX; ++iDX){
      std::string dumpfile("bondangles-" + DX[iDX].name + "-" + tr[it] + ".dat");
      std::ofstream fout;
      fout.open(dumpfile.c_str());

      std::cout << "Compound: " << format("%20s") % DX[iDX].name << std::endl;
      for (k=0; k<bondangle_list[iDX].costheta_ijk[it].size(); ++k){
	std::cout << format("%20.10f") % bondangle_list[iDX].costheta_ijk[it][k]  << "  "
	     << format("%20d")         % bondangle_list[iDX].ncostheta_ijk[it][k] << std::endl;
	fout << format("%20.10f")      % bondangle_list[iDX].costheta_ijk[it][k]  << "  "
	     << format("%20d")         % bondangle_list[iDX].ncostheta_ijk[it][k] << std::endl;
      }
      fout.close();
      fout.clear();
    }
  }

  std::cout << "----------------------------------------------------------" << std::endl;

  std::cout << "---------------" << std::endl;
  std::cout << "BOND DISTANCES:" << std::endl;
  std::cout << "---------------" << std::endl;

  std::string dumpfile("bonds-num.dat");
  std::ofstream fout;
  fout.open(dumpfile.c_str());

  ts1="Distance";
  ts2="Num_bonds";
  std::cout << "Species pair " << name1 << "-" << name2 << ":" << std::endl;
  for (iDX=0; iDX<sizeDX; ++iDX){
    std::cout << "Compound: " << format("%20s") % DX[iDX].name
	      << "  Potential energy to assign: " << format("%10.5e") % Ep_list[iDX] << std::endl;
    std::cout << format("%20s") % ts1 << "  " << format("%20s") % ts2 << std::endl;

    std::string dumpfile2("bonds-num-" + DX[iDX].name + ".dat");
    std::ofstream fout2;
    fout2.open(dumpfile2.c_str());


    if (name1 == name2){
      if (bond_list[iDX].dist.size()==0){
	std::cout << "No bond distances at all! Trying next compound ..." << std::endl;
	continue;
      }
      for (k=0; k<bond_list[iDX].dist.size(); ++k){
	if (bond_list[iDX].ndist[k]==0){
	  std::cout << "No bond distances of length " << bond_list[iDX].dist[k]
		    << " ! Trying next bond distance ..." << std::endl;
	  continue;
	}
	// nba = bond_list[iDX].ndist[k]/(0.5*(nat1[iDX]+nat2[iDX]));
	// nba = bond_list[iDX].ndist[k]/(nbonds_list[iDX]);
	std::cout << format("%20.10f") % bond_list[iDX].dist[k]  << "  "
		  << format("%20.10f") % bond_list[iDX].ndist[k] << std::endl;
	fout << format("%20.10f") % bond_list[iDX].dist[k] << "  "
	     << format("%20.10f") % bond_list[iDX].ndist[k]
	     << format("%20s")    % DX[iDX].name
	     << std::endl;
	fout2 << format("%20.10f") % bond_list[iDX].dist[k]  << "  "
	      << format("%20.10f") % bond_list[iDX].ndist[k] << std::endl;
      }
    }
    else {
      if (bond_list12[iDX].dist.size()==0 && bond_list21[iDX].dist.size()==0){
	std::cout << "No bond distances at all! Trying next compound ..." << std::endl;
	continue;
      }
      for (k=0; k<bond_list12[iDX].dist.size(); ++k){
	if (bond_list12[iDX].ndist[k]==0){
	  std::cout << "No bond distances of length " << bond_list12[iDX].dist[k]
		    << " ! Trying next bond distance ..." << std::endl;
	  continue;
	}
	std::cout << format("%20.10f") % bond_list12[iDX].dist[k]  << "  "
		  << format("%20.10f") % bond_list12[iDX].ndist[k] << std::endl;
	fout << format("%20.10f") % bond_list12[iDX].dist[k] << "  "
	     << format("%20.10f") % bond_list12[iDX].ndist[k]
	     << format("%20s")    % DX[iDX].name
	     << std::endl;
	fout2 << format("%20.10f") % bond_list12[iDX].dist[k]  << "  "
	      << format("%20.10f") % bond_list12[iDX].ndist[k] << std::endl;
      }
      for (k=0; k<bond_list21[iDX].dist.size(); ++k){
	if (bond_list21[iDX].ndist[k]==0){
	  std::cout << "No bond distances of length " << bond_list21[iDX].dist[k]
		    << " ! Trying next bond distance ..." << std::endl;
	  continue;
	}
	std::cout << format("%20.10f") % bond_list21[iDX].dist[k]  << "  "
		  << format("%20.10f") % bond_list21[iDX].ndist[k] << std::endl;
	fout << format("%20.10f") % bond_list21[iDX].dist[k] << "  "
	     << format("%20.10f") % bond_list21[iDX].ndist[k]
	     << format("%20s")    % DX[iDX].name
	     << std::endl;
	fout2 << format("%20.10f") % bond_list21[iDX].dist[k]  << "  "
	      << format("%20.10f") % bond_list21[iDX].ndist[k] << std::endl;
      }
    }




    fout2.close();
    fout2.clear();
  }
  fout.close();
  fout.clear();




  std::cout << "----------------------------------------------------------" << std::endl;


  std::cout << "---------------------------------" << std::endl;
  std::cout << "BOND DISTANCES AND BOND ENERGIES:" << std::endl;
  std::cout << "---------------------------------" << std::endl;

  for (iDX=0; iDX<sizeDX; ++iDX){
    if (name1==name2 && bond_list[iDX].dist.size()==0) continue;
    if (name1!=name2
	&& bond_list12[iDX].dist.size()==0
	&& bond_list21[iDX].dist.size()==0) continue;

    int nbonds = 0, nbonds_k=0;
    double rb = 0.0;

    if (name1==name2){
      for (k=0; k<bond_list[iDX].dist.size(); ++k){
	nbonds_k  = bond_list[iDX].ndist[k];
	nbonds   += nbonds_k;
	rb       += nbonds_k * bond_list[iDX].dist[k];
      }
    }
    else {
      for (k=0; k<bond_list12[iDX].dist.size(); ++k){
	nbonds_k  = bond_list12[iDX].ndist[k];
	nbonds   += nbonds_k;
	rb       += nbonds_k * bond_list12[iDX].dist[k];
      }
      for (k=0; k<bond_list21[iDX].dist.size(); ++k){
	nbonds_k  = bond_list21[iDX].ndist[k];
	nbonds   += nbonds_k;
	rb       += nbonds_k * bond_list21[iDX].dist[k];
      }
    }
    if (nbonds==0) continue;
    rb /= nbonds;
    double Vb = Ep_list[iDX] / nbonds;


    double nb;
    if (name1==name2)
      nb = nbonds_list[iDX] * 1.0/(nat1[iDX] + nat2[iDX]);
    else
      nb = (nbonds_list12[iDX] + nbonds_list21[iDX]) * 1.0/(nat1[iDX] + nat2[iDX]);

    /*
    if (name1==name2)
      printf("name1-name2: %s-%s  nbonds_list: %d  nat1: %d  nat2: %d\n",
	     name1.c_str(), name2.c_str(), nbonds_list[iDX], nat1[iDX], nat2[iDX]);
    else
      printf("name1-name2: %s-%s  nbonds_list12: %d  nbonds_list12: %d  nat1: %d  nat2: %d\n",
	     name1.c_str(), name2.c_str(), nbonds_list12[iDX], nbonds_list21[iDX], nat1[iDX], nat2[iDX]);
    */

    // std::cout << "sum_k bond_list[iDX].ndist[k] = " << nbonds << " nbonds_list[iDX] = " << nbonds_list[iDX] << std::endl;

    printf( "%20.10f  %20.10f  # Compound %20s : ave(rbond) Vbond   nbonds %d  ave(nbonds/nat_with_bonds) %20.10f\n",
	   rb, Vb, DX[iDX].name.c_str(), nbonds, nb );
  }


}




