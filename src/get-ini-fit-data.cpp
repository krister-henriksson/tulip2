
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
#include "funcfit-ls-gauss-newton.hpp"
#include "funcfit-ls-leve-marq.hpp"
#include "funcfit-ls-powelldogleg.hpp"
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

  int ivecij = param.p_potinfo->basepot_vecidx(typei, typej);
  double Rij = param.p_potinfo->pot_ABOP[ivecij].parval[8];
  double Dij = param.p_potinfo->pot_ABOP[ivecij].parval[9];
  double rcutij = Rij + Dij;

  int ivecii = param.p_potinfo->basepot_vecidx(typei, typei);
  double Rii = param.p_potinfo->pot_ABOP[ivecii].parval[8];
  double Dii = param.p_potinfo->pot_ABOP[ivecii].parval[9];
  double rcutii = Rii + Dii;

  int ivecjj = param.p_potinfo->basepot_vecidx(typej, typej);
  double Rjj = param.p_potinfo->pot_ABOP[ivecjj].parval[8];
  double Djj = param.p_potinfo->pot_ABOP[ivecjj].parval[9];
  double rcutjj = Rjj + Djj;

  std::string ts1, ts2;



  Vector<double>             Ep_list(sizeDX);
  Vector<int>                nat_with_bonds_list(sizeDX);
  Vector< Vector<BondData> > bond_list(sizeDX);
  Vector< Vector<BondAngleData> > bondangle_list(sizeDX);



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
    mds.omp_info     = param.p_potinfo->omp_info;
    mds.rcut         = mds.rcut_max = mds.p_potinfo->get_rcut_max( cmpfit.elemnames );
    mds.skint        = mds.specs.skint;
    
    // These settings are applied to inherited base class members:
    mds.use_def_xyz_fmt = param.p_potinfo->specs_prop.mds_specs_common.use_def_dump_xyz_fmt;
    mds.def_xyz_fmt     = param.p_potinfo->specs_prop.mds_specs_common.def_dump_xyz_fmt;
    
    // Create mds:
    double rm = mds.rcut + mds.skint;
    mds.create_from_structure(cmpfit, 2.0*rm); // removes old atoms
    // Single time step relax:
    mds.specs.dt = 1e-7;
    mds.specs.tend = mds.specs.tstart = 0.0;
    mds.relax();

    // Remove potential energy contributions from known interactions:
    int nat = mds.natoms();
    std::cout << "Energy from system: " << mds.calc_potential_energy()/nat << std::endl;
    Ep_list[iDX] = nat * cmpfit.prop_readin.Ecoh - mds.calc_potential_energy()/nat;


    // *******************************************************
    // Get bond distances and neighbor counts
    // *******************************************************
    // Bond distances d_i and their number n_i
    // Vbond will be calculated from Ec and number of bonds
    int nat_w_b=0;
    mds.get_bond_list( bond_list[iDX], name1, name2, nat_w_b, rcutij );
    //mds.get_bond_list( bond_dist_list[iDX], bond_num_list[iDX], name1[i], name2[j], n12 );

    nat_with_bonds_list[iDX] = nat_w_b;



    mds.get_bond_angle_list( bondangle_list[iDX], name1, name2, rcutii, rcutjj, rcutij );

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

  for (int it=0; it<8; ++it){
    if (name1 == name2 && it>=1) break;

    std::cout << "Species triplet " << tr[it] << std::endl;
    for (iDX=0; iDX<sizeDX; ++iDX){
      std::string dumpfile("bondangles-" + DX[iDX].name + "-" + tr[it] + ".dat");
      std::ofstream fout;
      fout.open(dumpfile.c_str());

      std::cout << "Compound: " << format("%20s") % DX[iDX].name << std::endl;
      for (k=0; k<bondangle_list[iDX].size(); ++k){
	std::cout << format("%20.10f") % bondangle_list[iDX][k].costheta_ijk[it] << "  "
	     << format("%20d") % bondangle_list[iDX][k].ntheta_ijk[it] << std::endl;
	fout << format("%20.10f") % bondangle_list[iDX][k].costheta_ijk[it] << "  "
	     << format("%20d") % bondangle_list[iDX][k].ntheta_ijk[it] << std::endl;
      }
      fout.close();
      fout.clear();
    }
  }

  std::cout << "----------------------------------------------------------" << std::endl;

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
    
    for (k=0; k<bond_list[iDX].size(); ++k){
      std::cout << format("%20.10f") % bond_list[iDX][k].dist << "  "
	   << format("%20d") % bond_list[iDX][k].nbonds << std::endl;
      fout << format("%20.10f") % bond_list[iDX][k].dist << "  "
	   << format("%20d") % bond_list[iDX][k].nbonds
	   << format("%20s") % DX[iDX].name
	   << std::endl;
      fout2 << format("%20.10f") % bond_list[iDX][k].dist << "  "
	    << format("%20d") % bond_list[iDX][k].nbonds
	    << std::endl;
    }
    fout2.close();
    fout2.clear();
  }
  fout.close();
  fout.clear();

  std::cout << "----------------------------------------------------------" << std::endl;

  for (iDX=0; iDX<sizeDX; ++iDX){

    int nbonds = 0, nbonds_k=0;
    double rb = 0.0;
    for (k=0; k<bond_list[iDX].size(); ++k){
      nbonds_k  = bond_list[iDX][k].nbonds;
      nbonds   += nbonds_k;
      rb       += nbonds_k * bond_list[iDX][k].dist;
    }
    rb /= nbonds;
    double Vb = Ep_list[iDX] / nbonds;

    printf( "%20.10f  %20.10f  # Compound %20s : ave(rbond) Vbond   nbonds %d  ave(nbonds/nat_with_bonds) %20.10f\n",
	   rb, Vb, DX[iDX].name.c_str(), nbonds, nbonds/(1.0*nat_with_bonds_list[iDX]) );
  }


}




