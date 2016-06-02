
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include <cmath>


#include "utils.hpp"
#include "utils-math.hpp"
#include "constants.hpp"
#include "utils-matrix.hpp"
#include "utils-string.hpp"
#include "utils-vector.hpp"
#include "utils-errors.hpp"

//#include "compound.hpp"
#include "compoundfit.hpp"
#include "physconst.hpp"
#include "errors.hpp"

#include "lattice-simple.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;



// Thanks to Morten for pointing out this construct for e.g. clang++:
extern "C" {
#include <spglib/spglib.h>
}




using namespace utils;
using namespace constants;
using namespace physconst;
using boost::format;


LatticeSimple::LatticeSimple(){
  nbasis   = 0;
  origin   = Vector3<double>(0.0);
  minpos   = Vector3<double>(0.0);
  avec     = Vector3<double>(0.0);
  bvec     = Vector3<double>(0.0);
  cvec     = Vector3<double>(0.0);
  pbc      = Vector3<bool>(true);
  csystem  = "unknown";
  csystem_sub = -1;
  csymaxis = "";
  pointgroup = "1";
  spacegroup = "1";
  spacegroup_number = 0;
}




void latsymm(Vector<CompoundStructureFit> & cmplist){

  int ic, nc=cmplist.size();
  int nbasis;
  int i, j;
  int imin;
  double sum,summin;
  MatrixSq3<double> A, Ainv;
  Vector3<double> b, x;
  Vector3<double> u1_vec, u2_vec, u3_vec;
  std::string name;
  Vector< Vector3<double> > pos, ipos;
  CompoundStructureFit cmp;


  for (ic=0; ic<nc; ++ic){
    // Local copy:
    cmp = cmplist[ic];

    if (cmp.csystem=="any") continue;

    if (!cmp.pbc[0] || !cmp.pbc[1] || !cmp.pbc[2]) continue;

    cmp.csystem="unknown";
    cmp.csystem_sub = -1;
    cmp.csymaxis = "";
    cmp.pointgroup = "1";
    cmp.spacegroup = "1";
    cmp.spacegroup_number = 0;





    name = cmp.name;
    nbasis = cmp.nbasis;
    pos  = cmp.basis_vecs;
    ipos = cmp.basis_vecs;


    std::string dumpfile("symana-" + name + ".out");
    std::ofstream fout;
    fout.open(dumpfile.c_str());
    fout << "Symmetry analysis log for compound " << name << std::endl;


    //cout << "nbasis " << nbasis << endl;
    //cout << "types " << types << endl;


    // Find atom with smallest (most negative) coordinates:
    for (i=0; i<nbasis; ++i){
      sum = pos[i][0] + pos[i][1] + pos[i][2];
    
      if (i==0 || (i>0 && sum<summin)){
	imin   = i;
	summin = sum;
      }
    }

    // Find internal positions:
    A.elem(0,0)=cmp.u1_vec[0];
    A.elem(1,0)=cmp.u1_vec[1];
    A.elem(2,0)=cmp.u1_vec[2];
    A.elem(0,1)=cmp.u2_vec[0];
    A.elem(1,1)=cmp.u2_vec[1];
    A.elem(2,1)=cmp.u2_vec[2];
    A.elem(0,2)=cmp.u3_vec[0];
    A.elem(1,2)=cmp.u3_vec[1];
    A.elem(2,2)=cmp.u3_vec[2];
    for (i=0; i<nbasis; ++i){
      b = pos[i];
      A.solve(b,x,Ainv);
      ipos[i] = x;
      for (j=0; j<3; ++j){
	while (ipos[i][j]< 0.0) ipos[i][j]++;
	while (ipos[i][j]>=1.0) ipos[i][j]--;
      }
    }


    LatticeSimple lattice;

    lattice.nbasis   = nbasis;
    lattice.minpos   = pos[imin];
    lattice.origin   = Vector3<double>(0);
    lattice.avec  = cmp.u1_vec;
    lattice.bvec  = cmp.u2_vec;
    lattice.cvec  = cmp.u3_vec;
    lattice.pbc   = cmp.pbc;
    lattice.types = cmp.basis_types;
    lattice.pos   = pos;
    lattice.ipos  = ipos;
    /*
    cout << "lattice nbasis " << lattice.nbasis << endl;
    cout << "lattice origin " << lattice.origin << endl;
    cout << "lattice avec " << lattice.avec << endl;
    cout << "lattice bvec " << lattice.bvec << endl;
    cout << "lattice cvec " << lattice.cvec << endl;
    cout << "lattice pos " << lattice.pos << endl;
    cout << "lattice pbc " << lattice.pbc << endl;
    cout << "lattice types " << lattice.types << endl;
    */




    // ##########################################################################
    // ##########################################################################
    // ##########################################################################
    // spglib stuff
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################

    double lat[3][3];
    lat[0][0] = lattice.avec[0];
    lat[1][0] = lattice.avec[1];
    lat[2][0] = lattice.avec[2];
    lat[0][1] = lattice.bvec[0];
    lat[1][1] = lattice.bvec[1];
    lat[2][1] = lattice.bvec[2];
    lat[0][2] = lattice.cvec[0];
    lat[1][2] = lattice.cvec[1];
    lat[2][2] = lattice.cvec[2];

    int *stypes = new int [nbasis];
    
    typedef double (P3)[3];
    double (*position)[3]; // 'position' is a pointer to an object, object is of type '3-element-array of double'
    //position = new double [nbasis][3];
    position = new P3 [nbasis];

    for (i=0; i<nbasis; ++i){
      stypes[i] = lattice.types[i];
      position[i][0] = lattice.ipos[i][0];
      position[i][1] = lattice.ipos[i][1];
      position[i][2] = lattice.ipos[i][2];
    }

    lattice.csystem_sub = 1;
    cmp.csystem_sub = 1;

    int num_spg, num_atom = nbasis;
    char symbol[21];
    double tol = 1.0e-5;//sqrt( eps_d() );

    SpglibDataset *p_spg_dataset;
    p_spg_dataset = spg_get_dataset(lat, position, stypes, nbasis, 1e-5);
    std::cout << "spglib: spacegroup number   : " << p_spg_dataset->spacegroup_number << std::endl;
    fout << "spglib: spacegroup number   : " << p_spg_dataset->spacegroup_number << std::endl;

    lattice.spacegroup_number = p_spg_dataset->spacegroup_number;
    cmp.spacegroup_number = p_spg_dataset->spacegroup_number;

    // http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
    int sp=p_spg_dataset->spacegroup_number;
    if (sp<=0 || sp>230) aborterror("Spacegroup could not be found for compound " + name );

    std::string csystem;
    if      (sp>=1   && sp<=2  ) csystem="triclinic";
    else if (sp>=3   && sp<=15 ) csystem="monoclinic";
    else if (sp>=16  && sp<=74 ) csystem="orthorombic";
    else if (sp>=75  && sp<=142) csystem="tetragonal";
    else if (sp>=143 && sp<=167) csystem="trigonal";
    else if (sp>=168 && sp<=194) csystem="hexagonal";
    else if (sp>=195 && sp<=230) csystem="cubic";

    lattice.csystem = csystem;
    cmp.csystem = csystem;

    std::cout << "Compound " << name << " has crystal system " << csystem << std::endl;
    fout << "Compound " << name << " has crystal system " << csystem << std::endl;

    std::cout << "spglib: hall number         : " << p_spg_dataset->hall_number << std::endl;
    std::cout << "spglib: international symbol: " << p_spg_dataset->international_symbol << std::endl;
    std::cout << "spglib: hall symbol         : " << p_spg_dataset->hall_symbol << std::endl;
    std::cout << "spglib: choice              : " << p_spg_dataset->choice << std::endl;
    fout << "spglib: hall number         : " << p_spg_dataset->hall_number << std::endl;
    fout << "spglib: international symbol: " << p_spg_dataset->international_symbol << std::endl;
    fout << "spglib: hall symbol         : " << p_spg_dataset->hall_symbol << std::endl;
    fout << "spglib: choice              : " << p_spg_dataset->choice << std::endl;
    
    SpglibSpacegroupType spg_type = spg_get_spacegroup_type(p_spg_dataset->hall_number);
    std::cout << "spglib: international symbol (full) : " << spg_type.international_full << std::endl;
    std::cout << "spglib: international symbol (short): " << spg_type.international_short << std::endl;
    std::cout << "spglib: schoenflies                 : " << spg_type.schoenflies << std::endl;
    fout << "spglib: international symbol (full) : " << spg_type.international_full << std::endl;
    fout << "spglib: international symbol (short): " << spg_type.international_short << std::endl;
    fout << "spglib: schoenflies                 : " << spg_type.schoenflies << std::endl;

    lattice.spacegroup = spg_type.international_full;
    cmp.spacegroup = spg_type.international_full;

    
    // triclinic
    if (lattice.spacegroup_number>=  1 && lattice.spacegroup_number<=  1) lattice.pointgroup = "1";
    if (lattice.spacegroup_number>=  2 && lattice.spacegroup_number<=  2) lattice.pointgroup = "bar(3)";
    // monoclinic
    if (lattice.spacegroup_number>=  3 && lattice.spacegroup_number<=  5) lattice.pointgroup = "2";
    if (lattice.spacegroup_number>=  6 && lattice.spacegroup_number<=  9) lattice.pointgroup = "m";
    if (lattice.spacegroup_number>= 10 && lattice.spacegroup_number<= 15) lattice.pointgroup = "2/m";
    // orthorombic
    if (lattice.spacegroup_number>= 16 && lattice.spacegroup_number<= 24) lattice.pointgroup = "222";
    if (lattice.spacegroup_number>= 25 && lattice.spacegroup_number<= 46) lattice.pointgroup = "mm2";
    if (lattice.spacegroup_number>= 47 && lattice.spacegroup_number<= 74) lattice.pointgroup = "2/m2/m2/m";
    // tetragonal
    if (lattice.spacegroup_number>= 75 && lattice.spacegroup_number<= 80) lattice.pointgroup = "4";
    if (lattice.spacegroup_number>= 81 && lattice.spacegroup_number<= 82) lattice.pointgroup = "bar(4)";
    if (lattice.spacegroup_number>= 83 && lattice.spacegroup_number<= 88) lattice.pointgroup = "4/m";
    if (lattice.spacegroup_number>= 99 && lattice.spacegroup_number<=110) lattice.pointgroup = "4mm";
    if (lattice.spacegroup_number>=111 && lattice.spacegroup_number<=122) lattice.pointgroup = "bar(4)2m";
    if (lattice.spacegroup_number>= 89 && lattice.spacegroup_number<= 98) lattice.pointgroup = "422";
    if (lattice.spacegroup_number>=123 && lattice.spacegroup_number<=142) lattice.pointgroup = "4/m2/m2/m";
    // trigonal
    if (lattice.spacegroup_number>=143 && lattice.spacegroup_number<=146) lattice.pointgroup = "3";
    if (lattice.spacegroup_number>=147 && lattice.spacegroup_number<=148) lattice.pointgroup = "bar(3)";
    if (lattice.spacegroup_number>=149 && lattice.spacegroup_number<=155) lattice.pointgroup = "32";
    if (lattice.spacegroup_number>=162 && lattice.spacegroup_number<=167) lattice.pointgroup = "bar(3)2/m";
    if (lattice.spacegroup_number>=156 && lattice.spacegroup_number<=161) lattice.pointgroup = "3m";
    // hexagonal
    if (lattice.spacegroup_number>=168 && lattice.spacegroup_number<=173) lattice.pointgroup = "6";
    if (lattice.spacegroup_number>=174 && lattice.spacegroup_number<=174) lattice.pointgroup = "bar(6)";
    if (lattice.spacegroup_number>=175 && lattice.spacegroup_number<=176) lattice.pointgroup = "6/m";
    if (lattice.spacegroup_number>=177 && lattice.spacegroup_number<=182) lattice.pointgroup = "622";
    if (lattice.spacegroup_number>=183 && lattice.spacegroup_number<=186) lattice.pointgroup = "6mm";
    if (lattice.spacegroup_number>=187 && lattice.spacegroup_number<=190) lattice.pointgroup = "bar(6)m2";
    if (lattice.spacegroup_number>=191 && lattice.spacegroup_number<=194) lattice.pointgroup = "6/m2/m2/m";
    // cubic
    if (lattice.spacegroup_number>=195 && lattice.spacegroup_number<=199) lattice.pointgroup = "23";
    if (lattice.spacegroup_number>=200 && lattice.spacegroup_number<=206) lattice.pointgroup = "2/mbar(3)";
    if (lattice.spacegroup_number>=207 && lattice.spacegroup_number<=214) lattice.pointgroup = "432";
    if (lattice.spacegroup_number>=215 && lattice.spacegroup_number<=220) lattice.pointgroup = "bar(4)3m";
    if (lattice.spacegroup_number>=221 && lattice.spacegroup_number<=230) lattice.pointgroup = "4/mbar(3)2/m";

    cmp.pointgroup = lattice.pointgroup;



    char c = p_spg_dataset->choice[0];
    bool ax_is_small = fp_is_small_tol(lattice.avec[0], tol);
    bool ay_is_small = fp_is_small_tol(lattice.avec[1], tol);
    bool az_is_small = fp_is_small_tol(lattice.avec[2], tol);
    bool bx_is_small = fp_is_small_tol(lattice.bvec[0], tol);
    bool by_is_small = fp_is_small_tol(lattice.bvec[1], tol);
    bool bz_is_small = fp_is_small_tol(lattice.bvec[2], tol);
    bool cx_is_small = fp_is_small_tol(lattice.cvec[0], tol);
    bool cy_is_small = fp_is_small_tol(lattice.cvec[1], tol);
    bool cz_is_small = fp_is_small_tol(lattice.cvec[2], tol);
    if      (c=='a'){
      if      (ax_is_small && ay_is_small) cmp.csymaxis = "z";
      else if (ax_is_small && az_is_small) cmp.csymaxis = "y";
      else if (ay_is_small && az_is_small) cmp.csymaxis = "x";
      else cmp.csymaxis = "a";
    }
    else if (c=='b'){
      if      (bx_is_small && by_is_small) cmp.csymaxis = "z";
      else if (bx_is_small && bz_is_small) cmp.csymaxis = "y";
      else if (by_is_small && bz_is_small) cmp.csymaxis = "x";
      else cmp.csymaxis = "b";
    }
    else if (c=='c'){
      if      (cx_is_small && cy_is_small) cmp.csymaxis = "z";
      else if (cx_is_small && cz_is_small) cmp.csymaxis = "y";
      else if (cy_is_small && cz_is_small) cmp.csymaxis = "x";
    }


    if (cmp.csystem=="monoclinic"){
      if ( ! (cmp.csymaxis == "y" || cmp.csymaxis == "z") ){
	aborterror("ERROR: Compound " + name +
		   " should have main symmetry axis either 'y' or 'z'." +
		   " Now it is " + cmp.csymaxis);
      }
    }


    spg_free_dataset(p_spg_dataset);
    delete [] stypes;
    delete [] position;

    // ##########################################################################
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################
    // ##########################################################################


    std::cout << "crystal system " << cmp.csystem
	 << " subclass " << cmp.csystem_sub
	 << " point group " << cmp.pointgroup
	 << " symmetry axis " << cmp.csymaxis
	 << std::endl;
    fout << "crystal system " << cmp.csystem
	 << " subclass " << cmp.csystem_sub
	 << " point group " << cmp.pointgroup
	 << " symmetry axis " << cmp.csymaxis
	 << std::endl;


    fout.close();
    fout.clear();


    // Copy back:
    cmplist[ic] = cmp;

  }




}




