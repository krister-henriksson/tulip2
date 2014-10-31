

#ifndef COMPOUND_HPP
#define COMPOUND_HPP



#include <iostream>
#include <string>

#include <cstdio>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"

#include "elem-iacs.hpp"
#include "specs-fit-prop-pot.hpp"
#include "mdsettings.hpp"
//#include "mdsystem.hpp"
//#include "param-pot.hpp"

using namespace std;
using namespace utils;


class MDSystem;
class ParamPot;





class CompoundPropertiesUse {
public:
  bool a;
  bool b;
  bool c;
  bool bpa;
  bool cpa;
  bool r0;
  bool angle_ab;
  bool angle_ac;
  bool angle_bc;
  bool Vatom;
  bool Ecoh;
  bool Emix;
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;

  CompoundPropertiesUse();

  void check_and_fix();
} ;


class CompoundProperties {
public:
  double a;
  double b;
  double c;
  double bpa;
  double cpa;
  double r0;
  double angle_ab;
  double angle_ac;
  double angle_bc;
  double Vatom;
  double Ecoh;
  double Emix;
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;

  CompoundProperties();
} ;


class CompoundPropertiesUseUncertainties {
public:
  bool a;
  bool b;
  bool c;
  bool bpa;
  bool cpa;
  bool r0;
  bool angle_ab;
  bool angle_ac;
  bool angle_bc;
  bool Vatom;
  bool Ecoh;
  bool Emix;
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;

  CompoundPropertiesUseUncertainties();
} ;

class CompoundPropertiesUseWeights {
public:
  bool a;
  bool b;
  bool c;
  bool bpa;
  bool cpa;
  bool r0;
  bool angle_ab;
  bool angle_ac;
  bool angle_bc;
  bool Vatom;
  bool Ecoh;
  bool Emix;
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;

  CompoundPropertiesUseWeights();
} ;



class CompoundPropertiesUncertainties {
public:
  double a;
  double b;
  double c;
  double bpa;
  double cpa;
  double r0;
  double angle_ab;
  double angle_ac;
  double angle_bc;
  double Vatom;
  double Ecoh;
  double Emix;
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;

  CompoundPropertiesUncertainties();
} ;


class CompoundPropertiesWeights {
public:
  double a;
  double b;
  double c;
  double bpa;
  double cpa;
  double r0;
  double angle_ab;
  double angle_ac;
  double angle_bc;
  double Vatom;
  double Ecoh;
  double Emix;
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;

  CompoundPropertiesWeights();
} ;



class CompoundStructure {
public:
  string filename;
  string name;
  string crystalname;
  int nelem;
  Vector<string> elemnames;
  Vector<bool> pbc;
  string csystem;

  double scalefactor;

  bool use_int;

  Vector<double> u1_vec;
  Vector<double> u2_vec;
  Vector<double> u3_vec;

  int nbasis;
  Vector<string>           basis_elems;
  Vector< Vector<double> > basis_vecs;


  CompoundStructure();
  //CompoundStructure(const CompoundStructure & sv);
  //CompoundStructure & operator=(const CompoundStructure & sv);


  void create_from_model(string name,
			 string elem1,
			 string elem2
			 );
  void read_structure(void);
  void finalize(const double ai, const double bi, const double ci);

} ;



// #####################################################################
// #####################################################################


class CompoundStructureFit
  :
  public CompoundStructure
{
public:
  CompoundPropertiesUse prop_use;
  CompoundProperties    prop_readin;
  CompoundProperties    prop_pred;

  MDSettings            mds_specs;

  CompoundPropertiesUseUncertainties use_u;
  CompoundPropertiesUseWeights       use_w;

  CompoundPropertiesUncertainties prop_u;
  CompoundPropertiesWeights       prop_w;


  CompoundStructureFit();
  // CompoundStructureFit(const CompoundStructureFit & sv);
  // CompoundStructureFit & operator=(const CompoundStructureFit & sv);

  void check_and_fix_uses();
  int NData();

  void getprop(ParamPot & param);
  void get_B_Bp(MDSystem             & mds,
		ParamPot             & param,
		Vector< Vector<double> > & pos_bak,
		Matrix<double>           & boxdir_bak,
		Vector<double>           & boxlen_bak,
		double & E0,
		double & V0);
  void get_Cij(MDSystem             & mds,
	       ParamPot             & param,
	       Vector< Vector<double> > & pos_bak,
	       Matrix<double>           & boxdir_bak,
	       Vector<double>           & boxlen_bak,
	       double & E0,
	       double & V0);


} ;
  



// #####################################################################
// #####################################################################




class CompoundListFit {
public:
  Elements elem;
  int ncompounds;
  Vector<CompoundStructureFit> compounds;

  // Read structures:
  CompoundListFit(const Elements & el,
		  MDSettings & mds_specs,
		  string filename);
  int NData();
} ;




#endif

