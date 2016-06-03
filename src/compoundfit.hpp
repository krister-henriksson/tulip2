


#ifndef COMPOUND_FIT_HPP
#define COMPOUND_FIT_HPP



#include <iostream>
#include <string>

#include <cstdio>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"

#include "elem-iacs.hpp"
#include "specs-fit-prop-pot.hpp"
//#include "mdsettings.hpp"

#include "compound.hpp"


#include "mdsystem.hpp"
#include "param-pot.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;




using namespace utils;





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
  bool Ecoh_delta;
  bool Emix;
  bool Eform;
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;
  bool frc;

  bool bondlen;
  bool bondangle;

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
  double Ecoh_delta;
  double Emix;
  double Eform;
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;
  Vector< Vector3<double> > frc;

  Vector< double > bondlen;
  Vector< int >    bondlen_a1;
  Vector< int >    bondlen_a2;
  Vector< double > bondangle;
  Vector< int >    bondangle_a1;
  Vector< int >    bondangle_a2;
  Vector< int >    bondangle_a3;

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
  bool Ecoh_delta;
  bool Emix;
  bool Eform;
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;
  bool frc;

  Vector<bool> bondlen;
  Vector<bool> bondangle;
  
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
  bool Ecoh_delta;
  bool Emix;
  bool Eform;
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;
  bool frc;

  Vector<bool> bondlen;
  Vector<bool> bondangle;

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
  double Ecoh_delta;
  double Emix;
  double Eform;
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;
  Vector< Vector3<double> > frc;

  Vector< double > bondlen;
  Vector< double > bondangle;

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
  double Ecoh_delta;
  double Emix;
  double Eform;
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;
  Vector< Vector3<double> > frc;

  Vector< double > bondlen;
  Vector< double > bondangle;

  CompoundPropertiesWeights();
} ;





class CompoundStructureFit
  :
  public CompoundStructure
{
public:
  CompoundPropertiesUse              prop_use;
  CompoundProperties                 prop_readin;
  CompoundProperties                 prop_pred;

  MDSettings                         mds_specs;

  CompoundPropertiesUseUncertainties use_u;
  CompoundPropertiesUseWeights       use_w;
  CompoundPropertiesUncertainties    prop_u;
  CompoundPropertiesWeights          prop_w;


  CompoundStructureFit();
  // CompoundStructureFit(const CompoundStructureFit & sv);
  // CompoundStructureFit & operator=(const CompoundStructureFit & sv);

  void read_forces(void);

  void check_and_fix_uses();
  int NData();

  void getprop(ParamPot                  & param);
  void get_B_Bp(MDSystem                 & mds,
		ParamPot                 & param,
		Vector< Vector3<double> > & pos_bak,
		MatrixSq3<double>           & boxdir_bak,
		Vector3<double>           & boxlen_bak,
		double                   & E0,
		double                   & V0);
  void get_Cij(MDSystem                  & mds,
	       ParamPot                  & param,
	       Vector< Vector3<double> >  & pos_bak,
	       MatrixSq3<double>            & boxdir_bak,
	       Vector3<double>            & boxlen_bak,
	       double                    & E0,
	       double                    & V0);


  void check_and_fix_Cij(void);
  void get_Cuse(Matrix<bool> & Cuse);
  void get_Cresolved(Matrix<double> & Clincomb,
		     Matrix<double> & Cfull);
					 


} ;
  



#endif

