

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
//#include "mdsettings.hpp"
//#include "mdsystem.hpp"
//#include "param-pot.hpp"


#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;

using namespace std;
using namespace utils;


//class MDSystem;
//class ParamPot;





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
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;
  bool frc;

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
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;
  Vector< Vector3<double> > frc;

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
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;
  bool frc;

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
  bool B;
  bool Bp;
  Matrix<bool> C;
  bool Fmax;
  bool Pmax;
  bool displmax;
  bool frc;

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
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;
  Vector< Vector3<double> > frc;

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
  double B;
  double Bp;
  Matrix<double> C;
  double Fmax;
  double Pmax;
  double displmax;
  Vector< Vector3<double> > frc;

  CompoundPropertiesWeights();
} ;



class CompoundStructure {
public:
  string filename;
  string filename_frc;
  string name;
  string crystalname;
  int nelem;
  Vector<string> elemnames;
  Vector3<bool> pbc;
  string csystem;
  int    csystem_sub;
  string csymaxis;
  string pointgroup;
  string spacegroup;
  int spacegroup_number;

  double scalefactor;
  double lpa, lpb, lpc;

  bool use_readin_structure;
  bool use_int;
  bool use_origin_spec;

  Vector3<double> origin;

  Vector3<double> u1_vec;
  Vector3<double> u2_vec;
  Vector3<double> u3_vec;

  int nbasis;
  Vector<int>              basis_types;
  Vector<string>           basis_elems;
  Vector< Vector3<double> > basis_vecs;

  // For construction of compound:
  int Ndesired[3];
  bool Neven_desired[3];
  bool Nodd_desired[3];

  bool Ecoh_delta_refcomp;


  CompoundStructure();
  //CompoundStructure(const CompoundStructure & sv);
  //CompoundStructure & operator=(const CompoundStructure & sv);


  // 1 ...
  void create_from_model(Elements & el,
			 string name, string elem1, string elem2,
			 double ai, double bi, double ci);
  // ... or 2
  void read_structure(Elements & el);


  // 1 and 2: Finalize basis vectors, etc:
  void finalize(const double ai, const double bi, const double ci);

  void origin_from_model(int & N1, int & N2, int & N3);

  //bool matrix_is_a_symm_op(const Matrix<double> & R);
  //void check_crystal_symm(void);

} ;





#endif

