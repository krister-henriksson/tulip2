
#ifndef LATSYMM_HPP
#define LATSYMM_HPP




#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"

#include "compound.hpp"
#include "compoundfit.hpp"



using namespace std;
using namespace utils;



class LatticeSimple {
public:
  int nbasis;
  Vector<double> origin;
  Vector<double> minpos;
  Vector<double> avec;
  Vector<double> bvec;
  Vector<double> cvec;
  Vector<bool> pbc;
  Vector<int>              types;
  Vector< Vector<double> > pos;
  Vector< Vector<double> > ipos;
  string csystem;
  int    csystem_sub;
  string csymaxis;
  string pointgroup;
  string spacegroup;
  int spacegroup_number;

  LatticeSimple();

  void get_ipos(void);
  void rotate(const Matrix<double> & R);
  void shift(const Vector<double> & svec);

  void dump_xyz(const string & filename);

  bool matrix_is_a_symm_op(const Matrix<double> & R);

} ;




void latsymm(Vector<CompoundStructureFit> & cmps);





#endif

