
#ifndef LATSYMM_HPP
#define LATSYMM_HPP




#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"

//#include "compound.hpp"
#include "compoundfit.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;



using namespace utils;



class LatticeSimple {
public:
  int nbasis;
  Vector3<double> origin;
  Vector3<double> minpos;
  Vector3<double> avec;
  Vector3<double> bvec;
  Vector3<double> cvec;
  Vector3<bool> pbc;
  Vector<int>              types;
  Vector< Vector3<double> > pos;
  Vector< Vector3<double> > ipos;
  std::string csystem;
  int    csystem_sub;
  std::string csymaxis;
  std::string pointgroup;
  std::string spacegroup;
  int spacegroup_number;

  LatticeSimple();

} ;




void latsymm(Vector<CompoundStructureFit> & cmps);





#endif

