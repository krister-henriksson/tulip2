

#ifndef ELEM_IACS_HPP
#define ELEM_IACS_HPP


#include <string>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;



using namespace utils;



class Elements {
private:
  Vector<std::string> melemnames;
  Vector<int>   matomtype;
  Vector<double> mmass;
  Vector<std::string> mreflat;
  Vector<double> mreflat_a;
  Vector<double> mreflat_b;
  Vector<double> mreflat_c;
  Vector<double> mreflat_bpa;
  Vector<double> mreflat_cpa;

public:

  void add_elem(std::string name);
  int nelem(void);
  int name2idx(std::string name);
  std::string idx2name(int typei);

  int   & atomtype(std::string s);
  double & mass(std::string s);
  std::string & reflat(std::string s);
  double & reflat_a(std::string s);
  double & reflat_b(std::string s);
  double & reflat_c(std::string s);
  double & reflat_bpa(std::string s);
  double & reflat_cpa(std::string s);

  int   & atomtype(int typei);
  double & mass(int typei);
  std::string & reflat(int typei);
  double & reflat_a(int typei);
  double & reflat_b(int typei);
  double & reflat_c(int typei);
  double & reflat_bpa(int typei);
  double & reflat_cpa(int typei);


  Vector<double> masses();

} ;

// ******************************************************************************
// ******************************************************************************

class Interactions {
private:
  Elements       mel;
  Matrix<std::string> miacnames;

public:
  void init(Elements el);
  std::string & name(std::string n1, std::string n2);
  std::string & name(int typei1, int typei2);

} ;



#endif





