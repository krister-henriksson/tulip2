

#ifndef ELEM_IACS_HPP
#define ELEM_IACS_HPP


#include <string>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"

using namespace std;
using namespace utils;



class Elements {
private:
  Vector<string> melemnames;
  Vector<int>   matomtype;
  Vector<double> mmass;
  Vector<string> mreflat;
  Vector<double> mreflat_a;
  Vector<double> mreflat_b;
  Vector<double> mreflat_c;
  Vector<double> mreflat_bpa;
  Vector<double> mreflat_cpa;

public:

  void add_elem(string name);
  int nelem(void);
  int name2idx(string name);
  string idx2name(int typei);

  int   & atomtype(string s);
  double & mass(string s);
  string & reflat(string s);
  double & reflat_a(string s);
  double & reflat_b(string s);
  double & reflat_c(string s);
  double & reflat_bpa(string s);
  double & reflat_cpa(string s);

  int   & atomtype(int typei);
  double & mass(int typei);
  string & reflat(int typei);
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
  Matrix<string> miacnames;

public:
  void init(Elements el);
  string & name(string n1, string n2);
  string & name(int typei1, int typei2);

} ;



#endif





