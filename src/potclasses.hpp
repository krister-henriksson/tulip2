

#ifndef POTCLASSES_HPP
#define POTCLASSES_HPP


#include <string>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "param.hpp"


using namespace utils;





enum PotentialClasses { Reppot, EAMpot, ABOPpot };




// ###############################################################
// Tabulated potentials
// ###############################################################


class Potential_Reppot {
public:
  std::string elemname1;
  std::string elemname2;

  int Nr_rep;
  Vector<double> r_rep;
  Vector<double> V_rep;
  Vector<double> d2_V_rep;
  // double bfermi;
  // double rfermi;

  Potential_Reppot();
  double rcut(void);
};


class Potential_EAM {
public:
  std::string elemname1;
  std::string elemname2;

  int Nr;
  int Nrho_s;
  int Nrho_p;
  int Nrho_d;
  Vector<double> r;
  Vector<double> V2;
  Vector<double> F_s;
  Vector<double> F_p;
  Vector<double> F_d;
  Vector<double> F_s_rho_s;
  Vector<double> F_p_rho_p;
  Vector<double> F_d_rho_d;
  Vector<double> rho_s;
  Vector<double> rho_p;
  Vector<double> rho_d;
  Vector<double> d2_V2;
  Vector<double> d2_F_s;
  Vector<double> d2_F_p;
  Vector<double> d2_F_d;
  Vector<double> d2_rho_s;
  Vector<double> d2_rho_p;
  Vector<double> d2_rho_d;
  double dr;
  double drho_s;
  double drho_p;
  double drho_d;
  double mrcut;

  Potential_EAM();
  double rcut(void);
} ;



// ###############################################################
// Parametrized potentials
// ###############################################################



class Potential_ABOP {
public:
  std::string elemname1;
  std::string elemname2;

  Vector<double> parval;
  Vector<std::string> parname;

  // for fittable versions:
  Vector<parametertype> partype;
  Vector<double> parmin;
  Vector<double> parmax;

  int maxindex;

  std::string rcut_fun;


  Potential_ABOP();
  /*
  ~Potential_ABOP();
  Potential_ABOP(const Potential_ABOP & sv);
  Potential_ABOP & operator=(const Potential_ABOP & sv);
  */

  int         parname2idx(std::string name);
  std::string paridx2name(int idx);
  double    & parname2val(std::string name);
  double      rcut(void);

  // for fittable versions:
  void init_lims(void);

};





#endif





