

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

// ###############################################################


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
// Helper classes for parametrized potentials
// ###############################################################


class CutoffScreeningPair {
public:
  bool tersoff;
  bool perriot_cut;
  bool perriot_scr;

  double R;
  double D;
  double pn;
  double pm;
  double prcut;
  double prmin;
  double prmax;

  double rcut;
  
  CutoffScreeningPair();
} ;


// ###############################################################

class CutoffScreening {
public:
  std::string mode;
  std::string name;
  bool use;

  Vector<double>      parval;
  Vector<std::string> parname;

  // for fittable versions:
  Vector<parametertype> partype;
  Vector<double> parmin;
  Vector<double> parmax;


  CutoffScreening();

  void set_name_mode(std::string na, std::string mo);
  void   set_parval(std::string name, double pv);
  double get_parval(std::string name);

  double      rcut(void);

  // for fittable versions:
  void init_lims(void);

  void   set_par_extr(std::string name, double min, std::string minmax);
  double get_par_extr(std::string name, std::string minmax);

  int           npar();
  void          set_par_types();
  void          set_par_type(std::string name);
  parametertype get_par_type(std::string name);

} ;





// ###############################################################
// Parametrized potentials
// ###############################################################

class Potential_ABOPPair {
public:
  std::string elemname1;
  std::string elemname2;

  double D0;
  double r0;
  double beta;
  double S;
  double p;
  double gamma;
  double c;
  double d;
  double h;
  double bfermi;
  double rfermi;

  Potential_ABOPPair();
} ;



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

  /*
  bool use_cutoff_only;
  bool use_screening;

  std::string rcut_fun;
  std::string rcut_scr;
  */
  CutoffScreening rcs;


  Potential_ABOP();
  /*
  ~Potential_ABOP();
  Potential_ABOP(const Potential_ABOP & sv);
  Potential_ABOP & operator=(const Potential_ABOP & sv);
  */

  void   set_parval(std::string name, double pv);
  double get_parval(std::string name);

  double      rcut(void);

  // for fittable versions:
  void init_lims(void);

  void   set_par_extr(std::string name, double min, std::string minmax);
  double get_par_extr(std::string name, std::string minmax);

  int           npar();
  void          set_par_types();
  void          set_par_type(std::string name);
  parametertype get_par_type(std::string name);

};





#endif





