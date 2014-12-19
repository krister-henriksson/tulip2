

#ifndef POTINFO_HPP
#define POTINFO_HPP


#include <iostream>
#include <string>

#include <cstdio>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"

#include "elem-iacs.hpp"
#include "specs-fit-prop-pot.hpp"
#include "potclasses.hpp"
#include "param.hpp"

#include "omp-basics.hpp"

using namespace std;
using namespace utils;



class PotentialInformation {
private:
  Matrix<string>  mbasepot;
  Matrix<int>     mbasepot_vecidx;
  Matrix<int>     mreppot_vecidx;
  Vector<bool>    muse_reppot;
  Matrix3<double> mabop_omega;

public:
  Elements       elem;
  Interactions   iacs;
  Vector<double> Ecoh_ref;

  Vector<Potential_ABOP>   pot_ABOP;
  Vector<Potential_EAM>    pot_EAM;
  Vector<Potential_Reppot> pot_Reppot;

  Matrix3<double> abop_alpha;
  // Matrix3<bool>   abop_omega_is_indep;
  // Note: abop_omega is set/get by methods.
  Matrix3<bool>   use_abop_alpha;
  Matrix3<bool>   use_abop_omega;

  Matrix<double>  abop_2mu;
  Matrix< bool>   use_abop_2mu;



  PotentialInformation();
  PotentialInformation(string filename_info);

  // -------------------------------------------------------------------------------
  // Called from constructor:
  // Read settings about potentials:
  void read_info(string file);

  // calls:
  //  void init(Elements & el, Interactions & ia);

  string basepot(string s1, string s2);
  string basepot(int i1, int i2);
  int basepot_vecidx(string s1, string s2);
  int basepot_vecidx(int i1, int i2);

  bool & use_reppot(string s1, string s2);
  bool & use_reppot(int i1, int i2);

  void reppot_finalize(void);
  int reppot_vecidx(string s1, string s2);
  int reppot_vecidx(int i1, int i2);

  void read_eampot(void);
  void read_reppot(void);

  // -------------------------------------------------------------------------------

  // General utilities for later use:
  double get_abop_omega(string s1, string s2, string s3);
  double get_abop_omega(int i1, int i2, int i3);
  void   set_abop_omega(string s1, string s2, string s3, double val);
  void   set_abop_omega(int i1, int i2, int i3, double val);

  double get_rcut_max(void);
  double get_rcut_max(Vector<string> elemnames);



} ;


// ###############################################################################
// ###############################################################################


class PotentialInformationFit
  : public PotentialInformation
{
private:
  Vector<bool>    mfit;        /* True when the pot. is to be fitted. */

public:
  SpecsFitProp specs_prop;
  SpecsFitPot  specs_pot;
  OMP_Info     omp_info;

  Matrix3<double>        abop_alpha_parmin;
  Matrix3<double>        abop_alpha_parmax;
  Matrix3<parametertype> abop_alpha_partype;

  Matrix3<double>        abop_omega_parmin;
  Matrix3<double>        abop_omega_parmax;
  Matrix3<parametertype> abop_omega_partype;

  Matrix<double>         abop_2mu_parmin;
  Matrix<double>         abop_2mu_parmax;
  Matrix<parametertype>  abop_2mu_partype;



  PotentialInformationFit();
  PotentialInformationFit(string filename_info,
			  string filename_specs);

  bool & is_fittable(string s1, string s2);
  bool & is_fittable(int i1, int i2);

  void read_info_fit(string filename_info);

  // Read settings about fitting structure properties and potential:
  void read_specs(string filename_specs);


  void limcheck(const string & pot,
		const string & parname,
		const string & elems,
		const parametertype & partype,
		const double & parmin,
		const double & parmax,
		const double & parval
		);

} ;



#endif





