


#ifndef SPECS_PROP_POT_FIT_HPP
#define SPECS_PROP_POT_FIT_HPP


#include <string>

#include "mdsettings.hpp"


using namespace std;




class SpecsFitProp
{
public:
  int seed;

  string fitmet;

  int   nitermax;
  int   nitermin;

  double functolabs;
  double functolrel;
  double gradtolabs;
  double steptolabs;
  double steptolrel;

  /* Dog-Leg method options: */
  double   dogleg_radius;
  double   dogleg_minradius;
  /* Simplex method options: */
  double   simplex_delta;
  /* Restricted optimization options: */
  double   barrier_scale;

  // ***********************************
  bool debug_fit_level0;
  bool debug_fit_level1;
  bool debug_fit_level2;
  bool debug_fit_level3;
  bool debug_fit_level4;
  bool report_iter;
  bool report_warn;
  bool report_error;
  bool report_conv;
  // ***********************************
  double   lattol;
  // ***********************************
  double   BM_fmin;
  double   BM_fmax;
  double   BM_Nf;
  double   C_fmin;
  double   C_fmax;
  double   C_Nf;
  // ***********************************
  MDSettings_common mds_specs_common;
  // ***********************************
  MDSettings        mds_specs;
  MDSettings        mds_specs_ref; // for reference compounds

  SpecsFitProp();
} ;
  


class SpecsFitPot
{
public:
  int seed;

  string fitmet;
  int   nitermax;
  int   nitermin;

  double functolabs;
  double functolrel;
  double gradtolabs;
  double steptolabs;
  double steptolrel;

  /* Dog-Leg method options: */
  double   dogleg_radius;
  double   dogleg_minradius;
  /* Simplex method options: */
  double   simplex_delta;
  /* Restricted optimization options: */
  double   barrier_scale;

  bool debug_fit_level0;
  bool debug_fit_level1;
  bool debug_fit_level2;
  bool debug_fit_level3;
  bool debug_fit_level4;
  bool report_iter;
  bool report_warn;
  bool report_error;
  bool report_conv;


  SpecsFitPot();
} ;




#endif

