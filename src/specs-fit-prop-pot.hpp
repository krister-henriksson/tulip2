


#ifndef SPECS_PROP_POT_FIT_HPP
#define SPECS_PROP_POT_FIT_HPP


#include <string>

#include "mdsettings.hpp"





class SpecsFitProp
{
public:
  int seed;

  std::string fitmet;

  int   nitermax;
  int   nitermin;
  int niterrestart;

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
  /* Simulated Annealing initial change in parameters (relative fraction): */
  double simann_delta_rel;
  /* MD: */
  double moldyn_min_dx;
  double moldyn_max_dx;
  /* Restricted optimization options: */
  double   barrier_scale;
  bool use_barrier_rescaling;

  bool use_data_scales;

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
  bool     BM_rel_sys;
  double   BM_fmin;
  double   BM_fmax;
  double   BM_Nf;
  double   BM_ef;

  bool     C_rel_sys;
  double   C_fmin;
  double   C_fmax;
  double   C_Nf;
  double   C_ef;
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

  std::string fitmet;
  int   nitermax;
  int   nitermin;
  int niterrestart;

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
  /* Simulated Annealing initial change in parameters (relative fraction): */
  double simann_delta_rel;
  /* MD: */
  double moldyn_min_dx;
  double moldyn_max_dx;
  /* Restricted optimization options: */
  double   barrier_scale;
  bool use_barrier_rescaling;

  bool use_data_scales;

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

