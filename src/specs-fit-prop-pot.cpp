

#include <limits>

#include <cmath>

#include "mdsettings.hpp"
#include "specs-fit-prop-pot.hpp"

using namespace std;




SpecsFitProp::SpecsFitProp(){
  double small = sqrt( std::numeric_limits<double>::epsilon() );

  nitermax = 200;
  nitermin = 2;

  functolabs =  small;
  functolrel = -small;
  gradtolabs =  small;
  steptolabs =  small;
  steptolrel = -small;

  fitmet = "LM";
  seed = 12345;
  dogleg_radius = 1.0;
  dogleg_minradius = 0.1;
  simplex_delta = 0.1;
  barrier_scale = 1.0;

  debug_fit_level0 = false;
  debug_fit_level1 = false;
  debug_fit_level2 = false;
  debug_fit_level3 = false;
  debug_fit_level4 = false;
  report_iter = false;
  report_warn = false;
  report_error = false;
  report_conv = false;
  lattol = small;
  BM_fmin = -1e-3;
  BM_fmax =  1e-3;
  BM_Nf   = 10;
  C_fmin = -1e-3;
  C_fmax =  1e-3;
  C_Nf   = 10;

  mds_specs_common = MDSettings_common();
  mds_specs        = MDSettings();
  mds_specs_ref    = MDSettings();

  mds_specs_ref.Tstart = 0.01;
  mds_specs_ref.btc_tau = 20.0;
  mds_specs_ref.btc_T0 = 0.01;
  mds_specs_ref.is_ref_comp = true;
  mds_specs_ref.heating_allowed = true;
  mds_specs_ref.fixed_geometry = false;
  mds_specs_ref.quench_always = false;


}


SpecsFitPot::SpecsFitPot(){
  double small = sqrt( std::numeric_limits<double>::epsilon() );

  seed = 12345;
  nitermax = 200;
  nitermin = 2;

  functolabs =  small;
  functolrel = -small;
  gradtolabs =  small;
  steptolabs =  small;
  steptolrel = -small;

  fitmet = "LM";
  dogleg_radius = 1.0;
  dogleg_minradius = 0.1;
  simplex_delta = 0.1;

  debug_fit_level0 = false;
  debug_fit_level1 = false;
  debug_fit_level2 = false;
  debug_fit_level3 = false;
  debug_fit_level4 = false;
  report_iter = true;
  report_warn = true;
  report_error = true;
  report_conv = true;
  barrier_scale = 1.0;
}



