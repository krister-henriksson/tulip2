

#include <limits>

#include <cmath>

#include "mdsettings.hpp"
#include "specs-fit-prop-pot.hpp"




SpecsFitProp::SpecsFitProp(){
  double small = sqrt( std::numeric_limits<double>::epsilon() );

  nitermax = 200;
  nitermin = 2;
  niterrestart = 20;

  functolabs =  small;
  functolrel = -small;
  gradtolabs =  small;
  steptolabs =  small;
  steptolrel = -small;

  fitmet = "LM";
  seed = 12345;
  dogleg_radius = 0.5;
  dogleg_minradius = 1e-5;
  simplex_delta = 0.2;
  simann_delta_rel = 0.2;
  moldyn_min_dx = 1.0e-10;
  moldyn_max_dx = 0.1;
  barrier_scale = 0.0;
  use_barrier_rescaling = false;
  use_data_scales = true;

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
  BM_rel_sys = false;
  BM_fmin = -0.01;
  BM_fmax =  0.01;
  BM_Nf   = 10;
  BM_ef   = 0.01;
  C_rel_sys = false;
  C_fmin = -0.01;
  C_fmax =  0.01;
  C_Nf   = 10;
  C_ef   = 0.01;

  mds_specs_common = MDSettings_common();

  mds_specs        = MDSettings();
  mds_specs_ref    = MDSettings();

  mds_specs_ref.is_ref_comp = true;
  mds_specs_ref.ext_relax = false;
  mds_specs_ref.quench_always = false;


}


SpecsFitPot::SpecsFitPot(){
  double small = sqrt( std::numeric_limits<double>::epsilon() );

  seed = 12345;
  nitermax = 200;
  nitermin = 2;
  niterrestart = 20;

  functolabs =  small;
  functolrel = -small;
  gradtolabs =  small;
  steptolabs =  small;
  steptolrel = -small;

  fitmet = "DL";
  dogleg_radius = 0.5;
  dogleg_minradius = 1e-5;
  simplex_delta = 0.2;
  simann_delta_rel = 0.2;
  moldyn_min_dx = 1.0e-10;
  moldyn_max_dx = 0.1;
  barrier_scale = 0.0;
  use_barrier_rescaling = false;
  use_data_scales = true;

  debug_fit_level0 = false;
  debug_fit_level1 = false;
  debug_fit_level2 = false;
  debug_fit_level3 = false;
  debug_fit_level4 = false;
  report_iter = true;
  report_warn = true;
  report_error = true;
  report_conv = true;
}



