

#include "compound.hpp"
#include "mdsettings.hpp"



// Set in main routine:
MDSettings_common::MDSettings_common(){
  debug_forces = false;
  debug_pressure = false;

  report_step  = false;

  use_def_dump_xyz_fmt = true;
  def_dump_xyz_fmt = "extxyz";

}




MDSettings::MDSettings(){
  // Default values:
  skint   = 1.0;
  seed    = 12345;
  ndump   = 10;
  tstart  = 0;
  tend    = 2000.0;

  Tstart  = 0.1;
  dt      = 3.0;
  max_dt  = 5.0;
  max_dE  = 1.0;
  max_dr  = 0.1;

  btc_tau = 20.0;
  btc_T0  = 0.0;

  bpc_tau   = 20.0;
  bpc_P0    = 0.0;
  bpc_scale = 80.0;

  use_quench   = false;
  use_Tcontrol = true;
  use_Pcontrol = true;

  quench_tstart = -1.0;
  quench_rate   = 0.0;

  use_error_T_gt = true;
  use_error_dt_lt = true;
  use_error_boxlen_gt = true;
  
  error_T_gt = 1.0e4;
  error_dt_lt = 0.01;
  error_boxlen_gt = 1.0e4;

  is_ref_comp = true;
  heating_allowed = true;
  fixed_geometry = false;
  quench_always = false;
}

