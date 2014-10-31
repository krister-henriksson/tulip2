

#include "compound.hpp"
#include "mdsettings.hpp"



// Set in main routine:
MDSettings_common::MDSettings_common(){
  debug_forces = false;
  debug_pressure = false;

  report_step  = false;
  quick_mode   = false;
}




MDSettings::MDSettings(){
  // Default values:
  skint   = 1.0;
  seed    = 12345;
  ndump   = 10;
  tstart  = 0;
  tend    = 2000.0;
  Tstart  = 0.01;
  dt      = 1.0;
  max_dt  = 1.0;
  max_dE  = 30.0;
  max_dr  = 0.2;

  use_quench   = false;
  use_Tcontrol = false;
  use_Pcontrol = false;

  btc_tau = 20.0;
  btc_T0  = 0.0;
  quench_tstart = -1.0;
  quench_rate   = 0.0;
  bpc_tau   = 100.0;
  bpc_P0    = 0.0;
  bpc_scale = 100.0;


  is_ref_comp = true;
  heating_allowed = true;
  fixed_geometry = false;
  quench_always = false;
}
