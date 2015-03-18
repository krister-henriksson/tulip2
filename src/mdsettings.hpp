
#ifndef MDSETTINGS_HPP
#define MDSETTINGS_HPP

#include <string>





// Set in main routine:
class MDSettings_common {
public:
  bool debug_forces;
  bool debug_pressure;
  bool report_step;
  bool   use_def_dump_xyz_fmt;
  std::string def_dump_xyz_fmt;


  MDSettings_common();
} ;



class MDSettings {
public:

  unsigned int seed;
  int ndump;
  double skint;

  double tstart;
  double tend;

  double Tstart;

  double dt;
  double max_dt;
  double max_dE;
  double max_dr;

  double btc_tau;
  double btc_T0;

  bool use_quench;
  bool use_Tcontrol;
  bool use_Pcontrol;

  double quench_tstart;
  double quench_rate;

  double bpc_tau;
  double bpc_P0;
  double bpc_scale;


  bool use_error_T_gt;
  bool use_error_dt_lt;
  bool use_error_boxlen_gt;

  double error_T_gt;
  double error_dt_lt;
  double error_boxlen_gt;


  bool is_ref_comp;
  bool heating_allowed;
  bool fixed_geometry;
  bool quench_always;

  MDSettings();
} ;

#endif

