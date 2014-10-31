
#ifndef MDSETTINGS_HPP
#define MDSETTINGS_HPP




// Set in main routine:
class MDSettings_common {
public:
  bool debug_forces;
  bool debug_pressure;

  bool report_step;
  bool quick_mode;

  MDSettings_common();
} ;



class MDSettings {
public:

  double skint;
  unsigned int seed;
  int ndump;

  double tstart;
  double tend;

  double Tstart;

  double dt;
  double max_dt;
  double max_dE;
  double max_dr;

  bool use_quench;
  bool use_Tcontrol;
  bool use_Pcontrol;

  double btc_tau;
  double btc_T0;

  double quench_tstart;
  double quench_rate;

  double bpc_tau;
  double bpc_P0;
  double bpc_scale;

  bool is_ref_comp;
  bool heating_allowed;
  bool fixed_geometry;
  bool quench_always;

  MDSettings();
} ;

#endif

