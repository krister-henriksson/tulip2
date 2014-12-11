

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include <boost/format.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "potinfo.hpp"
#include "utils-string.hpp"
#include "utils-streamio.hpp"
#include "exiterrors.hpp"
#include "utils-math.hpp"

using namespace std;
using namespace utils;
using namespace exiterrors;
using boost::format;



void PotentialInformationFit::read_specs(string filename){
  ifstream fp;
  string line;
  vector<string> args;
  istringstream strbuf;
  int ns;
  double eps = std::numeric_limits<double>::epsilon();


  cout << "Reading specifications ... " << endl;



  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");


  int c;
  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;

    c=0;
    if      (args[0]=="prop") c=1;
    else if (args[0]=="pot")  c=2;


    if (c==1 || c==2){

      if (args[1]=="seed"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.seed;
	else if (c==2) strbuf >> specs_pot.seed;
	strbuf.clear();
      }
      else if (args[1]=="fitmet"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.fitmet;
	else if (c==2) strbuf >> specs_pot.fitmet;
	strbuf.clear();
      }
      else if (args[1]=="nitermax"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.nitermax;
	else if (c==2) strbuf >> specs_pot.nitermax;
	strbuf.clear();
      }
      else if (args[1]=="nitermin"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.nitermin;
	else if (c==2) strbuf >> specs_pot.nitermin;
	strbuf.clear();
      }
      else if (args[1]=="functolabs"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.functolabs;
	else if (c==2) strbuf >> specs_pot.functolabs;
	strbuf.clear();
      }
      else if (args[1]=="functolrel"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.functolrel;
	else if (c==2) strbuf >> specs_pot.functolrel;
	strbuf.clear();
      }
      else if (args[1]=="gradtolabs"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.gradtolabs;
	else if (c==2) strbuf >> specs_pot.gradtolabs;
	strbuf.clear();
      }
      else if (args[1]=="steptolabs"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.steptolabs;
	else if (c==2) strbuf >> specs_pot.steptolabs;
	strbuf.clear();
      }
      else if (args[1]=="steptolrel"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.steptolrel;
	else if (c==2) strbuf >> specs_pot.steptolrel;
	strbuf.clear();
      }
      else if (args[1]=="dogleg_radius"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.dogleg_radius;
	else if (c==2) strbuf >> specs_pot.dogleg_radius;
	strbuf.clear();
      }
      else if (args[1]=="dogleg_minradius"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.dogleg_minradius;
	else if (c==2) strbuf >> specs_pot.dogleg_minradius;
	strbuf.clear();
      }
      else if (args[1]=="simplex_delta"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.simplex_delta;
	else if (c==2) strbuf >> specs_pot.simplex_delta;
	strbuf.clear();
      }
      else if (args[1]=="barrier_scale"){
	strbuf.str(args[2]);
	if      (c==1) strbuf >> specs_prop.barrier_scale;
	else if (c==2) strbuf >> specs_pot.barrier_scale;
	strbuf.clear();
      }

    }



    if (args[0]=="prop" && args[1]=="ref" && args[2]=="mds" &&
	args[3]=="prop" && args[4]=="mds"){
      // prop:ref:mds = prop:mds
      cout << "" << endl;
      cout << "NOTE: Copying mds settings for reference compounds from the general mds settings" << endl;
      cout << "" << endl;
      specs_prop.mds_specs_ref = specs_prop.mds_specs;
    }



    // Only for fitting of properties:
    if (c==1){

      if (args[1]=="lattol"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.lattol;
	strbuf.clear();
      }

      else if (args[1]=="BM_rel_sys"){
	string ri;
	strbuf.str(args[2]); strbuf >> ri; strbuf.clear();
	if (ri[0]=='y' || ri[0]=='Y' || ri[0]=='t' || ri[0]=='T')
	  specs_prop.BM_rel_sys=true;
	else
	  specs_prop.BM_rel_sys=false;
      }
      else if (args[1]=="BM_fmin"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.BM_fmin;
	strbuf.clear();
      }
      else if (args[1]=="BM_fmax"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.BM_fmax;
	strbuf.clear();
      }
      else if (args[1]=="BM_Nf"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.BM_Nf;
	strbuf.clear();
      }
      else if (args[1]=="BM_ef"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.BM_ef;
	strbuf.clear();
      }
      
      else if (args[1]=="C_rel_sys"){
	string ri;
	strbuf.str(args[2]); strbuf >> ri; strbuf.clear();
	if (ri[0]=='y' || ri[0]=='Y' || ri[0]=='t' || ri[0]=='T')
	  specs_prop.C_rel_sys=true;
	else
	  specs_prop.C_rel_sys=false;
      }
      else if (args[1]=="C_fmin"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.C_fmin;
	strbuf.clear();
      }
      else if (args[1]=="C_fmax"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.C_fmax;
	strbuf.clear();
      }
      else if (args[1]=="C_Nf"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.C_Nf;
	strbuf.clear();
      }
      else if (args[1]=="C_ef"){
	strbuf.str(args[2]);
	strbuf >> specs_prop.C_ef;
	strbuf.clear();
      }






      // Only for handling of reference compounds:
      if ( (args[1]=="ref" &&
	    args[2][0]=='m' && args[2][1]=='d' && args[2][2]=='s')
	   ||
	   (args[1][0]=='m' && args[1][1]=='d' && args[1][2]=='s') ){

	int t=1, id=1;
	if (args[1]=="ref"){ t=0; id=2; }

	if      (args[id]=="mds_skint"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.skint;
	  else      strbuf >> specs_prop.mds_specs_ref.skint;
	  strbuf.clear();
	}
	else if (args[id]=="mds_seed"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.seed;
	  else      strbuf >> specs_prop.mds_specs_ref.seed;
	  strbuf.clear();
	}
	else if (args[id]=="mds_ndump"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.ndump;
	  else      strbuf >> specs_prop.mds_specs_ref.ndump;
	  strbuf.clear();
	}
	else if (args[id]=="mds_tstart"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.tstart;
	  else      strbuf >> specs_prop.mds_specs_ref.tstart;
	  strbuf.clear();
	}
	else if (args[id]=="mds_tend"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.tend;
	  else      strbuf >> specs_prop.mds_specs_ref.tend;
	  strbuf.clear();
	}
	else if (args[id]=="mds_Tstart"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.Tstart;
	  else      strbuf >> specs_prop.mds_specs_ref.Tstart;
	  strbuf.clear();
	}
	else if (args[id]=="mds_dt"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.dt;
	  else      strbuf >> specs_prop.mds_specs_ref.dt;
	  strbuf.clear();
	}
	else if (args[id]=="mds_max_dt"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.max_dt;
	  else      strbuf >> specs_prop.mds_specs_ref.max_dt;
	  strbuf.clear();
	}
	else if (args[id]=="mds_max_dE"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.max_dE;
	  else      strbuf >> specs_prop.mds_specs_ref.max_dE;
	  strbuf.clear();
	}
	else if (args[id]=="mds_max_dr"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.max_dr;
	  else      strbuf >> specs_prop.mds_specs_ref.max_dr;
	  strbuf.clear();
	}

	else if (args[id]=="mds_btc_tau"){
	  strbuf.str(args[id+1]);
	  if (t==1){
	    strbuf >> specs_prop.mds_specs.btc_tau;
	    if (specs_prop.mds_specs.btc_tau<0.0 ||
		abs(specs_prop.mds_specs.btc_tau)<eps)
	      specs_prop.mds_specs.use_Tcontrol = false;
	    else
	      specs_prop.mds_specs.use_Tcontrol = true;
	  }
	  else {
	    strbuf >> specs_prop.mds_specs_ref.btc_tau;
	    if (specs_prop.mds_specs_ref.btc_tau<0.0 ||
		abs(specs_prop.mds_specs_ref.btc_tau)<eps)
	      specs_prop.mds_specs_ref.use_Tcontrol = false;
	    else
	      specs_prop.mds_specs_ref.use_Tcontrol = true;
	  }
	  strbuf.clear();
	}

	else if (args[id]=="mds_btc_T0"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.btc_T0;
	  else      strbuf >> specs_prop.mds_specs_ref.btc_T0;
	  strbuf.clear();
	  if (t==1) specs_prop.mds_specs.use_Tcontrol = true;
	  else      specs_prop.mds_specs_ref.use_Tcontrol = true;
	}
	else if (args[id]=="mds_bpc_tau"){
	  strbuf.str(args[id+1]);
	  if (t==1){
	    strbuf >> specs_prop.mds_specs.bpc_tau;
	    if (specs_prop.mds_specs.bpc_tau<0.0 ||
		abs(specs_prop.mds_specs.bpc_tau)<eps)
	      specs_prop.mds_specs.use_Pcontrol = false;
	    else
	      specs_prop.mds_specs.use_Pcontrol = true;
	  }
	  else {
	    strbuf >> specs_prop.mds_specs_ref.bpc_tau;
	    if (specs_prop.mds_specs_ref.bpc_tau<0.0 ||
		abs(specs_prop.mds_specs_ref.bpc_tau)<eps)
	      specs_prop.mds_specs_ref.use_Pcontrol = true;
	  }
	  strbuf.clear();
	}

	else if (args[id]=="mds_bpc_P0"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.bpc_P0;
	  else      strbuf >> specs_prop.mds_specs_ref.bpc_P0;
	  strbuf.clear();
	  if (t==1) specs_prop.mds_specs.use_Pcontrol = true;
	  else      specs_prop.mds_specs_ref.use_Pcontrol = true;
	}
	else if (args[id]=="mds_bpc_scale"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.bpc_scale;
	  else      strbuf >> specs_prop.mds_specs_ref.bpc_scale;
	  strbuf.clear();
	  if (t==1) specs_prop.mds_specs.use_Pcontrol = true;
	  else      specs_prop.mds_specs_ref.use_Pcontrol = true;
	}
	else if (args[id]=="mds_quench_tstart"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.quench_tstart;
	  else      strbuf >> specs_prop.mds_specs_ref.quench_tstart;
	  strbuf.clear();
	  if (t==1) specs_prop.mds_specs.use_quench = true;
	  else      specs_prop.mds_specs_ref.use_quench = true;
	}
	else if (args[id]=="mds_quench_rate"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.quench_rate;
	  else      strbuf >> specs_prop.mds_specs_ref.quench_rate;
	  strbuf.clear();
	  if (t==1) specs_prop.mds_specs.use_quench = true;
	  else      specs_prop.mds_specs_ref.use_quench = true;
	}

	else if (args[id]=="mds_error_T_gt"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.error_T_gt;
	  else      strbuf >> specs_prop.mds_specs_ref.error_T_gt;
	  strbuf.clear();
	  if (t==1) specs_prop.mds_specs.use_error_T_gt = true;
	  else      specs_prop.mds_specs_ref.use_error_T_gt = true;
	}
	else if (args[id]=="mds_error_dt_lt"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.error_dt_lt;
	  else      strbuf >> specs_prop.mds_specs_ref.error_dt_lt;
	  strbuf.clear();
	  if (t==1) specs_prop.mds_specs.use_error_dt_lt = true;
	  else      specs_prop.mds_specs_ref.use_error_dt_lt = true;
	}
	else if (args[id]=="mds_error_boxlen_gt"){
	  strbuf.str(args[id+1]);
	  if (t==1) strbuf >> specs_prop.mds_specs.error_boxlen_gt;
	  else      strbuf >> specs_prop.mds_specs_ref.error_boxlen_gt;
	  strbuf.clear();
	  if (t==1) specs_prop.mds_specs.use_error_boxlen_gt = true;
	  else      specs_prop.mds_specs_ref.use_error_boxlen_gt = true;
	}


      }

    }


    if (! fp) break; 
  }
  fp.close();
  fp.clear();


  cout << "Read-in of specifications completed." << endl;

}

