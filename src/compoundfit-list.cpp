

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include <boost/format.hpp>

#include <cstdio>

#include "constants.hpp"
#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "utils-string.hpp"
#include "utils-streamio.hpp"
#include "utils-errors.hpp"

#include "elem-iacs.hpp"
//#include "specs-fit-prop-pot.hpp"
//#include "compound.hpp"
#include "physconst.hpp"


#include "compoundfit-list.hpp"
#include "lattice-simple.hpp"


#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;



using namespace utils;
using namespace constants;
using boost::format;





CompoundListFit::CompoundListFit(Elements & el,
				 MDSettings & mds_specs_general,
				 std::string filename)
  : elem(el), ncompounds(0)
{

  std::ifstream fp;
  std::ofstream fpo;
  std::string line;
  std::vector<std::string> args, opts, opts_special;
  std::istringstream strbuf;
  int ns, nlats, ilat, nelem;
  double td;
  Vector< Vector3<double> > bv;
  int itr, its, ot;
  bool reading_latinfo, elem_ok;
  int i, j, k, p;

  double eps = std::numeric_limits<double>::epsilon();



  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");


  nlats = 0;
  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;


    if (args[0]=="LAT") nlats++;


    if (! fp) break;
  }
  fp.close();
  fp.clear();

  std::cout << "Detected " << nlats << " compounds." << std::endl;

  compounds.resize(nlats);
  ncompounds = compounds.size();


  opts.resize(0);
  opts.push_back("a");
  opts.push_back("b");
  opts.push_back("c");
  opts.push_back("bpa");
  opts.push_back("cpa");
  opts.push_back("r0");
  opts.push_back("angle_ab");
  opts.push_back("angle_ac");
  opts.push_back("angle_bc");
  opts.push_back("Vatom");
  opts.push_back("Vat");
  opts.push_back("V0");
  opts.push_back("Ecoh");
  opts.push_back("Ec");
  opts.push_back("Ecoh_delta");
  opts.push_back("Ec_delta");
  opts.push_back("Ec");
  opts.push_back("E0");
  opts.push_back("Eatom");
  opts.push_back("Eatom_delta");
  opts.push_back("Eat");
  opts.push_back("Eat_delta");
  opts.push_back("Emix");
  opts.push_back("B");
  opts.push_back("B");
  opts.push_back("Bp");
  opts.push_back("dB/dP");
  opts.push_back("Fmax");
  opts.push_back("F_max");
  opts.push_back("Pmax");
  opts.push_back("P_max");
  opts.push_back("displmax");
  opts.push_back("displ_max");

  opts_special.resize(0);
  opts_special.push_back("C");



  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  ilat = 0;
  reading_latinfo = false;

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;


    if (args[0]=="LAT"){
      if (reading_latinfo){
	ilat++;
      }
      else {
	reading_latinfo = true;
	k = compounds.size();
	if (k<=ilat)
	  compounds.resize(k+1);
      }


      // Initialization:
      compounds[ilat].mds_specs = mds_specs_general;

      // Since we are reading this compound from a file, it is not a reference compound:
      compounds[ilat].mds_specs.is_ref_comp = false;


    }


    // *******************************************************************
    // Compound properties, needed for e.g. MD cell construction in
    // conjunction with the information in the compound's LAT file:
    // *******************************************************************

    else if (args[0]=="name"){
      compounds[ilat].name = args[1];
    }
    else if (args[0]=="crystalname"){
      compounds[ilat].crystalname = args[1];
    }
    else if (args[0]=="file"){
      compounds[ilat].filename = args[1];
    }
    else if (args[0]=="elements"){
      compounds[ilat].elemnames.resize(0);
      for (i=1; i<ns; ++i){
	if (args[i][0]=='#') break;
	compounds[ilat].elemnames.push_back(args[i]);
      }
      nelem = compounds[ilat].elemnames.size();
      compounds[ilat].nelem = nelem;


      // Check each supplied element, if it is known from before:
      for (i=0; i<nelem; ++i){
	elem_ok=false;
	for (j=0; j<elem.nelem(); ++j){
	  if (compounds[ilat].elemnames[i] == elem.idx2name(j)){
	    // found the element
	    elem_ok=true;
	    break;
	  }
	}
	if (! elem_ok){
	  aborterror("Element " + compounds[ilat].elemnames[i] + " was not found "
		     + "in the list of known element names. Exiting.");
	}
      }
    }
    else if (args[0]=="csystem"){
      if      ( args[1][0] == 'c' || args[1][0] == 'C' )
	compounds[ilat].csystem = "cubic";
      else if ( args[1][0] == 'h' || args[1][0] == 'H' )
	compounds[ilat].csystem = "hexagonal";
      else if ( args[1][0] == 'o' || args[1][0] == 'O' )
	compounds[ilat].csystem = "orthorombic";
      else if ( (args[1][0] == 't' || args[1][0] == 'T') &&
		(args[1][1] == 'r' || args[1][1] == 'R') &&
		(args[1][2] == 'i' || args[1][2] == 'I') &&
		(args[1][3] == 'c' || args[1][3] == 'C') )
	compounds[ilat].csystem = "triclinic";
      else if ( (args[1][0] == 't' || args[1][0] == 'T') &&
		(args[1][1] == 'r' || args[1][1] == 'R') &&
		(args[1][2] == 'i' || args[1][2] == 'I') &&
		(args[1][3] == 'g' || args[1][3] == 'G') )
	compounds[ilat].csystem = "trigonal";
      else if ( (args[1][0] == 't' || args[1][0] == 'T') &&
		(args[1][1] == 'e' || args[1][1] == 'E') )
	compounds[ilat].csystem = "tetragonal";
      else if ( args[1][0] == 'm' || args[1][0] == 'M' )
	compounds[ilat].csystem = "monoclinic";
      else {
	aborterror("ERROR: Unknown crystal system " + args[1] + ". Exiting.");
      }
    }

    else if (args[0]=="frc_file"){
      compounds[ilat].filename_frc = args[1];
      compounds[ilat].prop_use.frc = true;
    }
    else if (args[0]=="frc_use"){
      if (args[1][0]=='y' || args[1][0]=='Y' ||
	  args[1][0]=='T' || args[1][0]=='T'){
	compounds[ilat].prop_use.frc = true;
      }
      else
	compounds[ilat].prop_use.frc = false;
    }
    else if (args[0]=="frc_use_u"){
      compounds[ilat].use_u.frc = true;
      compounds[ilat].use_w.frc = false;
    }
    else if (args[0]=="frc_use_w"){
      compounds[ilat].use_u.frc = false;
      compounds[ilat].use_w.frc = true;
    }

    else if (args[0]=="Ndesired"){
      strbuf.str(args[1]); strbuf >> compounds[ilat].Ndesired[0]; strbuf.clear();
      strbuf.str(args[2]); strbuf >> compounds[ilat].Ndesired[1]; strbuf.clear();
      strbuf.str(args[3]); strbuf >> compounds[ilat].Ndesired[2]; strbuf.clear();
    }
    else if (args[0]=="Neven_desired"){
      compounds[ilat].Neven_desired[0] = bool_in_string(args[1]);
      compounds[ilat].Neven_desired[1] = bool_in_string(args[2]);
      compounds[ilat].Neven_desired[2] = bool_in_string(args[3]);
    }
    else if (args[0]=="Nodd_desired"){
      compounds[ilat].Nodd_desired[0] = bool_in_string(args[1]);
      compounds[ilat].Nodd_desired[1] = bool_in_string(args[2]);
      compounds[ilat].Nodd_desired[2] = bool_in_string(args[3]);
    }

    else if (args[0]=="Ecoh_delta_ref" || args[0]=="Ecoh_delta_refcomp" || 
	     args[0]=="Ecoh_delta_ref_comp"){
      compounds[ilat].Ecoh_delta_refcomp = bool_in_string(args[1]);
    }




    // *******************************************************************
    // Compound-specific conditions for e.g. MD relaxation:
    // *******************************************************************

    else if (args[0]=="option"){

      if      (args[1]=="no_heating"){
	compounds[ilat].mds_specs.heating_allowed = false;
      }
      else if (args[1]=="fixed_geometry"){
	compounds[ilat].mds_specs.fixed_geometry = true;
      }
      else if (args[1]=="quench_always"){
	compounds[ilat].mds_specs.quench_always = true;
      }

    }

    else if (args[0][0]=='m' && args[0][1]=='d' && args[0][2]=='s'){

      if      (args[0]=="mds_skint"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.skint;
	strbuf.clear();
      }
      else if (args[0]=="mds_seed"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.seed;
	strbuf.clear();
      }
      else if (args[0]=="mds_ndump"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.ndump;
	strbuf.clear();
      }
      else if (args[0]=="mds_tstart"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.tstart;
	strbuf.clear();
      }
      else if (args[0]=="mds_tend"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.tend;
	strbuf.clear();
      }
      else if (args[0]=="mds_Tstart"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.Tstart;
	strbuf.clear();
      }
      else if (args[0]=="mds_dt"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.dt;
	strbuf.clear();
      }
      else if (args[0]=="mds_max_dt"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.max_dt;
	strbuf.clear();
      }
      else if (args[0]=="mds_max_dE"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.max_dE;
	strbuf.clear();
      }
      else if (args[0]=="mds_max_dr"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.max_dr;
	strbuf.clear();
      }
      else if (args[0]=="mds_btc_tau"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.btc_tau;
	if (compounds[ilat].mds_specs.btc_tau<0.0 ||
	    abs(compounds[ilat].mds_specs.btc_tau)<eps)
	  compounds[ilat].mds_specs.use_Tcontrol = false;
	else
	  compounds[ilat].mds_specs.use_Tcontrol = true;
	strbuf.clear();
      }

      else if (args[0]=="mds_btc_T0"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.btc_T0;
	strbuf.clear();
      }
      else if (args[0]=="mds_bpc_tau"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.bpc_tau;

	if (compounds[ilat].mds_specs.bpc_tau<0.0 ||
	    abs(compounds[ilat].mds_specs.bpc_tau)<eps)
	  compounds[ilat].mds_specs.use_Pcontrol = false;
	else
	  compounds[ilat].mds_specs.use_Pcontrol = true;
	strbuf.clear();

      }
      else if (args[0]=="mds_bpc_P0"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.bpc_P0;
	strbuf.clear();
      }
      else if (args[0]=="mds_bpc_scale"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.bpc_scale;
	strbuf.clear();
      }
      else if (args[0]=="mds_quench_tstart"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.quench_tstart;
	strbuf.clear();
      }
      else if (args[0]=="mds_quench_rate"){
	strbuf.str(args[1]);
	strbuf >> compounds[ilat].mds_specs.quench_rate;
	strbuf.clear();
      }

      else if (args[0]=="mds_error_T_gt"){
	strbuf.str(args[1]); strbuf >> compounds[ilat].mds_specs.error_T_gt; strbuf.clear();
	compounds[ilat].mds_specs.use_error_T_gt = true;
      }
      else if (args[0]=="mds_error_dt_lt"){
	strbuf.str(args[1]); strbuf >> compounds[ilat].mds_specs.error_dt_lt; strbuf.clear();
	compounds[ilat].mds_specs.use_error_dt_lt = true;
      }
      else if (args[0]=="mds_error_boxlen_gt"){
	strbuf.str(args[1]); strbuf >> compounds[ilat].mds_specs.error_boxlen_gt; strbuf.clear();
	compounds[ilat].mds_specs.use_error_boxlen_gt = true;
      }

      // OTHER OPTIONS POSSIBLE TO ADD HERE


    }


    // *******************************************************************
    // Compound properties to use as targets when fitting:
    // *******************************************************************

    else {
      // args[0] is being investigated.
      // It could be e.g. 'a', 'w_a', or 'u_a'

      itr = -1;
      its = -1;
      ot  = -1;

      // Regular options, matches an entire string:
      for (unsigned int i=0; i<opts.size(); ++i){
	if      (args[0]==opts[i]){
	  itr = i; ot = 1; break;
	}
	else if (args[0]==("w_" + opts[i])){
	  itr = i; ot = 2; break;
	}
	else if (args[0]==("u_" + opts[i])){
	  itr = i; ot = 3; break;
	}
      }
      // Special options, starting with a single character and then followed by
      // others:
      if (itr<0){
	for (unsigned int i=0; i<opts_special.size(); ++i){
	  if      (args[0][0]==opts_special[i][0]){
	    its = i; ot = 1; break;
	  }
	  else if (args[0][0]=='w' &&
		   args[0][1]=='_' &&
		   args[0][2]==opts_special[i][0]){
	    its = i; ot = 2; break;
	  }
	  else if (args[0][0]=='u' &&
		   args[0][1]=='_' &&
		   args[0][2]==opts_special[i][0]){
	    its = i; ot = 3; break;
	  }
	}
      }




      std::string match="n";
      if (itr>=0) match = opts[itr];
      if (its>=0) match = opts_special[its];
      if (itr>=0 || its>=0){
	strbuf.str(args[1]); strbuf >> td; strbuf.clear();
      }

      //cout << "match = " << match << " match[0] = " << match[0] << endl;

      
      if (match=="a"){
	if (ot==1){
	  compounds[ilat].prop_readin.a = td;
	  compounds[ilat].prop_use.a    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.a = td; compounds[ilat].use_w.a = true;  compounds[ilat].use_u.a = false; }
	else if (ot==3) { compounds[ilat].prop_u.a = td; compounds[ilat].use_w.a = false; compounds[ilat].use_u.a = true; }
      }
      else if (match=="b"){
	if (ot==1){
	  compounds[ilat].prop_readin.b = td;
	  compounds[ilat].prop_use.b    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.b = td; compounds[ilat].use_w.b = true;  compounds[ilat].use_u.b = false; }
	else if (ot==3) { compounds[ilat].prop_u.b = td; compounds[ilat].use_w.b = false; compounds[ilat].use_u.b = true; }
      }
      else if (match=="c"){
	if (ot==1){
	  compounds[ilat].prop_readin.c = td;
	  compounds[ilat].prop_use.c    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.c = td; compounds[ilat].use_w.c = true;  compounds[ilat].use_u.c = false; }
	else if (ot==3) { compounds[ilat].prop_u.c = td; compounds[ilat].use_w.c = false; compounds[ilat].use_u.c = true; }
      }
      else if (match=="bpa"){
	if (ot==1){
	  compounds[ilat].prop_readin.bpa = td;
	  compounds[ilat].prop_use.bpa    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.bpa = td; compounds[ilat].use_w.bpa = true;  compounds[ilat].use_u.bpa = false; }
	else if (ot==3) { compounds[ilat].prop_u.bpa = td; compounds[ilat].use_w.bpa = false; compounds[ilat].use_u.bpa = true; }
      }
      else if (match=="cpa"){
	if (ot==1){
	  compounds[ilat].prop_readin.cpa = td;
	  compounds[ilat].prop_use.cpa    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.cpa = td; compounds[ilat].use_w.cpa = true;  compounds[ilat].use_u.cpa = false; }
	else if (ot==3) { compounds[ilat].prop_u.cpa = td; compounds[ilat].use_w.cpa = false; compounds[ilat].use_u.cpa = true; }
      }
      else if (match=="r0"){
	if (ot==1){
	  compounds[ilat].prop_readin.r0 = td;
	  compounds[ilat].prop_use.r0    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.r0 = td; compounds[ilat].use_w.r0 = true;  compounds[ilat].use_u.r0 = false; }
	else if (ot==3) { compounds[ilat].prop_u.r0 = td; compounds[ilat].use_w.r0 = false; compounds[ilat].use_u.r0 = true; }
      }
      else if (match=="angle_ab"){
	if (ot==1){
	  compounds[ilat].prop_readin.angle_ab = td;
	  compounds[ilat].prop_use.angle_ab    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.angle_ab = td; compounds[ilat].use_w.angle_ab = true;  compounds[ilat].use_u.angle_ab = false; }
	else if (ot==3) { compounds[ilat].prop_u.angle_ab = td; compounds[ilat].use_w.angle_ab = false; compounds[ilat].use_u.angle_ab = true; }
      }
      else if (match=="angle_ac"){
	if (ot==1){
	  compounds[ilat].prop_readin.angle_ac = td;
	  compounds[ilat].prop_use.angle_ac    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.angle_ac = td; compounds[ilat].use_w.angle_ac = true;  compounds[ilat].use_u.angle_ac = false; }
	else if (ot==3) { compounds[ilat].prop_u.angle_ac = td; compounds[ilat].use_w.angle_ac = false; compounds[ilat].use_u.angle_ac = true; }
      }
      else if (match=="angle_bc"){
	if (ot==1){
	  compounds[ilat].prop_readin.angle_bc = td;
	  compounds[ilat].prop_use.angle_bc    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.angle_bc = td; compounds[ilat].use_w.angle_bc = true;  compounds[ilat].use_u.angle_bc = false; }
	else if (ot==3) { compounds[ilat].prop_u.angle_bc = td; compounds[ilat].use_w.angle_bc = false; compounds[ilat].use_u.angle_bc = true; }
      }
      else if (match=="V0" || match=="Vat" || match=="Vatom"){
	if (ot==1){
	  compounds[ilat].prop_readin.Vatom = td;
	  compounds[ilat].prop_use.Vatom    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Vatom = td; compounds[ilat].use_w.Vatom = true;  compounds[ilat].use_u.Vatom = false; }
	else if (ot==3) { compounds[ilat].prop_u.Vatom = td; compounds[ilat].use_w.Vatom = false; compounds[ilat].use_u.Vatom = true; }
      }
      else if (match=="Ec" || match=="Ecoh" || match=="Eat" || match=="Eatom"){
	if (ot==1){
	  compounds[ilat].prop_readin.Ecoh = td;
	  compounds[ilat].prop_use.Ecoh    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Ecoh = td; compounds[ilat].use_w.Ecoh = true;  compounds[ilat].use_u.Ecoh = false; }
	else if (ot==3) { compounds[ilat].prop_u.Ecoh = td; compounds[ilat].use_w.Ecoh = false; compounds[ilat].use_u.Ecoh = true; }
      }
      else if (match=="Ec_delta" || match=="Ecoh_delta" || match=="Eat_delta" || match=="Eatom_delta"){
	if (ot==1){
	  compounds[ilat].prop_readin.Ecoh_delta = td;
	  compounds[ilat].prop_use.Ecoh_delta    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Ecoh_delta = td; compounds[ilat].use_w.Ecoh_delta = true;  compounds[ilat].use_u.Ecoh_delta = false; }
	else if (ot==3) { compounds[ilat].prop_u.Ecoh_delta = td; compounds[ilat].use_w.Ecoh_delta = false; compounds[ilat].use_u.Ecoh_delta = true; }
      }

      else if (match=="Emix"){
	if (ot==1){
	  compounds[ilat].prop_readin.Emix = td;
	  compounds[ilat].prop_use.Emix    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Emix = td; compounds[ilat].use_w.Emix = true;  compounds[ilat].use_u.Emix = false; }
	else if (ot==3) { compounds[ilat].prop_u.Emix = td; compounds[ilat].use_w.Emix = false; compounds[ilat].use_u.Emix = true; }
      }
      else if (match=="B"){
	if (ot==1){
	  compounds[ilat].prop_readin.B = td;
	  compounds[ilat].prop_use.B    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.B = td; compounds[ilat].use_w.B = true;  compounds[ilat].use_u.B = false; }
	else if (ot==3) { compounds[ilat].prop_u.B = td; compounds[ilat].use_w.B = false; compounds[ilat].use_u.B = true; }
      }
      else if (match=="Bp" || match=="dB/dP"){
	if (ot==1){
	  compounds[ilat].prop_readin.Bp = td;
	  compounds[ilat].prop_use.Bp    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Bp = td; compounds[ilat].use_w.Bp = true;  compounds[ilat].use_u.Bp = false; }
	else if (ot==3) { compounds[ilat].prop_u.Bp = td; compounds[ilat].use_w.Bp = false; compounds[ilat].use_u.Bp = true; }
      }

      else if (match=="Fmax" || match=="F_max"){
	if (ot==1){
	  compounds[ilat].prop_readin.Fmax = td;
	  compounds[ilat].prop_use.Fmax    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Fmax = td; compounds[ilat].use_w.Fmax = true;  compounds[ilat].use_u.Fmax = false; }
	else if (ot==3) { compounds[ilat].prop_u.Fmax = td; compounds[ilat].use_w.Fmax = false; compounds[ilat].use_u.Fmax = true; }
      }
      else if (match=="Pmax" || match=="P_max"){
	if (ot==1){
	  compounds[ilat].prop_readin.Pmax = td;
	  compounds[ilat].prop_use.Pmax    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.Pmax = td; compounds[ilat].use_w.Pmax = true;  compounds[ilat].use_u.Pmax = false; }
	else if (ot==3) { compounds[ilat].prop_u.Pmax = td; compounds[ilat].use_w.Pmax = false; compounds[ilat].use_u.Pmax = true; }
      }
      else if (match=="displmax" || match=="displ_max"){
	if (ot==1){
	  compounds[ilat].prop_readin.displmax = td;
	  compounds[ilat].prop_use.displmax    = true;
	}
	else if (ot==2) { compounds[ilat].prop_w.displmax = td; compounds[ilat].use_w.displmax = true;  compounds[ilat].use_u.displmax = false; }
	else if (ot==3) { compounds[ilat].prop_u.displmax = td; compounds[ilat].use_w.displmax = false; compounds[ilat].use_u.displmax = true; }
      }


      // Special parameters:
      else if (its>=0 && match=="C"){
	char k=0,p=0;
	if (ot==1){
	  k = args[0][1];
	  p = args[0][2];
	}
	else if (ot==2 || ot==3){
	  k = args[0][3];
	  p = args[0][4];
	}
	// Now get integers:
	int u, v;
	std::stringstream sstream;
	sstream.clear(); sstream << k; sstream >> u;
	sstream.clear(); sstream << p; sstream >> v;
	//cout << "u v " << u << " " << v << endl;

	if (u>=1 && u<=6 && v>=1 && v<=6){
	  //cout << "C (again)" << u << v << endl;

	  if (ot==1){
	    compounds[ilat].prop_readin.C.elem(u-1,v-1) = td;
	    compounds[ilat].prop_use.C.elem(u-1,v-1)    = true;
	    //cout << "use C" << u << v << " : " << compounds[ilat].prop_use.C.elem(u-1,v-1) << endl;
	  }
	  else if (ot==2) {
	    compounds[ilat].prop_w.C.elem(u-1,v-1) = td;
	    compounds[ilat].use_w.C.elem(u-1,v-1) = true;
	    compounds[ilat].use_u.C.elem(u-1,v-1) = false;
	  }
	  else if (ot==3) {
	    compounds[ilat].prop_u.C.elem(u-1,v-1) = td;
	    compounds[ilat].use_w.C.elem(u-1,v-1) = false;
	    compounds[ilat].use_u.C.elem(u-1,v-1) = true;
	  }
	}
      }


    }





    if (!fp) break; 
  }
  fp.clear();
  fp.close();
  
  std::cout << "Compounds information collected." << std::endl;





  /*
  if (compounds.size() != ilat)
    compounds.resize(ilat+1);
  */
  // ncompounds = compounds.size();



  std::cout << "Number of compounds read: " << compounds.size() << std::endl;
  std::cout << "Debugging compounds ..." << std::endl;



  int nEcohref=0;
  int nEcohdelta=0;
  /* -----------------------------------------------------------------------------
     Debugging of settings:
     ----------------------------------------------------------------------------- */
  for (ilat=0; ilat<compounds.size(); ++ilat){


    if (compounds[ilat].Ecoh_delta_refcomp) nEcohref++;
    if (compounds[ilat].prop_use.Ecoh_delta) nEcohdelta++;

    if (compounds[ilat].prop_use.Ecoh &&
	compounds[ilat].prop_use.Ecoh_delta)
      aborterror("ERROR: Both Ecoh and Ecoh_delta are used for compound "
		 + compounds[ilat].name + ". At most one of these can be used!");

    if (compounds[ilat].mds_specs.fixed_geometry){
      compounds[ilat].prop_use.displmax   = false;
      compounds[ilat].prop_use.a = false;
      compounds[ilat].prop_use.b = false;
      compounds[ilat].prop_use.c = false;
      compounds[ilat].prop_use.bpa = false;
      compounds[ilat].prop_use.cpa = false;
      compounds[ilat].prop_use.angle_ab = false;
      compounds[ilat].prop_use.angle_ac = false;
      compounds[ilat].prop_use.angle_bc = false;
      compounds[ilat].prop_use.r0 = false;
      compounds[ilat].prop_use.Vatom = false;
      compounds[ilat].prop_use.frc = false;
    }

    if (compounds[ilat].prop_use.frc &&
	compounds[ilat].filename_frc=="none"){
      aborterror("ERROR: No file specified for reading in the forces for compound "
		 + compounds[ilat].name);
    }
  }



  for (ilat=0; ilat<compounds.size(); ++ilat){

    // Make sure we are not using the same properties in different
    // disguises, since this will lead to problems in some fitting
    // algorithms:
    compounds[ilat].check_and_fix_uses();

  }



  if (nEcohdelta>0 && nEcohref==0)
    aborterror("ERROR: Ecoh_delta values are used, but there is no reference Ecoh compound!");

  if (nEcohref>=1 && nEcohdelta==0)
    std::cout << "Warning: Reference Ecoh compound is used, but no Ecoh_delta values." << std::endl;








  /* Some notes:

     - If a,b,c is mentioned in the geometry file they will be used. No checks or
     comparisons is made with the length of the primitive vectors in the LAT files.
   */


    
    
    

  /* ##############################################################################
     ##############################################################################
     Get lattice data for the read-in structures.
     ##############################################################################
     ##############################################################################    
  */





  int nc = compounds.size();
  std::cout << "Reading structural info for compounds ... " << std::endl;
  for (ilat=0; ilat<nc; ++ilat){
    compounds[ilat].read_structure(el);
    // calls finalize() internally to insert a,b,c

    //compounds[ilat].check_crystal_symm();
    

  }


  /* -----------------------------------------------------------------------------
     Some debugging
     ----------------------------------------------------------------------------- */


  for (ilat=0; ilat<nc; ++ilat){

    if (compounds[ilat].name == "none")
      aborterror("Error: No compound name specified for compound " + tostring(ilat) + ". Exiting.");

    if (compounds[ilat].filename == "none")
      aborterror("Error: No file name specified for compound " + tostring(ilat) + ". Exiting.");

    if (compounds[ilat].nelem == 0)
      aborterror("Error: No elements specified for compound " + tostring(ilat) + ". Exiting.");

    if (compounds[ilat].nbasis == 0)
      aborterror("Error: No basis specified for compound " + tostring(ilat) + ". Exiting.");



    if (!compounds[ilat].pbc[0] ||
	!compounds[ilat].pbc[1] ||
	!compounds[ilat].pbc[2]){
      compounds[ilat].prop_use.B  = false;
      compounds[ilat].prop_use.Bp = false;
      compounds[ilat].prop_use.Pmax = false;
      for (k=0; k<6; ++k)
	for (p=0; p<6; ++p)
	  compounds[ilat].prop_use.C.elem(k,p) = false;
    }


    for (j=0; j<compounds[ilat].nbasis; ++j){
      elem_ok=false;
      for (i=0; i<compounds[ilat].nelem; ++i){
	if (compounds[ilat].basis_elems[j] == compounds[ilat].elemnames[i]){
	  // found the element
	  elem_ok=true;
	  compounds[ilat].basis_types[j] = i;
	  break;
	}
      }
      if (! elem_ok){
	aborterror("Error: Compound " + tostring(ilat) + " basis vector " + tostring(j) + " "
		   "with element " + compounds[ilat].basis_elems[j] + " was not found "
		   + "in the list of known element names for this compound. Exiting.");
      }
    }
  }


  /* ##############################################################################
     ##############################################################################
     Get forces for the read-in structures.
     ##############################################################################
     ##############################################################################    
  */

  std::cout << "Reading forces for compounds ... " << std::endl;
  for (ilat=0; ilat<nc; ++ilat){
    if (compounds[ilat].prop_use.frc){
      std::cout << "Forces used for compound " << ilat+1 << " of " << nc << std::endl;
      compounds[ilat].read_forces();
    }
    else {
      std::cout << "Forces NOT used for compound " << ilat+1 << " of " << nc << std::endl;
    }
  }




  /* ##############################################################################
     ##############################################################################
     Lattice symmetries, elastic constants stuff
     ##############################################################################
     ##############################################################################    
  */


  latsymm(compounds);


  for (ilat=0; ilat<nc; ++ilat){

    compounds[ilat].check_and_fix_Cij();

  }







  // ####################################################################################
  // Normalize weights of all read-in compound properties:
  // ####################################################################################

  double wsum=0.0;
  nc = compounds.size();
  for (ilat=0; ilat<nc; ++ilat){

    if (compounds[ilat].prop_use.a && compounds[ilat].use_w.a)
      wsum += square( compounds[ilat].prop_w.a );

    if (compounds[ilat].prop_use.b && compounds[ilat].use_w.b)
      wsum += square( compounds[ilat].prop_w.b );

    if (compounds[ilat].prop_use.c && compounds[ilat].use_w.c)
      wsum += square( compounds[ilat].prop_w.c );

    if (compounds[ilat].prop_use.bpa && compounds[ilat].use_w.bpa)
      wsum += square( compounds[ilat].prop_w.bpa );

    if (compounds[ilat].prop_use.cpa && compounds[ilat].use_w.cpa)
      wsum += square( compounds[ilat].prop_w.cpa );

    if (compounds[ilat].prop_use.r0 && compounds[ilat].use_w.r0)
      wsum += square( compounds[ilat].prop_w.r0 );

    if (compounds[ilat].prop_use.angle_ab && compounds[ilat].use_w.angle_ab)
      wsum += square( compounds[ilat].prop_w.angle_ab );

    if (compounds[ilat].prop_use.angle_ac && compounds[ilat].use_w.angle_ac)
      wsum += square( compounds[ilat].prop_w.angle_ac );

    if (compounds[ilat].prop_use.angle_bc && compounds[ilat].use_w.angle_bc)
      wsum += square( compounds[ilat].prop_w.angle_bc );

    if (compounds[ilat].prop_use.Vatom && compounds[ilat].use_w.Vatom)
      wsum += square( compounds[ilat].prop_w.Vatom );

    if (compounds[ilat].prop_use.Ecoh && compounds[ilat].use_w.Ecoh)
      wsum += square( compounds[ilat].prop_w.Ecoh );

    if (compounds[ilat].prop_use.Ecoh_delta && compounds[ilat].use_w.Ecoh_delta)
      wsum += square( compounds[ilat].prop_w.Ecoh_delta );

    if (compounds[ilat].prop_use.Emix && compounds[ilat].use_w.Emix)
      wsum += square( compounds[ilat].prop_w.Emix );

    if (compounds[ilat].prop_use.B && compounds[ilat].use_w.B)
      wsum += square( compounds[ilat].prop_w.B );

    if (compounds[ilat].prop_use.Bp && compounds[ilat].use_w.Bp)
      wsum += square( compounds[ilat].prop_w.Bp );

    for (int k=0; k<6; ++k){
      for (int p=0; p<6; ++p){
	if (compounds[ilat].prop_use.C.elem(k,p) && compounds[ilat].use_w.C.elem(k,p))
	  wsum += square( compounds[ilat].prop_w.C.elem(k,p) );
      }
    }

    if (compounds[ilat].prop_use.Fmax && compounds[ilat].use_w.Fmax)
      wsum += square( compounds[ilat].prop_w.Fmax );

    if (compounds[ilat].prop_use.Pmax && compounds[ilat].use_w.Pmax)
      wsum += square( compounds[ilat].prop_w.Pmax );

    if (compounds[ilat].prop_use.displmax && compounds[ilat].use_w.displmax)
      wsum += square( compounds[ilat].prop_w.displmax );

    if (compounds[ilat].prop_use.frc){
      int nb = compounds[ilat].basis_elems.size();
      for (int iat=0; iat<nb; ++iat){
	for (int k=0; k<3; ++k){
	  if (compounds[ilat].use_w.frc)
	    wsum += square( compounds[ilat].prop_w.frc[iat][k] );
	}
      }
    }

  }
  wsum = 1.0/sqrt(wsum);
  for (ilat=0; ilat<nc; ++ilat){

    if (compounds[ilat].prop_use.a && compounds[ilat].use_w.a)
      compounds[ilat].prop_w.a *= wsum;

    if (compounds[ilat].prop_use.b && compounds[ilat].use_w.b)
      compounds[ilat].prop_w.b *= wsum;

    if (compounds[ilat].prop_use.c && compounds[ilat].use_w.c)
      compounds[ilat].prop_w.c *= wsum;

    if (compounds[ilat].prop_use.bpa && compounds[ilat].use_w.bpa)
      compounds[ilat].prop_w.bpa *= wsum;

    if (compounds[ilat].prop_use.cpa && compounds[ilat].use_w.cpa)
      compounds[ilat].prop_w.cpa *= wsum;

    if (compounds[ilat].prop_use.r0 && compounds[ilat].use_w.r0)
      compounds[ilat].prop_w.r0 *= wsum;

    if (compounds[ilat].prop_use.angle_ab && compounds[ilat].use_w.angle_ab)
      compounds[ilat].prop_w.angle_ab *= wsum;

    if (compounds[ilat].prop_use.angle_ac && compounds[ilat].use_w.angle_ac)
      compounds[ilat].prop_w.angle_ac *= wsum;

    if (compounds[ilat].prop_use.angle_bc && compounds[ilat].use_w.angle_bc)
      compounds[ilat].prop_w.angle_bc *= wsum;

    if (compounds[ilat].prop_use.Vatom && compounds[ilat].use_w.Vatom)
      compounds[ilat].prop_w.Vatom *= wsum;

    if (compounds[ilat].prop_use.Ecoh && compounds[ilat].use_w.Ecoh)
      compounds[ilat].prop_w.Ecoh *= wsum;

    if (compounds[ilat].prop_use.Ecoh_delta && compounds[ilat].use_w.Ecoh_delta)
      compounds[ilat].prop_w.Ecoh_delta *= wsum;

    if (compounds[ilat].prop_use.Emix && compounds[ilat].use_w.Emix)
      compounds[ilat].prop_w.Emix *= wsum;

    if (compounds[ilat].prop_use.B && compounds[ilat].use_w.B)
      compounds[ilat].prop_w.B *= wsum;

    if (compounds[ilat].prop_use.Bp && compounds[ilat].use_w.Bp)
      compounds[ilat].prop_w.Bp *= wsum;

    for (int k=0; k<6; ++k){
      for (int p=0; p<6; ++p){
	if (compounds[ilat].prop_use.C.elem(k,p) && compounds[ilat].use_w.C.elem(k,p))
	  compounds[ilat].prop_w.C.elem(k,p) *= wsum;
      }
    }

    if (compounds[ilat].prop_use.Fmax && compounds[ilat].use_w.Fmax)
      compounds[ilat].prop_w.Fmax *= wsum;

    if (compounds[ilat].prop_use.Pmax && compounds[ilat].use_w.Pmax)
      compounds[ilat].prop_w.Pmax *= wsum;

    if (compounds[ilat].prop_use.displmax && compounds[ilat].use_w.displmax)
      compounds[ilat].prop_w.displmax *= wsum;

    if (compounds[ilat].prop_use.frc){
      int nb = compounds[ilat].basis_elems.size();
      for (int iat=0; iat<nb; ++iat){
	for (int k=0; k<3; ++k){
	  if (compounds[ilat].use_w.frc)
	    compounds[ilat].prop_w.frc[iat][k] *= wsum;
	}
      }
    }

  }
  // ####################################################################################
  // ####################################################################################







  std::cout << "Compound information read-in completed." << std::endl;

  return;
}







int CompoundListFit::NData(){
  int N=0, i;

  for (i=0; i<compounds.size(); ++i)
    N += compounds[i].NData();

  return N;
}


