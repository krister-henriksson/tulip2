



#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>

#include <boost/format.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "utils-string.hpp"
#include "utils-streamio.hpp"
#include "utils-errors.hpp"
#include "utils-math.hpp"

#include "potinfo.hpp"
#include "param.hpp"
#include "helpfuns.hpp"


using namespace utils;
using boost::format;




// #################################################################################
// #################################################################################
// #################################################################################
// #################################################################################




void PotentialInformation::read_info(string filename){
  std::ifstream fp;
  std::ofstream fpo;
  std::string line, ts, ts1, ts2, ts3, tsm, tsi, s1, s2, tso, tsoi, ts0;
  std::vector<std::string> args;
  std::istringstream strbuf;
  double td;
  int tl, i1,i2,i3, ivec;
  std::string potname;
  int ns;

  std::cout << "Reading general information about potentials ..." << std::endl;


  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information about elements only.
     The rest of the data will be obtained from subsequent readings of the file.

     ############################################################################
     ############################################################################ */



  std::cout << "Pass 1 ... " << std::endl;

  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  line.resize(0);

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=" );
    if (ns==0 && fp) continue;


    /*
    for (int i=0; i<ns; ++i)
      std::cout << args[i] << std::endl;
    */

    if (args[0]=="elem"){
      strbuf.str(args[2]); strbuf >> ts; strbuf.clear();
      //std::cout << ts << std::endl;
      elem.add_elem(ts);
    }
    /*
    else if (args[0]=="atomtype"){
      strbuf.str(args[1]); strbuf >> ts; strbuf.clear();
      strbuf.str(args[2]); strbuf >> tl; strbuf.clear();
      elem.atomtype(ts) = tl;
    }
    */
    else if (args[0]=="mass"){
      strbuf.str(args[1]); strbuf >> ts; strbuf.clear();
      strbuf.str(args[2]); strbuf >> td; strbuf.clear();
      elem.mass(ts) = td;
    }
    else if (args[0]=="lat"){
      /* Specify reference lattice. TWO species expected. */
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[3]); strbuf >> ts2; strbuf.clear();
      elem.reflat(ts1) = ts2;
    }
    else if (args[0]=="a"){
      /* Specify reference lattice parameter a. TWO species expected. */
      strbuf.str(args[1]); strbuf >> ts; strbuf.clear();
      strbuf.str(args[3]); strbuf >> td; strbuf.clear();
      elem.reflat_a(ts) = td;
      elem.reflat_b(ts) = td; // lattices by default are cubic, unless b,c, etc specified later
      elem.reflat_c(ts) = td; // lattices by default are cubic, unless b,c, etc specified later
    }
    else if (args[0]=="b"){
      /* Specify reference lattice parameter b. TWO species expected. */
      strbuf.str(args[1]); strbuf >> ts; strbuf.clear();
      strbuf.str(args[3]); strbuf >> td; strbuf.clear();
      elem.reflat_b(ts) = td;
    }
    else if (args[0]=="c"){
      /* Specify reference lattice parameter c. TWO species expected. */
      strbuf.str(args[1]); strbuf >> ts; strbuf.clear();
      strbuf.str(args[3]); strbuf >> td; strbuf.clear();
      elem.reflat_c(ts) = td;
    }
    else if (args[0]=="bpa"){
      /* Specify reference lattice parameter ratio b/a. TWO species expected. */
      strbuf.str(args[1]); strbuf >> ts; strbuf.clear();
      strbuf.str(args[3]); strbuf >> td; strbuf.clear();
      elem.reflat_bpa(ts) = td;
    }
    else if (args[0]=="cpa"){
      /* Specify reference lattice parameter ratio c/a. TWO species expected. */
      strbuf.str(args[1]); strbuf >> ts; strbuf.clear();
      strbuf.str(args[3]); strbuf >> td; strbuf.clear();
      elem.reflat_cpa(ts) = td;
    }


    if (!fp) break;
  }
  fp.close();
  fp.clear();








  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information about interactions.
     The rest of the data will be obtained from subsequent readings of the file.

     ############################################################################
     ############################################################################ */

  std::cout << "Pass 2 ... " << std::endl;


  iacs.init( elem );

  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;


    if (args[0]=="iac"){
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> ts;  strbuf.clear();

      iacs.name(ts1, ts2) = ts;

    }

    if (! fp) break; 
  }
  fp.close();
  fp.clear();



  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Some debugging:
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for (int i=0; i<elem.nelem(); i++){
    for (int j=i; j<elem.nelem(); j++){
      std::string ts1 = elem.idx2name(i);
      std::string ts2 = elem.idx2name(j);
      if (iacs.name(ts1, ts2)=="none"){
	aborterror("ERROR: Interaction type for species " + ts1 + "-" + ts2 +
		   " has not been specified. Exiting.");
      }
    }
  }




  /* ############################################################################
     ############################################################################

     Initialize the PotentialInformation object with elements and interactions

     ############################################################################
     ############################################################################ */

  //  init(el, ia);
  
  std::cout << "Initializing interactions ..." << std::endl;

  mbasepot.resize(       elem.nelem(), elem.nelem());
  mbasepot_vecidx.resize(elem.nelem(), elem.nelem());

  int ipair_EAM    = 0;
  int ipair_ABOP   = 0;
  int npairs = 0;

  // Establish vectors for basic potentials:
  for (int i=0; i<elem.nelem(); ++i){
    for (int j=i; j<elem.nelem(); ++j){

      std::string si = elem.idx2name(i);
      std::string sj = elem.idx2name(j);

      mbasepot.elem(i,j) = "none";

      if (iacs.name(si,sj)=="EAM-s"){
	mbasepot.elem(i,j) = "EAM";
	mbasepot_vecidx.elem(i,j) = ipair_EAM++;
	pot_EAM.resize( ipair_EAM );
      }
      else if (iacs.name(si,sj)=="EAM-p"){
	mbasepot.elem(i,j) = "EAM";
	mbasepot_vecidx.elem(i,j) = ipair_EAM++;
	pot_EAM.resize( ipair_EAM );
      }
      else if (iacs.name(si,sj)=="EAM-d"){
	mbasepot.elem(i,j) = "EAM";
	mbasepot_vecidx.elem(i,j) = ipair_EAM++;
	pot_EAM.resize( ipair_EAM );
      }
      else if (iacs.name(si,sj)=="EAM-sd"){
	mbasepot.elem(i,j) = "EAM";
	mbasepot_vecidx.elem(i,j) = ipair_EAM++;
	pot_EAM.resize( ipair_EAM );
      }
      else if (iacs.name(si,sj)=="EAM-pd"){
	mbasepot.elem(i,j) = "EAM";
	mbasepot_vecidx.elem(i,j) = ipair_EAM++;
	pot_EAM.resize( ipair_EAM );
      }
      else if (iacs.name(si,sj)=="EAM-spd"){
	mbasepot.elem(i,j) = "EAM";
	mbasepot_vecidx.elem(i,j) = ipair_EAM++;
	pot_EAM.resize( ipair_EAM );
      }
      else if (iacs.name(si,sj)=="ABOP"){
	mbasepot.elem(i,j) = "ABOP";
	mbasepot_vecidx.elem(i,j) = ipair_ABOP++;
	pot_ABOP.resize( ipair_ABOP );
      }

      if (mbasepot.elem(i,j) == "none"){
	aborterror("ERROR: No interaction specified for "
		   + elem.idx2name(i) + "-" + elem.idx2name(j) + ". Exiting.");
      }

      npairs++;
    }
  }

  // Reppot:
  // ------------------------------------------------------
  muse_reppot.resize(npairs);
  // default: reppot not used
  for (int i=0; i<npairs; ++i) muse_reppot[i] = false;
  mreppot_vecidx.resize(elem.nelem(), elem.nelem());


  // ABOP alpha, omega, 2mu matrices:
  // ------------------------------------------------------
  int n = elem.nelem();
  abop_alpha.resize(n,n,n);
  mabop_omega.resize(n,n,n); // private
  use_abop_alpha.resize(n,n,n);
  use_abop_omega.resize(n,n,n);

  abop_2mu.resize(n,n);
  use_abop_2mu.resize(n,n);

  /*
    ----------------------------------------------------------------------
    Notes:
    ----------------------------------------------------------------------
    ABOP, full form: use abop_alpha
                   : use abop_omega as indepependent parameter
		   : do NOT use abop_2mu

    ABOP, 2mu form : do NOT use abop_alpha
                   : do NOT use abop_omega
		   : use abop_2mu

    Brenner        : use abop_alpha
                   : use abop_omega expression explicitly, given by
		       omega_ijk = exp(-alpha_ijk*(r0_ij - r0_ik))
		   : do NOT use abop_2mu
    ----------------------------------------------------------------------
    There is no default version.
    ----------------------------------------------------------------------
  */


  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      for (int k=0; k<n; k++){
	// default: alpha = 0
	abop_alpha.elem(i,j,k)          = 0.0;
	use_abop_alpha.elem(i,j,k)      = false;
	use_abop_omega.elem(i,j,k)      = false;
	// default: omega = 1
	mabop_omega.elem(i,j,k)         = 1.0;
      }
    }
  }
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      abop_2mu.elem(i,j)     = 0.0;
      use_abop_2mu.elem(i,j) = false;
    }
  }






  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information about reppots.

     ############################################################################
     ############################################################################ */

  std::cout << "Pass 3 ..." << std::endl;


  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;

    // Use reppot?
    if (args[0]=="use_rep_core"){
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> ts;  strbuf.clear();

      use_reppot(ts1, ts2) = get_boolean_choice(ts);

    }



    if (! fp) break; 
  }
  fp.close();
  fp.clear();



  // Establish the vector of reppots based on read-in information:
  reppot_finalize();

  // Read the reppot potential data:
  read_reppot();

  // Read EAM data:
  read_eampot();





  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information about analytical parametrizations.

     ############################################################################
     ############################################################################ */

  std::cout << "Pass 4 ... " << std::endl;


  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;



    // Parameters of potentials
    if (args[0]=="potpar"){
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> ts;  strbuf.clear();
      
      potname = basepot(ts1, ts2);
      ivec    = basepot_vecidx(ts1, ts2);
      i1 = elem.name2idx(ts1);
      i2 = elem.name2idx(ts2);
      if (i1>i2){
	int k = i2;
	i2 = i1;
	i1 = k;
      }
      s1 = elem.idx2name(i1);
      s2 = elem.idx2name(i2);



      // *******************************************************************
      // ABOP
      // *******************************************************************
      if (potname=="ABOP"){
	pot_ABOP[ivec].elemname1 = s1;
	pot_ABOP[ivec].elemname2 = s2;

	if (ts=="option"){
	  // ********************************************************
	  strbuf.str(args[4]); strbuf >> tso;  strbuf.clear();
	  strbuf.str(args[5]); strbuf >> tsoi; strbuf.clear();

	  if (tso=="rcut_fun"){
	    if (tsoi=="tersoff"){
	      pot_ABOP[ivec].rcut_fun = "tersoff";
	    }
	    else if (tsoi=="perriot"){
	      pot_ABOP[ivec].rcut_fun = "perriot";
	    }
	  }

	  /*
	  if (tso=="albepot")
	    pot_ABOP[ivec].albepot = get_boolean_choice(tsoi);
	  else if (tso=="exponentone")
	    pot_ABOP[ivec].exponentone = get_boolean_choice(tsoi);
	  */

	}
	else {
	  // ********************************************************
	  strbuf.str(args[4]); strbuf >> td;  strbuf.clear();

	  int ip = pot_ABOP[ivec].parname2idx(ts);
	  if (ip>=0 && ip<=pot_ABOP[ivec].maxindex)
	    pot_ABOP[ivec].parval[ip] = td;

	  if      (ts=="D" || ts=="R"){
	    pot_ABOP[ivec].rcut_fun = "tersoff";
	  }
	  else if (ts=="pm" || ts=="pn" || ts=="prcut" ||
		   ts=="prmin" || ts =="prmax" ){
	    pot_ABOP[ivec].rcut_fun = "perriot";
	  }

	  
	}


      }
      
    }
    else if (args[0]=="abop_alpha"){
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> ts3; strbuf.clear();
      strbuf.str(args[4]); strbuf >> td; strbuf.clear();

      i1 = elem.name2idx(ts1);
      i2 = elem.name2idx(ts2);
      i3 = elem.name2idx(ts3);
      abop_alpha.elem(i1, i2, i3) = td;
      use_abop_alpha.elem(i1, i2, i3) = true;
    }
    else if (args[0]=="abop_omega"){
      // **************************************************************************
      // ABOP omega's which are specified are taken as independent parameters.
      // **************************************************************************
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> ts3; strbuf.clear();
      strbuf.str(args[4]); strbuf >> td; strbuf.clear();

      i1 = elem.name2idx(ts1);
      i2 = elem.name2idx(ts2);
      i3 = elem.name2idx(ts3);
      mabop_omega.elem(i1, i2, i3) = td;
      use_abop_omega.elem(i1, i2, i3) = true;
    }
    else if (args[0]=="abop_2mu"){
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> td; strbuf.clear();

      i1 = elem.name2idx(ts1);
      i2 = elem.name2idx(ts2);
      abop_2mu.elem(i1, i2)    = td;
      use_abop_2mu.elem(i1,i2) = true; // ABOP, 2mu form
    }



    if (! fp) break;
  }







  /* ############################################################################
     ############################################################################

     All info has been read in.

     ############################################################################
     ############################################################################ */








  /* ############################################################################
     ############################################################################

     Debug read information

     ############################################################################
     ############################################################################ */


  std::cout << "#########################################################################" << std::endl;
  std::cout << "Debugging ABOP cutoff parameters ..." << std::endl;
  std::cout << "#########################################################################" << std::endl;

  for (int i=0; i<elem.nelem(); i++){
    for (int j=i; j<elem.nelem(); j++){
      std::string s1 = elem.idx2name(i);
      std::string s2 = elem.idx2name(j);

      if (basepot(s1,s2)=="none")
	aborterror("Interaction type for " + s1 + "-" + s2 + " is unknown. Exiting.");


      if (basepot(s1,s2)=="ABOP"){
	ivec = basepot_vecidx(s1,s2);

	if      (pot_ABOP[ivec].rcut_fun=="tersoff"){
	  std::cout << "Note: Using Tersoff cutoff for " << s1 << "-" << s2 << "." << std::endl;
	}
	else if (pot_ABOP[ivec].rcut_fun=="perriot"){
	  std::cout << "Note: Using Perriot cutoff for " << s1 << "-" << s2 << "." << std::endl;
	}
	else
	  aborterror("ERROR: No cutoff given for " + s1 + "-" + s2 + ". Exiting.");
      }


    }
  }


  std::cout << "#########################################################################" << std::endl;
  std::cout << "Debugging ABOP alpha/omega/2mu parameters ..." << std::endl;
  std::cout << "#########################################################################" << std::endl;

  for (int i1=0; i1<elem.nelem(); i1++){
    for (int i2=0; i2<elem.nelem(); i2++){
      for (int i3=0; i3<elem.nelem(); i3++){
	std::string s1 = elem.idx2name(i1);
	std::string s2 = elem.idx2name(i2);
	std::string s3 = elem.idx2name(i3);
	
	std::cout << "Using abop_alpha(" << s1 << "," << s2 << "," << s3 << ")?: " << use_abop_alpha.elem(i1, i2, i3) << std::endl;
	std::cout << "Using abop_omega(" << s1 << "," << s2 << "," << s3 << ")?: " << use_abop_omega.elem(i1, i2, i3) << std::endl;
      }
      std::cout << "Using abop_2mu(" << s1 << "," << s2 << ")?: " << use_abop_2mu.elem(i1, i2) << std::endl;
    }
  }

  std::cout << "#########################################################################" << std::endl;
  std::cout << "Debugging ABOP alpha/omega/2mu parameter clashes ..." << std::endl;
  std::cout << "#########################################################################" << std::endl;


  for (int i1=0; i1<elem.nelem(); i1++){
    for (int i3=0; i3<elem.nelem(); i3++){

      int n_abop_alpha=0;
      int n_abop_omega=0;
      int n_abop_alpha_abop_omega=0;
      int n_abop_alpha_Brenner_omega=0;
      for (int i2=0; i2<elem.nelem(); i2++){
	std::string s1 = elem.idx2name(i1);
	std::string s2 = elem.idx2name(i2);
	std::string s3 = elem.idx2name(i3);

	if (use_abop_alpha.elem(i1,i2,i3)==true) n_abop_alpha++;
	if (use_abop_omega.elem(i1,i2,i3)==true) n_abop_omega++;

	if      (use_abop_alpha.elem(i1,i2,i3)==true && use_abop_omega.elem(i1,i2,i3)==true){
	  n_abop_alpha_abop_omega++;
	  std::cout << "NOTE: Using ABOP alpha and ABOP omega for combination (ijk) "
	       << s1 << "-" << s2 << "-" << s3 << std::endl;
	}
	else if (use_abop_alpha.elem(i1,i2,i3)==true && use_abop_omega.elem(i1,i2,i3)==false){
	  n_abop_alpha_Brenner_omega++;
	  std::cout << "NOTE: Using ABOP alpha and *Brenner omega* for combination (ijk) "
	       << s1 << "-" << s2 << "-" << s3 << std::endl;
	}
      }

      if (n_abop_alpha_abop_omega==0 && n_abop_alpha_Brenner_omega==0 &&
	  use_abop_2mu.elem(i1,i3)==true){
	std::cout << "NOTE: Using ABOP 2mu(i,k) for combination "
	     << s1 << "-" << s2 << std::endl;
      }

      if (use_abop_2mu.elem(i1,i3)==true && (n_abop_alpha>0 || n_abop_omega>0)){
	std::cout << "*** ERROR *** There is a clash between usage of alpha/omega and 2mu parameters! Exiting." << std::endl;
	exit(EXIT_FAILURE);
      }


    }
  }


  std::cout << "#########################################################################" << std::endl;
  std::cout << "Read-in of general information about potentials completed." << std::endl;
  std::cout << "#########################################################################" << std::endl;


  return;
}
















// ###############################################################################
// ###############################################################################
// ###############################################################################
// ###############################################################################




void PotentialInformationFit::read_info_fit(string filename){
  std::ifstream fp;
  std::ofstream fpo;
  std::string line, ts, ts1, ts2, ts3, tsm, tsi, s1, s2;
  std::vector<std::string> args;
  std::istringstream strbuf;
  double td;
  int ns, i1,i2,i3, ivec;
  std::string potname;





  int npairs = 0;

  std::cout << "Initializing some info about fittable potentials ..." << std::endl;


  for (int i=0; i<elem.nelem(); ++i){
    for (int j=i; j<elem.nelem(); ++j){
      npairs++;
    }
  }

  // Fitting indicator:
  // ------------------------------------------------------
  mfit.resize(npairs);
  // default: fixed potential, not fittable
  for (int i=0; i<npairs; ++i) mfit[i] = false;


  // Initialize ABOP omega and alpha matrices for fittable ABOP parametrizations:
  int n = elem.nelem();

  abop_alpha_parmin.resize(n,n,n);
  abop_alpha_parmax.resize(n,n,n);
  abop_alpha_partype.resize(n,n,n);

  abop_omega_parmin.resize(n,n,n);
  abop_omega_parmax.resize(n,n,n);
  abop_omega_partype.resize(n,n,n);

  abop_2mu_parmin.resize(n,n);
  abop_2mu_parmax.resize(n,n);
  abop_2mu_partype.resize(n,n);

  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      for (int k=0; k<n; k++){
	// default: alpha_ijk is fixed
	abop_alpha_parmin.elem(i,j,k) = 1;
	abop_alpha_parmax.elem(i,j,k) = 1;
	abop_alpha_partype.elem(i,j,k) = PARAM_FIXED;

	// default: omega_ijk is fixed
	abop_omega_parmin.elem(i,j,k) = 1;
	abop_omega_parmax.elem(i,j,k) = 1;
	abop_omega_partype.elem(i,j,k) = PARAM_FIXED;
      }
    }
  }
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      abop_2mu_parmin.elem(i,j) = 1;
      abop_2mu_parmax.elem(i,j) = 1;
      abop_2mu_partype.elem(i,j) = PARAM_FIXED;
    }
  }



  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information if an interaction is to be fitted,
     if reppots are used.
     The rest of the data will be obtained from subsequent readings of the file.

     ############################################################################
     ############################################################################ */


  std::cout << "Reading information about fittable potentials ..." << std::endl;

  std::cout << "Pass 1 ... " << std::endl;

  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;


    // Fit?
    if (args[0]=="fit"){
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> ts;  strbuf.clear();

      is_fittable(ts1, ts2) = get_boolean_choice(ts);
      
      if ( is_fittable(ts1, ts2) ){

	// Allocate space for parameter limits:
	potname = basepot(ts1, ts2);
	ivec    = basepot_vecidx(ts1, ts2);
	if (potname=="ABOP") pot_ABOP[ivec].init_lims();
      }

      
    }


    if (! fp) break; 
  }
  fp.close();
  fp.clear();




  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Some delayed debugging for same-element interaction:
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (int i=0; i<elem.nelem(); ++i){
    std::string ts = elem.idx2name(i);

    if ( is_fittable(ts,ts) ) continue;

    if (elem.reflat(ts)=="skip"){
      std::cout << "WARNING: Skipping this reference lattice. Hopefully you know what you are doing ..."
		<< std::endl;
      continue;
    }


    if (elem.reflat(ts)=="none"){
      aborterror("ERROR: Reference lattice for species " + ts +
		 " has not been specified. Exiting.");
    }
    if (elem.reflat_a(ts)<0)
      aborterror("Error: Lattice parameter a for reference lattice " + ts + " " +
		 "is missing. Exiting.");
    
    if (elem.reflat_bpa(ts)>0) elem.reflat_b(ts) = elem.reflat_bpa(ts) * elem.reflat_a(ts);
    if (elem.reflat_cpa(ts)>0) elem.reflat_c(ts) = elem.reflat_cpa(ts) * elem.reflat_a(ts);

    if (elem.reflat_bpa(ts)<0) elem.reflat_bpa(ts) = elem.reflat_b(ts)/elem.reflat_a(ts);
    if (elem.reflat_cpa(ts)<0) elem.reflat_cpa(ts) = elem.reflat_c(ts)/elem.reflat_a(ts);
  }






  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information about analytical parametrizations.

     ############################################################################
     ############################################################################ */

  std::cout << "Pass 2 ... " << std::endl;

  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;


    // Parameters of potentials that are to be fitted
    // Limits:
    if (args[0]=="min" || args[0]=="max"){

      strbuf.str(args[1]); strbuf >> tsi; strbuf.clear();

      // Parameter:
      if (tsi=="potpar"){
	// --------------------------------------------------------------------------
	strbuf.str(args[2]); strbuf >> ts1; strbuf.clear();
	strbuf.str(args[3]); strbuf >> ts2; strbuf.clear();
	strbuf.str(args[4]); strbuf >> ts;  strbuf.clear();
	strbuf.str(args[5]); strbuf >> td;  strbuf.clear();
      
	potname = basepot(ts1, ts2);
	ivec    = basepot_vecidx(ts1, ts2);
	i1 = elem.name2idx(ts1);
	i2 = elem.name2idx(ts2);
	if (i1>i2){
	  int k = i2;
	  i2 = i1;
	  i1 = k;
	}
	s1 = elem.idx2name(i1);
	s2 = elem.idx2name(i2);

	// Do not allocate any parameter limits for potentials that
	// are not to be fitted:
	if (! is_fittable(ts1,ts2)) continue;


	// *******************************************************************
	// ABOP
	// *******************************************************************

	// Allocate lower/upper limits of fittable potentials:
	if (potname=="ABOP"){
	  int ip = pot_ABOP[ivec].parname2idx(ts);
	  
	  if (args[0]=="min"){
	    if (ip>=0 && ip<=pot_ABOP[ivec].maxindex) pot_ABOP[ivec].parmin[ip] = td;
	    
	    /*
	    if      (ts=="D0") pot_ABOP[ivec].parmin->D0 = td;
	    else if (ts=="r0") pot_ABOP[ivec].parmin->r0 = td;
	    else if (ts=="beta") pot_ABOP[ivec].parmin->beta = td;
	    else if (ts=="S") pot_ABOP[ivec].parmin->S = td;
	    else if (ts=="gamma") pot_ABOP[ivec].parmin->gamma = td;
	    else if (ts=="c") pot_ABOP[ivec].parmin->c = td;
	    else if (ts=="d") pot_ABOP[ivec].parmin->d = td;
	    else if (ts=="h") pot_ABOP[ivec].parmin->h = td;
	    else if (ts=="R") pot_ABOP[ivec].parmin->R = td;
	    else if (ts=="D") pot_ABOP[ivec].parmin->D = td;
	    */
	  }
	  else if (args[0]=="max"){
	    if (ip>=0 && ip<=pot_ABOP[ivec].maxindex) pot_ABOP[ivec].parmax[ip] = td;

	    /*
	    if      (ts=="D0") pot_ABOP[ivec].parmax->D0 = td;
	    else if (ts=="r0") pot_ABOP[ivec].parmax->r0 = td;
	    else if (ts=="beta") pot_ABOP[ivec].parmax->beta = td;
	    else if (ts=="S") pot_ABOP[ivec].parmax->S = td;
	    else if (ts=="gamma") pot_ABOP[ivec].parmax->gamma = td;
	    else if (ts=="c") pot_ABOP[ivec].parmax->c = td;
	    else if (ts=="d") pot_ABOP[ivec].parmax->d = td;
	    else if (ts=="h") pot_ABOP[ivec].parmax->h = td;
	    else if (ts=="R") pot_ABOP[ivec].parmax->R = td;
	    else if (ts=="D") pot_ABOP[ivec].parmax->D = td;
	    */
	  }

	}
      }

      else if (tsi=="abop_alpha"){
	// --------------------------------------------------------------------------
	strbuf.str(args[2]); strbuf >> ts1; strbuf.clear();
	strbuf.str(args[3]); strbuf >> ts2; strbuf.clear();
	strbuf.str(args[4]); strbuf >> ts3; strbuf.clear();
	strbuf.str(args[5]); strbuf >> td;  strbuf.clear();

	i1 = elem.name2idx(ts1);
	i2 = elem.name2idx(ts2);
	i3 = elem.name2idx(ts3);

	if (use_abop_alpha.elem(i1,i2,i3)){
	  if      (args[0]=="min") abop_alpha_parmin.elem(i1, i2, i3) = td;
	  else if (args[0]=="max") abop_alpha_parmax.elem(i1, i2, i3) = td;
	}
      }

      else if (tsi=="abop_omega"){
	strbuf.str(args[2]); strbuf >> ts1; strbuf.clear();
	strbuf.str(args[3]); strbuf >> ts2; strbuf.clear();
	strbuf.str(args[4]); strbuf >> ts3; strbuf.clear();
	strbuf.str(args[5]); strbuf >> td;  strbuf.clear();

	i1 = elem.name2idx(ts1);
	i2 = elem.name2idx(ts2);
	i3 = elem.name2idx(ts3);

	if (use_abop_omega.elem(i1,i2,i3)){
	  if      (args[0]=="min") abop_omega_parmin.elem(i1, i2, i3) = td;
	  else if (args[0]=="max") abop_omega_parmax.elem(i1, i2, i3) = td;
	}
      }

      else if (tsi=="abop_2mu"){
	strbuf.str(args[2]); strbuf >> ts1; strbuf.clear();
	strbuf.str(args[3]); strbuf >> ts2; strbuf.clear();
	strbuf.str(args[4]); strbuf >> td;  strbuf.clear();
	
	i1 = elem.name2idx(ts1);
	i2 = elem.name2idx(ts2);

	if (use_abop_2mu.elem(i1,i2)){
	  if      (args[0]=="min") abop_2mu_parmin.elem(i1, i2) = td;
	  else if (args[0]=="max") abop_2mu_parmax.elem(i1, i2) = td;
	}
      }



    }


    if (! fp) break;
  }




  /* ############################################################################
     ############################################################################

     All info has been read in.

     ############################################################################
     ############################################################################ */








  /* ############################################################################
     ############################################################################

     Debug read information

     ############################################################################
     ############################################################################ */



  std::cout << "Some debugging of fittable potentials ..." << std::endl;

  /* ############################################################################
     Debug the potential settings:
     ############################################################################ */

  for (int i=0; i<elem.nelem(); i++){
    for (int j=i; j<elem.nelem(); j++){
      std::string s1 = elem.idx2name(i);
      std::string s2 = elem.idx2name(j);

      if (basepot(s1,s2)=="none"){
	aborterror("Interaction type for " + s1 + "-" + s2 + " is unknown. Exiting.");
      }


      if (is_fittable(s1,s2)){

	if (basepot(s1,s2)=="ABOP"){
	  int iv = basepot_vecidx(s1,s2);

#if 0
	  int i1,i2,i3;
	  if (pot_ABOP[iv].use_cutoff_tersoff){
	    // Disable all other cutoff parameters:
	    // Perriot:
	    i1 = pot_ABOP[iv].parname2idx("n");
	    i2 = pot_ABOP[iv].parname2idx("m");
	    i3 = pot_ABOP[iv].parname2idx("rc");
	    pot_ABOP[iv].parmin[i1]=1; pot_ABOP[iv].parmax[i1]=1;
	    pot_ABOP[iv].parmin[i2]=1; pot_ABOP[iv].parmax[i2]=1;
	    pot_ABOP[iv].parmin[i3]=1; pot_ABOP[iv].parmax[i3]=1;
	  }
	  if (pot_ABOP[iv].use_cutoff_perriot){
	    // Disable all other cutoff parameters:
	    // Tersoff:
	    i1 = pot_ABOP[iv].parname2idx("D");
	    i2 = pot_ABOP[iv].parname2idx("R");
	    pot_ABOP[iv].parmin[i1]=1; pot_ABOP[iv].parmax[i1]=1;
	    pot_ABOP[iv].parmin[i2]=1; pot_ABOP[iv].parmax[i2]=1;
	  }
#endif


	  for (int k=0; k<pot_ABOP[iv].parval.size(); ++k){

	    set_param_type(pot_ABOP[iv].parmin[k], pot_ABOP[iv].parmax[k], pot_ABOP[iv].partype[k]);

	    limcheck("ABOP", pot_ABOP[iv].parname[k], s1 + "-" + s2,
		     pot_ABOP[iv].partype[k],
		     pot_ABOP[iv].parmin[k], pot_ABOP[iv].parmax[k],
		     pot_ABOP[iv].parval[k]);
	    

	  }

	}

      }

    }
  }
  

  for (int i=0; i<elem.nelem(); i++){
    for (int j=0; j<elem.nelem(); j++){
      for (int k=0; k<elem.nelem(); k++){

	std::string s1 = elem.idx2name(i);
	std::string s2 = elem.idx2name(j);
	std::string s3 = elem.idx2name(k);

	if (use_abop_alpha.elem(i,j,k)){

	  set_param_type(abop_alpha_parmin.elem(i,j,k), abop_alpha_parmax.elem(i,j,k), abop_alpha_partype.elem(i,j,k));

	  limcheck("ABOP", "alpha", s1 + "-" + s2 + "-" + s3,
		   abop_alpha_partype.elem(i,j,k),
		   abop_alpha_parmin.elem(i,j,k),
		   abop_alpha_parmax.elem(i,j,k),
		   abop_alpha.elem(i,j,k));
	}

      }
    }
  }



  for (int i=0; i<elem.nelem(); i++){
    for (int j=0; j<elem.nelem(); j++){
      for (int k=0; k<elem.nelem(); k++){

	std::string s1 = elem.idx2name(i);
	std::string s2 = elem.idx2name(j);
	std::string s3 = elem.idx2name(k);

	if (use_abop_omega.elem(i,j,k)){

	  set_param_type(abop_omega_parmin.elem(i,j,k), abop_omega_parmax.elem(i,j,k), abop_omega_partype.elem(i,j,k));

	  limcheck("ABOP", "omega", s1 + "-" + s2 + "-" + s3,
		   abop_alpha_partype.elem(i,j,k),
		   abop_alpha_parmin.elem(i,j,k),
		   abop_alpha_parmax.elem(i,j,k),
		   abop_alpha.elem(i,j,k));
	}


      }
    }
  }



  for (int i=0; i<elem.nelem(); i++){
    for (int j=0; j<elem.nelem(); j++){

      std::string s1 = elem.idx2name(i);
      std::string s2 = elem.idx2name(j);

      if (use_abop_2mu.elem(i,j)){

	set_param_type(abop_2mu_parmin.elem(i,j), abop_2mu_parmax.elem(i,j), abop_2mu_partype.elem(i,j));

	limcheck("ABOP", "2mu", s1 + "-" + s2,
		   abop_2mu_partype.elem(i,j),
		   abop_2mu_parmin.elem(i,j),
		   abop_2mu_parmax.elem(i,j),
		   abop_2mu.elem(i,j));
      }

    }
  }


  std::cout << "##############################################################################" << std::endl;

  for (int i=0; i<elem.nelem(); i++){
    for (int j=0; j<elem.nelem(); j++){
      for (int k=0; k<elem.nelem(); k++){
	std::string s1 = elem.idx2name(i);
	std::string s2 = elem.idx2name(j);
	std::string s3 = elem.idx2name(k);

	std::cout << "ABOP alpha(" << s1 << "," << s2 << "," << s3 << "): ";
	if (use_abop_alpha.elem(i,j,k)) std::cout << "used";
	else std::cout << "NOT used";
	std::cout << std::endl;
      }
    }
  }
  for (int i=0; i<elem.nelem(); i++){
    for (int j=0; j<elem.nelem(); j++){
      for (int k=0; k<elem.nelem(); k++){
	std::string s1 = elem.idx2name(i);
	std::string s2 = elem.idx2name(j);
	std::string s3 = elem.idx2name(k);

	std::cout << "ABOP omega(" << s1 << "," << s2 << "," << s3 << "): ";
	if (use_abop_omega.elem(i,j,k)) std::cout << "used";
	else std::cout << "NOT used";
	std::cout << std::endl;
      }
    }
  }
  for (int i=0; i<elem.nelem(); i++){
    for (int j=0; j<elem.nelem(); j++){
      std::string s1 = elem.idx2name(i);
      std::string s2 = elem.idx2name(j);

      std::cout << "ABOP 2mu(" << s1 << "," << s2 << "): ";
      if (use_abop_2mu.elem(i,j)) std::cout << "used";
      else std::cout << "NOT used";
      std::cout << std::endl;
    }
  }

  std::cout << "##############################################################################" << std::endl;
  std::cout << "Read-in of information about fittable potentials completed." << std::endl;
  std::cout << "##############################################################################" << std::endl;


}





void PotentialInformationFit::limcheck(const std::string        & pot,
				       const std::string        & parname,
				       const std::string        & elems,
				       const parametertype & partype,
				       const double        & parmin,
				       const double        & parmax,
				       const double        & parval
				       ){
  double twoeps = std::numeric_limits<double>::epsilon();
  std::string intro = "ERROR: Potential " + pot + ": " + parname + ": "
    + elems + ": ";

  if (partype==PARAM_FREE_WITH_LIMITS){
    if (parmin > parmax)
      aborterror(intro
		 + "Lower parameter limit " + tostring(parmin)
		 + " is larger than upper limit " + tostring(parmax)
		 + ". Exiting.");
    if (fp_are_equal(parval, parmin, twoeps))
      aborterror(intro
		 + "Parameter value " + tostring(parval)
		 + " is too close to lower limit " + tostring(parmin)
		 + ". Exiting.");
    if (fp_are_equal(parval, parmax, twoeps))
      aborterror(intro
		 + "Parameter value " + tostring(parval)
		 + " is too close to upper limit " + tostring(parmax)
		 + ". Exiting.");
    if (parval < parmin)
      aborterror(intro
		 + "Parameter value " + tostring(parval)
		 + " is smaller than lower limit " + tostring(parmax)
		 + ". Exiting.");
    if (parval > parmax)
      aborterror(intro
		 + "Parameter value " + tostring(parval)
		 + " is larger than upper limit " + tostring(parmax)
		 + ". Exiting.");
  }
}



