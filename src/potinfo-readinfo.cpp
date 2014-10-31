



#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

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
#include "exiterrors.hpp"
#include "potinfo.hpp"
#include "param.hpp"
#include "helpfuns.hpp"

using namespace std;
using namespace utils;
using namespace exiterrors;
using boost::format;




// #################################################################################
// #################################################################################
// #################################################################################
// #################################################################################




void PotentialInformation::read_info(string filename){
  ifstream fp;
  ofstream fpo;
  string line, ts, ts1, ts2, ts3, tsm, tsi, s1, s2, tso, tsoi, ts0;
  vector<string> args;
  istringstream strbuf;
  double td;
  int tl, i1,i2,i3, ivec;
  string potname;
  int ns;

  cout << "Reading general information about potentials ..." << endl;


  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information about elements only.
     The rest of the data will be obtained from subsequent readings of the file.

     ############################################################################
     ############################################################################ */



  cout << "Pass 1 ... " << endl;

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
      cout << args[i] << endl;
    */

    if (args[0]=="elem"){
      strbuf.str(args[2]); strbuf >> ts; strbuf.clear();
      //cout << ts << endl;
      elem.add_elem(ts);
    }
    else if (args[0]=="atomtype"){
      strbuf.str(args[1]); strbuf >> ts; strbuf.clear();
      strbuf.str(args[2]); strbuf >> tl; strbuf.clear();
      elem.atomtype(ts) = tl;
    }
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

  cout << "Pass 2 ... " << endl;


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
      string ts1 = elem.idx2name(i);
      string ts2 = elem.idx2name(j);
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
  
  cout << "Initializing interactions ..." << endl;

  mbasepot.resize(       elem.nelem(), elem.nelem());
  mbasepot_vecidx.resize(elem.nelem(), elem.nelem());

  int ipair_EAM    = 0;
  int ipair_ABOP   = 0;
  int npairs = 0;

  // Establish vectors for basic potentials:
  for (int i=0; i<elem.nelem(); ++i){
    for (int j=i; j<elem.nelem(); ++j){

      string si = elem.idx2name(i);
      string sj = elem.idx2name(j);

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


  // ABOP omega and alpha matrices:
  // ------------------------------------------------------
  int n = elem.nelem();
  abop_alpha.resize(n,n,n);
  abop_omega_is_free.resize(n,n,n);
  mabop_omega.resize(n,n,n);

  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      for (int k=0; k<n; k++){
	// default: alpha = 0
	abop_alpha.elem(i,j,k)     = 0.0;
	// default: omega is independent, i.e. is not given by the form
	//          omega_ijk = exp(-alpha_ijk*(r0_ij - r0_ik))
	// note: special syntax needed in input file to make omega_ijk dependent
	abop_omega_is_free.elem(i,j,k) = true;
	// default: omega = 1
	mabop_omega.elem(i,j,k)    = 1.0;
      }
    }
  }






  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information about reppots.

     ############################################################################
     ############################################################################ */

  cout << "Pass 3 ..." << endl;


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

  cout << "Pass 4 ... " << endl;


  fp.open(filename.c_str());
  if (!fp)
    aborterror("Error: Could not find file " + filename + ". Exiting.");

  while (true){
    utils::get_line(fp, line);
    ns = utils::get_substrings( line, args, "\t :,()[]=");
    if (ns==0 && fp) continue;


    // Reppot parameters
    if (args[0]=="bfermi"){
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> td;  strbuf.clear();

      if (use_reppot(ts1, ts2))
	pot_Reppot[ reppot_vecidx(ts1, ts2) ].bfermi = td;
    }
    else if (args[0]=="rfermi"){
      strbuf.str(args[1]); strbuf >> ts1; strbuf.clear();
      strbuf.str(args[2]); strbuf >> ts2; strbuf.clear();
      strbuf.str(args[3]); strbuf >> td;  strbuf.clear();

      if (use_reppot(ts1, ts2))
	pot_Reppot[ reppot_vecidx(ts1, ts2) ].rfermi = td;
    }

    // Parameters of potentials
    else if (args[0]=="potpar"){
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
	  strbuf.str(args[4]); strbuf >> tsoi; strbuf.clear();

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
	  if (ip>=0 && ip<10)
	    pot_ABOP[ivec].parval[ip] = td;

	  /*
	  if      (ts=="D0") pot_ABOP[ivec].D0 = td;
	  else if (ts=="r0") pot_ABOP[ivec].r0 = td;
	  else if (ts=="beta") pot_ABOP[ivec].beta = td;
	  else if (ts=="S") pot_ABOP[ivec].S = td;
	  else if (ts=="gamma") pot_ABOP[ivec].gamma = td;
	  else if (ts=="c") pot_ABOP[ivec].c = td;
	  else if (ts=="d") pot_ABOP[ivec].d = td;
	  else if (ts=="h") pot_ABOP[ivec].h = td;
	  else if (ts=="R") pot_ABOP[ivec].R = td;
	  else if (ts=="D") pot_ABOP[ivec].D = td;
	  */
	  
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
      abop_omega_is_free.elem(i1, i2, i3) = true;
    }

    else if (args[0]=="rule"){

      // # Handle entries like:
      // rule : abop_omega( Fe, Cr, C ) : use_Brenner_form
      // # All omega permutations of Fe, Cr, C will be specified to be dependent
      // # parameters, i.e. omega_ijk = exp(-alpha_ijk*(r0_ij - r0_ik))
      if (args[1]=="abop_omega"){
	strbuf.str(args[2]); strbuf >> ts1; strbuf.clear();
	strbuf.str(args[3]); strbuf >> ts2; strbuf.clear();
	strbuf.str(args[4]); strbuf >> ts3; strbuf.clear();

	i1 = elem.name2idx(ts1);
	i2 = elem.name2idx(ts2);
	i3 = elem.name2idx(ts3);

	if (args[5]=="use_Brenner_form"){
	  for (int i=0; i<elem.nelem(); ++i)
	    for (int j=0; j<elem.nelem(); ++j)
	      for (int k=0; k<elem.nelem(); ++k)
		abop_omega_is_free.elem(i, j, k) = false;
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


  for (int i=0; i<elem.nelem(); i++){
    for (int j=i; j<elem.nelem(); j++){
      string s1 = elem.idx2name(i);
      string s2 = elem.idx2name(j);

      if (basepot(s1,s2)=="none")
	aborterror("Interaction type for " + s1 + "-" + s2 + " is unknown. Exiting.");


    }
  }

  

  cout << "Read-in of general information about potentials completed." << endl;

  return;
}
















// ###############################################################################
// ###############################################################################
// ###############################################################################
// ###############################################################################




void PotentialInformationFit::read_info_fit(string filename){
  ifstream fp;
  ofstream fpo;
  string line, ts, ts1, ts2, ts3, tsm, tsi, s1, s2;
  vector<string> args;
  istringstream strbuf;
  double td;
  int ns, i1,i2,i3, ivec;
  string potname;





  int npairs = 0;

  cout << "Initializing some info about fittable potentials ..." << endl;


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



  /* ############################################################################
     ############################################################################

     Read the potinfo file to get information if an interaction is to be fitted,
     if reppots are used.
     The rest of the data will be obtained from subsequent readings of the file.

     ############################################################################
     ############################################################################ */


  cout << "Reading information about fittable potentials ..." << endl;

  cout << "Pass 1 ... " << endl;

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
    string ts = elem.idx2name(i);

    if ( is_fittable(ts,ts) ) continue;

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

  cout << "Pass 2 ... " << endl;

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
	    if (ip>=0 && ip<10) pot_ABOP[ivec].parmin[ip] = td;
	    
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
	    if (ip>=0 && ip<10) pot_ABOP[ivec].parmax[ip] = td;

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

	if      (args[0]=="min") abop_alpha_parmin.elem(i1, i2, i3) = td;
	else if (args[0]=="max") abop_alpha_parmax.elem(i1, i2, i3) = td;
      }

      else if (tsi=="abop_omega"){
	strbuf.str(args[2]); strbuf >> ts1; strbuf.clear();
	strbuf.str(args[3]); strbuf >> ts2; strbuf.clear();
	strbuf.str(args[4]); strbuf >> ts3; strbuf.clear();
	strbuf.str(args[5]); strbuf >> td;  strbuf.clear();

	i1 = elem.name2idx(ts1);
	i2 = elem.name2idx(ts2);
	i3 = elem.name2idx(ts3);

	if (abop_omega_is_free.elem(i1,i2,i3)){
	  if      (args[0]=="min") abop_omega_parmin.elem(i1, i2, i3) = td;
	  else if (args[0]=="max") abop_omega_parmax.elem(i1, i2, i3) = td;
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



  cout << "Some debugging of fittable potentials ..." << endl;

  /* ############################################################################
     Debug the potential settings:
     ############################################################################ */

  for (int i=0; i<elem.nelem(); i++){
    for (int j=i; j<elem.nelem(); j++){
      string s1 = elem.idx2name(i);
      string s2 = elem.idx2name(j);

      if (basepot(s1,s2)=="none"){
	aborterror("Interaction type for " + s1 + "-" + s2 + " is unknown. Exiting.");
      }


      if (is_fittable(s1,s2)){

	if (basepot(s1,s2)=="ABOP"){
	  int iv = basepot_vecidx(s1,s2);

	  for (int k=0; k<pot_ABOP[iv].parval.size(); ++k){
	    
	    if (pot_ABOP[iv].parmin[k] > pot_ABOP[iv].parmax[k])
	      limerr1("ABOP", s1, s2, pot_ABOP[iv].parname[k]);

	    set_param_type(pot_ABOP[iv].parmin[k], pot_ABOP[iv].parmax[k], pot_ABOP[iv].partype[k]);

	    if (pot_ABOP[iv].partype[k] == PARAM_FREE_WITH_LIMITS
		&&
		(pot_ABOP[iv].parmin[k] > pot_ABOP[iv].parval[k]))
	      limerr2("ABOP", s1, s2, pot_ABOP[iv].parname[k]);

	    if (pot_ABOP[iv].partype[k] == PARAM_FREE_WITH_LIMITS
		&&
		(pot_ABOP[iv].parval[k] > pot_ABOP[iv].parmax[k]))
	      limerr3("ABOP", s1, s2, pot_ABOP[iv].parname[k]);

	  }

	    /*
	      if (pot_ABOP[iv].parmin->D0    > pot_ABOP[iv].parmax->D0)    limerr1("ABOP", s1, s2, "D0");
	      if (pot_ABOP[iv].parmin->r0    > pot_ABOP[iv].parmax->r0)    limerr1("ABOP", s1, s2, "r0");
	      if (pot_ABOP[iv].parmin->beta  > pot_ABOP[iv].parmax->beta)  limerr1("ABOP", s1, s2, "beta");
	      if (pot_ABOP[iv].parmin->S     > pot_ABOP[iv].parmax->S)     limerr1("ABOP", s1, s2, "S");
	      if (pot_ABOP[iv].parmin->gamma > pot_ABOP[iv].parmax->gamma) limerr1("ABOP", s1, s2, "gamma");
	      if (pot_ABOP[iv].parmin->c > pot_ABOP[iv].parmax->c) limerr1("ABOP", s1, s2, "c");
	      if (pot_ABOP[iv].parmin->d > pot_ABOP[iv].parmax->d) limerr1("ABOP", s1, s2, "d");
	      if (pot_ABOP[iv].parmin->h > pot_ABOP[iv].parmax->h) limerr1("ABOP", s1, s2, "h");
	      if (pot_ABOP[iv].parmin->R > pot_ABOP[iv].parmax->R) limerr1("ABOP", s1, s2, "R");
	      if (pot_ABOP[iv].parmin->D > pot_ABOP[iv].parmax->D) limerr1("ABOP", s1, s2, "D");

	      set_param_type(pot_ABOP[iv].parmin->D0,    pot_ABOP[iv].parmax->D0,    pot_ABOP[iv].partype->D0);
	      set_param_type(pot_ABOP[iv].parmin->r0,    pot_ABOP[iv].parmax->r0,    pot_ABOP[iv].partype->r0);
	      set_param_type(pot_ABOP[iv].parmin->beta,  pot_ABOP[iv].parmax->beta,  pot_ABOP[iv].partype->beta);
	      set_param_type(pot_ABOP[iv].parmin->S,     pot_ABOP[iv].parmax->S,     pot_ABOP[iv].partype->S);
	      set_param_type(pot_ABOP[iv].parmin->gamma, pot_ABOP[iv].parmax->gamma, pot_ABOP[iv].partype->gamma);
	      set_param_type(pot_ABOP[iv].parmin->c, pot_ABOP[iv].parmax->c, pot_ABOP[iv].partype->c);
	      set_param_type(pot_ABOP[iv].parmin->d, pot_ABOP[iv].parmax->d, pot_ABOP[iv].partype->d);
	      set_param_type(pot_ABOP[iv].parmin->h, pot_ABOP[iv].parmax->h, pot_ABOP[iv].partype->h);
	      set_param_type(pot_ABOP[iv].parmin->R, pot_ABOP[iv].parmax->R, pot_ABOP[iv].partype->R);
	      set_param_type(pot_ABOP[iv].parmin->D, pot_ABOP[iv].parmax->D, pot_ABOP[iv].partype->D);

	      if (pot_ABOP[ivec].partype->D0 == PARAM_FREE_WITH_LIMITS    && (pot_ABOP[iv].parmin->D0 >    pot_ABOP[iv].D0))    limerr2("ABOP", s1, s2, "D0");
	      if (pot_ABOP[ivec].partype->r0 == PARAM_FREE_WITH_LIMITS    && (pot_ABOP[iv].parmin->r0 >    pot_ABOP[iv].r0))    limerr2("ABOP", s1, s2, "r0");
	      if (pot_ABOP[ivec].partype->beta == PARAM_FREE_WITH_LIMITS  && (pot_ABOP[iv].parmin->beta >  pot_ABOP[iv].beta))  limerr2("ABOP", s1, s2, "beta");
	      if (pot_ABOP[ivec].partype->S == PARAM_FREE_WITH_LIMITS     && (pot_ABOP[iv].parmin->S >     pot_ABOP[iv].S))     limerr2("ABOP", s1, s2, "S");
	      if (pot_ABOP[ivec].partype->gamma == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].parmin->gamma > pot_ABOP[iv].gamma)) limerr2("ABOP", s1, s2, "gamma");
	      if (pot_ABOP[ivec].partype->c == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].parmin->c > pot_ABOP[iv].c)) limerr2("ABOP", s1, s2, "c");
	      if (pot_ABOP[ivec].partype->d == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].parmin->d > pot_ABOP[iv].d)) limerr2("ABOP", s1, s2, "d");
	      if (pot_ABOP[ivec].partype->h == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].parmin->h > pot_ABOP[iv].h)) limerr2("ABOP", s1, s2, "h");
	      if (pot_ABOP[ivec].partype->R == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].parmin->R > pot_ABOP[iv].R)) limerr2("ABOP", s1, s2, "D");
	      if (pot_ABOP[ivec].partype->D == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].parmin->D > pot_ABOP[iv].D)) limerr2("ABOP", s1, s2, "R");

	      if (pot_ABOP[ivec].partype->D0 == PARAM_FREE_WITH_LIMITS    && (pot_ABOP[iv].D0 >    pot_ABOP[iv].parmax->D0))   limerr3("ABOP", s1, s2, "D0");
	      if (pot_ABOP[ivec].partype->r0 == PARAM_FREE_WITH_LIMITS    && (pot_ABOP[iv].r0 >    pot_ABOP[iv].parmax->r0))   limerr3("ABOP", s1, s2, "r0");
	      if (pot_ABOP[ivec].partype->beta == PARAM_FREE_WITH_LIMITS  && (pot_ABOP[iv].beta >  pot_ABOP[iv].parmax->beta)) limerr3("ABOP", s1, s2, "beta");
	      if (pot_ABOP[ivec].partype->S == PARAM_FREE_WITH_LIMITS     && (pot_ABOP[iv].S >     pot_ABOP[iv].parmax->S))    limerr3("ABOP", s1, s2, "S");
	      if (pot_ABOP[ivec].partype->gamma == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].gamma > pot_ABOP[iv].parmax->gamma)) limerr3("ABOP", s1, s2, "gamma");
	      if (pot_ABOP[ivec].partype->c == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].c > pot_ABOP[iv].parmax->c)) limerr3("ABOP", s1, s2, "c");
	      if (pot_ABOP[ivec].partype->d == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].d > pot_ABOP[iv].parmax->d)) limerr3("ABOP", s1, s2, "d");
	      if (pot_ABOP[ivec].partype->h == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].h > pot_ABOP[iv].parmax->h)) limerr3("ABOP", s1, s2, "h");
	      if (pot_ABOP[ivec].partype->R == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].R > pot_ABOP[iv].parmax->R)) limerr3("ABOP", s1, s2, "D");
	      if (pot_ABOP[ivec].partype->D == PARAM_FREE_WITH_LIMITS && (pot_ABOP[iv].D > pot_ABOP[iv].parmax->D)) limerr3("ABOP", s1, s2, "R");
	    */


	}

      }

    }
  }
  

  for (int i=0; i<elem.nelem(); i++){
    for (int j=0; j<elem.nelem(); j++){
      for (int k=0; k<elem.nelem(); k++){

	string s1 = elem.idx2name(i);
	string s2 = elem.idx2name(j);
	string s3 = elem.idx2name(k);

	if (abop_alpha_parmin.elem(i,j,k) > abop_alpha_parmax.elem(i,j,k))
	  limerr1_abop("alpha", s1,s2,s3);

	set_param_type(abop_alpha_parmin.elem(i,j,k), abop_alpha_parmax.elem(i,j,k), abop_alpha_partype.elem(i,j,k));

	if (abop_alpha_partype.elem(i,j,k) == PARAM_FREE_WITH_LIMITS && (abop_alpha_parmin.elem(i,j,k) > abop_alpha.elem(i,j,k)))
	  limerr2_abop("alpha", s1,s2,s3);
	if (abop_alpha_partype.elem(i,j,k) == PARAM_FREE_WITH_LIMITS && (abop_alpha_parmax.elem(i,j,k) < abop_alpha.elem(i,j,k)))
	  limerr3_abop("alpha", s1,s2,s3);

      }
    }
  }



  for (int i=0; i<elem.nelem(); i++){
    for (int j=0; j<elem.nelem(); j++){
      for (int k=0; k<elem.nelem(); k++){

	string s1 = elem.idx2name(i);
	string s2 = elem.idx2name(j);
	string s3 = elem.idx2name(k);


	if (abop_omega_is_free.elem(i,j,k)){

	  if (abop_omega_parmin.elem(i,j,k) > abop_omega_parmax.elem(i,j,k))
	    limerr1_abop("omega", s1,s2,s3);

	  set_param_type(abop_omega_parmin.elem(i,j,k), abop_omega_parmax.elem(i,j,k), abop_omega_partype.elem(i,j,k));

	  if (abop_omega_partype.elem(i,j,k) == PARAM_FREE_WITH_LIMITS && (abop_omega_parmin.elem(i,j,k) > get_abop_omega(s1,s2,s3)))
	    limerr2_abop("omega", s1,s2,s3);
	  if (abop_omega_partype.elem(i,j,k) == PARAM_FREE_WITH_LIMITS && (abop_omega_parmax.elem(i,j,k) < get_abop_omega(s1,s2,s3)))
	    limerr3_abop("omega", s1,s2,s3);

	}


      }
    }
  }




  cout << "Read-in of information about fittable potentials completed." << endl;

}







void PotentialInformationFit::limerr1(string s0, string s1, string s2, string sp){
  aborterror("Error: " + s0 + ": " + s1 + "-" + s2 + ": " + sp +
	     " limits are reversed ?! Exiting.");
}
void PotentialInformationFit::limerr2(string s0, string s1, string s2, string sp){
  aborterror("Error: " + s0 + ": " + s1 + "-" + s2 + ": " + sp +
	     " value is smaller than lower limit ?! Exiting.");
}
void PotentialInformationFit::limerr3(string s0, string s1, string s2, string sp){
  aborterror("Error: " + s0 + ": " + s1 + "-" + s2 + ": " + sp +
	     " value is larger than upper limit ?! Exiting.");
}

void PotentialInformationFit::limerr1_abop(string s0, string s1, string s2, string s3){
  aborterror("Error: ABOP " + s0 + ": " + s1 + "-" + s2 + "-" + s3 + ": " +
	     " limits are reversed ?! Exiting.");
}
void PotentialInformationFit::limerr2_abop(string s0, string s1, string s2, string s3){
  aborterror("Error: ABOP " + s0 + ": " + s1 + "-" + s2 + "-" + s3 + ": " +
	     " value is smaller than lower limit ?! Exiting.");
}
void PotentialInformationFit::limerr3_abop(string s0, string s1, string s2, string s3){
  aborterror("Error: ABOP " + s0 + ": " + s1 + "-" + s2 + "-" + s3 + ": " +
	     " value is larger than upper limit ?! Exiting.");
}






