


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
#include "utils-errors.hpp"

#include "potinfo.hpp"
#include "helpfuns.hpp"
#include "utils-streamio.hpp"


using namespace utils;
using boost::format;


void PotentialInformation::read_eampot(void){


  std::cout << "Checking for EAM potentials to read in ..." << std::endl;
  int nread=0;

  for (int i=0; i<elem.nelem(); ++i){
    for (int j=i; j<elem.nelem(); ++j){

      int i1 = i;
      int i2 = j;
      
      int k;
      if (i1>i2){
	k = i2;
	i2 = i1;
	i1 = k;
      }

      std::string s1 = elem.idx2name(i1);
      std::string s2 = elem.idx2name(i2);

      if (basepot(s1,s2) != "EAM") continue;

      int ivec = basepot_vecidx(s1, s2);


      std::string filename;

      double r, dr;
      double rho, drho;

      std::ifstream fp;
      std::ofstream fpo;
      std::string line;
      std::vector<std::string> args;
      std::istringstream strbuf;

  
      pot_EAM[ivec].elemname1 = elem.idx2name(i1);
      pot_EAM[ivec].elemname2 = elem.idx2name(i2);

      int eam = -1;

      /* Initialize. */
      pot_EAM[ivec].Nr     = 0;
      pot_EAM[ivec].Nrho_s = 0;
      pot_EAM[ivec].Nrho_p = 0;
      pot_EAM[ivec].Nrho_d = 0;
      pot_EAM[ivec].dr     = 0.0;
      pot_EAM[ivec].drho_s = 0.0;
      pot_EAM[ivec].drho_p = 0.0;
      pot_EAM[ivec].drho_d = 0.0;
  


      if (iacs.name(s1,s2)=="EAM-d") eam=1;
      if (iacs.name(s1,s2)=="EAM-s") eam=2;
      if (iacs.name(s1,s2)=="EAM-p") eam=3;
      if (iacs.name(s1,s2)=="EAM-sd") eam=4;
      if (iacs.name(s1,s2)=="EAM-pd") eam=5;
      if (iacs.name(s1,s2)=="EAM-spd") eam=6;

  
      /*
	EAM-d:

	Comment line
	Z1element Z2element mass1 mass2 latpar1 latpar2 latname1 latname2
	Nrhod drhod Nr dr rcut 
	(Nrhod points of Fd)
	(Nr points of V2)
	(Nr points of rhod)

	EAM-s: as for EAM-d, but d => s

	EAM-sd:

	Comment line
	Z1element Z2element mass1 mass2 latpar1 latpar2 latname1 latname2
	Nrhod drhod Nr dr rcut Nrhos drhos
	(Nrhod points of Fd)
	(Nr points of V2)
	(Nr points of rhod)
	(Nrhos points of Fs)
	(Nr points of rhos)

	EAM-pd: as for EAM-sd, but s => p

	EAM-spd:

	Comment line
	Z1element Z2element mass1 mass2 latpar1 latpar2 latname1 latname2
	Nrhod drhod Nr dr rcut Nrhos drhos Nrhop drhop
	(Nrhod points of Fd)
	(Nr points of V2)
	(Nr points of rhod)
	(Nrhos points of Fs)
	(Nr points of rhos)
	(Nrhop points of Fp)
	(Nr points of rhop)
    
      */

      filename = "eam." + s1 + "." + s2 + ".in";
      fp.open(filename.c_str());
      if (!fp){
	std::cout << "ERROR: Could not find file " << filename << ". Trying symmetric version ..." << std::endl;
	filename = "eam." + s2 + "." + s1 + ".in";
	aborterror("ERROR: Could not find file " + filename + ". Exiting.");
      }



      utils::get_line(fp, line);
      utils::get_substrings( line, args, "\t :,()[]=");


      utils::get_line(fp, line);
      utils::get_substrings( line, args, "\t :,()[]=");


      if (eam==1 || eam==4 || eam==5 || eam==6){
	utils::get_line(fp, line);
	utils::get_substrings( line, args, "\t :,()[]=");

	strbuf.str(args[0]); strbuf >> pot_EAM[ivec].Nrho_d; strbuf.clear();
	strbuf.str(args[1]); strbuf >> pot_EAM[ivec].drho_d; strbuf.clear();
	strbuf.str(args[2]); strbuf >> pot_EAM[ivec].Nr;   strbuf.clear();
	strbuf.str(args[3]); strbuf >> pot_EAM[ivec].dr;   strbuf.clear();
	strbuf.str(args[4]); strbuf >> pot_EAM[ivec].mrcut; strbuf.clear();
      }
      else if (eam==2){
	utils::get_line(fp, line);
	utils::get_substrings( line, args, "\t :,()[]=");

	strbuf.str(args[0]); strbuf >> pot_EAM[ivec].Nrho_s; strbuf.clear();
	strbuf.str(args[1]); strbuf >> pot_EAM[ivec].drho_s; strbuf.clear();
	strbuf.str(args[2]); strbuf >> pot_EAM[ivec].Nr;   strbuf.clear();
	strbuf.str(args[3]); strbuf >> pot_EAM[ivec].dr;   strbuf.clear();
	strbuf.str(args[4]); strbuf >> pot_EAM[ivec].mrcut; strbuf.clear();
      }
      else if (eam==3){
	utils::get_line(fp, line);
	utils::get_substrings( line, args, "\t :,()[]=");

	strbuf.str(args[0]); strbuf >> pot_EAM[ivec].Nrho_p; strbuf.clear();
	strbuf.str(args[1]); strbuf >> pot_EAM[ivec].drho_p; strbuf.clear();
	strbuf.str(args[2]); strbuf >> pot_EAM[ivec].Nr;   strbuf.clear();
	strbuf.str(args[3]); strbuf >> pot_EAM[ivec].dr;   strbuf.clear();
	strbuf.str(args[4]); strbuf >> pot_EAM[ivec].mrcut; strbuf.clear();
      }


		    
      if (eam==4 || eam==6){
	strbuf.str(args[5]); strbuf >> pot_EAM[ivec].Nrho_s; strbuf.clear();
	strbuf.str(args[6]); strbuf >> pot_EAM[ivec].drho_s; strbuf.clear();
      }
      else if (eam==5){
	strbuf.str(args[5]); strbuf >> pot_EAM[ivec].Nrho_p; strbuf.clear();
	strbuf.str(args[6]); strbuf >> pot_EAM[ivec].drho_p; strbuf.clear();
      }
		    
      if (eam==6){
	strbuf.str(args[7]); strbuf >> pot_EAM[ivec].Nrho_p; strbuf.clear();
	strbuf.str(args[8]); strbuf >> pot_EAM[ivec].drho_p; strbuf.clear();
      }




      /* 
	 pot_EAM[ivec].Nr as read from the EAM data file gives us positions
     
	 r = 0, dr, ..., (Nr-1)*dr
     
	 The cutoff is
     
	 rcut
     
	 We want r to go from 0 to rcut, interval end points included.
     
      */


      /*
	if (pot_EAM[ivec].Nr * pot_EAM[ivec].dr <= pot_EAM[ivec].rcut){
	std::cout << "ERROR: The number of points for (r, V2) is too small for interaction "
	<< i1 << "-" << i2 << ". Exiting." << std::endl;
	exit(EXIT_FAILURE);
	}
      */





      pot_EAM[ivec].r.resize(pot_EAM[ivec].Nr);
      pot_EAM[ivec].V2.resize(pot_EAM[ivec].Nr);
      pot_EAM[ivec].d2_V2.resize(pot_EAM[ivec].Nr);
		

      for (k=0; k<pot_EAM[ivec].Nr; k++)
	pot_EAM[ivec].r[k] = k * pot_EAM[ivec].dr;






      if (eam==1 || eam==4 || eam==5 || eam==6){
	pot_EAM[ivec].F_d.resize(       pot_EAM[ivec].Nrho_d);
	pot_EAM[ivec].F_d_rho_d.resize( pot_EAM[ivec].Nrho_d);
	for (k=0; k<pot_EAM[ivec].Nrho_d; k++)
	  pot_EAM[ivec].F_d_rho_d[k] = k * pot_EAM[ivec].drho_d;

	pot_EAM[ivec].d2_F_d.resize(   pot_EAM[ivec].Nrho_d);
	pot_EAM[ivec].rho_d.resize(    pot_EAM[ivec].Nr);
	pot_EAM[ivec].d2_rho_d.resize( pot_EAM[ivec].Nr);
      }

      if (eam==2 || eam==4 || eam==6){
	pot_EAM[ivec].F_s.resize(       pot_EAM[ivec].Nrho_s);
	pot_EAM[ivec].F_s_rho_s.resize( pot_EAM[ivec].Nrho_s);
	for (k=0; k<pot_EAM[ivec].Nrho_s; k++)
	  pot_EAM[ivec].F_s_rho_s[k] = k * pot_EAM[ivec].drho_s;

	pot_EAM[ivec].d2_F_s.resize(   pot_EAM[ivec].Nrho_s);
	pot_EAM[ivec].rho_s.resize(    pot_EAM[ivec].Nr);
	pot_EAM[ivec].d2_rho_s.resize( pot_EAM[ivec].Nr);
      }

      if (eam==3 || eam==5 || eam==6){
	pot_EAM[ivec].F_p.resize(       pot_EAM[ivec].Nrho_p);
	pot_EAM[ivec].F_p_rho_p.resize( pot_EAM[ivec].Nrho_p);
	for (k=0; k<pot_EAM[ivec].Nrho_p; k++)
	  pot_EAM[ivec].F_p_rho_p[k] = k * pot_EAM[ivec].drho_p;

	pot_EAM[ivec].d2_F_p.resize(   pot_EAM[ivec].Nrho_p);
	pot_EAM[ivec].rho_p.resize(    pot_EAM[ivec].Nr);
	pot_EAM[ivec].d2_rho_p.resize( pot_EAM[ivec].Nr);
      }






      /* --------------------------------------------------------------
	 Read secondary data.
	 -------------------------------------------------------------- */

      if (eam==1 || eam==4 || eam==5 || eam==6){

	for (k=0; k<pot_EAM[ivec].Nrho_d; ++k){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].F_d[k];
	  strbuf.clear();
	}

      }
      else if (eam==2){
	for (k=0; k<pot_EAM[ivec].Nrho_s; ++k){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].F_s[k];
	  strbuf.clear();
	}

      }
      else if (eam==3){

	for (k=0; k<pot_EAM[ivec].Nrho_p; ++k){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].F_p[k];
	  strbuf.clear();
	}

      }














      for (k=0; k<pot_EAM[ivec].Nr; ++k){
	utils::get_line(fp, line);
	utils::get_substrings( line, args, "\t :,()[]=");
	strbuf.str(args[0]); strbuf >> pot_EAM[ivec].V2[k];
	strbuf.clear();
      }





      if (eam==1 || eam==4 || eam==5 || eam==6){
	for (k=0; k<pot_EAM[ivec].Nr; ++k){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].rho_d[k];
	  strbuf.clear();
	}
      }
      else if (eam==2){
	for (k=0; k<pot_EAM[ivec].Nr; k++){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].rho_s[k];
	  strbuf.clear();
	}
      }
      else if (eam==3){
	for (k=0; k<pot_EAM[ivec].Nr; k++){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].rho_p[k];
	  strbuf.clear();
	}
      }





      if (eam==4 || eam==6){
	for (k=0; k<pot_EAM[ivec].Nrho_s; k++){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].F_s[k];
	  strbuf.clear();
	}
      }
      else if (eam==5){
	for (k=0; k<pot_EAM[ivec].Nrho_p; k++){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].F_p[k];
	  strbuf.clear();
	}
      }

      if (eam==6){
	for (k=0; k<pot_EAM[ivec].Nrho_p; k++){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].F_p[k];
	  strbuf.clear();
	}
      }





      if (eam==4 || eam==6){
	for (k=0; k<pot_EAM[ivec].Nr; k++){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].rho_s[k];
	  strbuf.clear();
	}
      }
      else if (eam==5){
	for (k=0; k<pot_EAM[ivec].Nr; k++){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].rho_p[k];
	  strbuf.clear();
	}
      }

      if (eam==6){
	for (k=0; k<pot_EAM[ivec].Nr; k++){
	  utils::get_line(fp, line);
	  utils::get_substrings( line, args, "\t :,()[]=");
	  strbuf.str(args[0]); strbuf >> pot_EAM[ivec].rho_p[k];
	  strbuf.clear();
	}
      }
      fp.close();
      fp.clear();




      std::cout << "Interaction " << s1 << "-" << s2 << ": Read in " << pot_EAM[ivec].Nr << " "
	   << "data points for V2 in EAM, with rcut " << pot_EAM[ivec].rcut() << std::endl;
      if (pot_EAM[ivec].Nrho_s > 0)
	std::cout << "Interaction " << s1 << "-" << s2 << ": Read in " << pot_EAM[ivec].Nrho_s << " "
	     << "data points for Fs in EAM, with rcut " << pot_EAM[ivec].rcut() << std::endl;
      if (pot_EAM[ivec].Nrho_p > 0)
	std::cout << "Interaction " << s1 << "-" << s2 << ": Read in " << pot_EAM[ivec].Nrho_p << " "
	     << "data points for Fp in EAM, with rcut " << pot_EAM[ivec].rcut() << std::endl;
      if (pot_EAM[ivec].Nrho_d > 0)
	std::cout << "Interaction " << s1 << "-" << s2 << ": Read in " << pot_EAM[ivec].Nrho_d << " "
	     << "data points for Fd in EAM, with rcut " << pot_EAM[ivec].rcut() << std::endl;






      /* --------------------------------------------------------------
	 Get the spline tables.
	 -------------------------------------------------------------- */
      spline(pot_EAM[ivec].r,
	     pot_EAM[ivec].V2, 
	     pot_EAM[ivec].d2_V2);
		  
      if (eam==2 || eam==4 || eam==6){
	spline(pot_EAM[ivec].r,
	       pot_EAM[ivec].rho_s, 
	       pot_EAM[ivec].d2_rho_s);
	spline(pot_EAM[ivec].F_s_rho_s,
	       pot_EAM[ivec].F_s,
	       pot_EAM[ivec].d2_F_s);
      }

      if (eam==3 || eam==5 || eam==6){
	spline(pot_EAM[ivec].r,
	       pot_EAM[ivec].rho_p, 
	       pot_EAM[ivec].d2_rho_p);
	spline(pot_EAM[ivec].F_p_rho_p,
	       pot_EAM[ivec].F_p,
	       pot_EAM[ivec].d2_F_p);
      }


      if (eam==1 || eam==4 || eam==5 || eam==6){
	spline(pot_EAM[ivec].r,
	       pot_EAM[ivec].rho_d, 
	       pot_EAM[ivec].d2_rho_d);
	spline(pot_EAM[ivec].F_d_rho_d,
	       pot_EAM[ivec].F_d,
	       pot_EAM[ivec].d2_F_d);
      }





	  

      /* Print out data for manual checking: */

      filename = "eampot." + s1 + "." + s2 + ".V2.test";
      fpo.open(filename.c_str());

      dr = pot_EAM[ivec].r[ pot_EAM[ivec].Nr ]/100.0;
      r = 0.0;
      while (r <= pot_EAM[ivec].r[ pot_EAM[ivec].Nr-1 ]){
	fpo << r
	    << " "
	    << splint(pot_EAM[ivec].r,
		      pot_EAM[ivec].V2, 
		      pot_EAM[ivec].d2_V2,
		      r)
	    << " "
	    << splint_dy(pot_EAM[ivec].r,
			 pot_EAM[ivec].V2, 
			 pot_EAM[ivec].d2_V2,
			 r)
	    << std::endl;
	r+=dr;
      }
      fpo.close();
      fpo.clear();

  


      // --------------------------------------------------------------------
      // --------------------------------------------------------------------
      if (eam==2 || eam==4 || eam==6){
	filename = "eampot." + s1 + "." + s2 + ".rho_s.test";
	fpo.open(filename.c_str());

	dr = pot_EAM[ivec].r[ pot_EAM[ivec].Nr ]/100.0;
	r = 0.0;
	while (r <= pot_EAM[ivec].r[ pot_EAM[ivec].Nr-1 ]){
	  fpo << r
	      << " "
	      << splint(pot_EAM[ivec].r,
			pot_EAM[ivec].rho_s, 
			pot_EAM[ivec].d2_rho_s,
			r)
	      << " "
	      << splint_dy(pot_EAM[ivec].r,
			   pot_EAM[ivec].rho_s, 
			   pot_EAM[ivec].d2_rho_s,
			   r)
	      << std::endl;
	  r+=dr;
	}
	fpo.close();
	fpo.clear();

	// --------------------------------------------------------------------
	filename = "eampot." + s1 + "." + s2 + ".F_s.test";
	fpo.open(filename.c_str());

	drho = pot_EAM[ivec].rho_s[ pot_EAM[ivec].Nrho_s ]/100.0;
	rho = 0.0;
	while (rho <= pot_EAM[ivec].rho_s[ pot_EAM[ivec].Nrho_s-1 ]){
	  fpo << rho
	      << " "
	      << splint(pot_EAM[ivec].F_s_rho_s,
			pot_EAM[ivec].F_s,
			pot_EAM[ivec].d2_F_s,
			rho)
	      << " "
	      << splint_dy(pot_EAM[ivec].F_s_rho_s,
			   pot_EAM[ivec].F_s, 
			   pot_EAM[ivec].d2_F_s,
			   r)
	      << std::endl;
	  rho+=drho;
	}
	fpo.close();
	fpo.clear();

      }
      // --------------------------------------------------------------------
      // --------------------------------------------------------------------
      if (eam==3 || eam==5 || eam==6){
	filename = "eampot." + s1 + "." + s2 + ".rho_p.test";
	fpo.open(filename.c_str());

	dr = pot_EAM[ivec].r[ pot_EAM[ivec].Nr ]/100.0;
	r = 0.0;
	while (r <= pot_EAM[ivec].r[ pot_EAM[ivec].Nr-1 ]){
	  fpo << r
	      << " "
	      << splint(pot_EAM[ivec].r,
			pot_EAM[ivec].rho_p, 
			pot_EAM[ivec].d2_rho_p,
			r)
	      << " "
	      << splint_dy(pot_EAM[ivec].r,
			   pot_EAM[ivec].rho_p, 
			   pot_EAM[ivec].d2_rho_p,
			   r)
	      << std::endl;
	  r+=dr;
	}
	fpo.close();
	fpo.clear();

	// --------------------------------------------------------------------
	filename = "eampot." + s1 + "." + s2 + ".F_p.test";
	fpo.open(filename.c_str());

	drho = pot_EAM[ivec].rho_p[ pot_EAM[ivec].Nrho_p ]/100.0;
	rho = 0.0;
	while (rho <= pot_EAM[ivec].rho_p[ pot_EAM[ivec].Nrho_p-1 ]){
	  fpo << rho
	      << " "
	      << splint(pot_EAM[ivec].F_p_rho_p,
			pot_EAM[ivec].F_p,
			pot_EAM[ivec].d2_F_p,
			rho)
	      << " "
	      << splint_dy(pot_EAM[ivec].F_p_rho_p,
			   pot_EAM[ivec].F_p, 
			   pot_EAM[ivec].d2_F_p,
			   r)
	      << std::endl;
	  rho+=drho;
	}
	fpo.close();
	fpo.clear();
      }
      // --------------------------------------------------------------------
      // --------------------------------------------------------------------
      if (eam==1 || eam==4 || eam==5 || eam==6){
	filename = "eampot." + s1 + "." + s2 + ".rho_d.test";
	fpo.open(filename.c_str());

	dr = pot_EAM[ivec].r[ pot_EAM[ivec].Nr ]/100.0;
	r = 0.0;
	while (r <= pot_EAM[ivec].r[ pot_EAM[ivec].Nr-1 ]){
	  fpo << r
	      << " "
	      << splint(pot_EAM[ivec].r,
			pot_EAM[ivec].rho_d, 
			pot_EAM[ivec].d2_rho_d,
			r)
	      << " "
	      << splint_dy(pot_EAM[ivec].r,
			   pot_EAM[ivec].rho_d, 
			   pot_EAM[ivec].d2_rho_d,
			   r)
	      << std::endl;
	  r+=dr;
	}
	fpo.close();
	fpo.clear();
	// --------------------------------------------------------------------
	filename = "eampot." + s1 + "." + s2 + ".F_d.test";
	fpo.open(filename.c_str());

	drho = pot_EAM[ivec].rho_d[ pot_EAM[ivec].Nrho_d ]/100.0;
	rho = 0.0;
	while (rho <= pot_EAM[ivec].rho_d[ pot_EAM[ivec].Nrho_d-1 ]){
	  fpo << rho
	      << " "
	      << splint(pot_EAM[ivec].F_d_rho_d,
			pot_EAM[ivec].F_d,
			pot_EAM[ivec].d2_F_d,
			rho)
	      << " "
	      << splint_dy(pot_EAM[ivec].F_d_rho_d,
			   pot_EAM[ivec].F_d, 
			   pot_EAM[ivec].d2_F_d,
			   r)
	      << std::endl;
	  rho+=drho;
	}
	fpo.close();
	fpo.clear();
      }


      nread++;
    }
  }


  std::cout << "Read in " << nread << " EAM potentials." << std::endl;
}
  
