


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
#include "utils-streamio.hpp"
#include "utils-errors.hpp"

#include "helpfuns.hpp"
#include "potinfo.hpp"
#include "potclasses.hpp"
#include "utils-string.hpp"

using namespace utils;
using boost::format;


void PotentialInformation::read_reppot(void){


  std::cout << "Reading reppots ..." << std::endl;

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
      std::cout << "Checking " << s1 << "-" << s2 << " ..." << std::endl;
      

      if (! use_reppot(s1,s2)){
	cout << "No reppot can be read for combination " << s1 << " " << s2 << std::endl;
	continue;
      }
      

      int ivec = reppot_vecidx(s1, s2);


      std::string filename = "reppot." + s1 + "." + s2 + ".in";

      double rmin, r, dr;

      std::ifstream fp;
      std::ofstream fpo;
      std::string line;
      std::vector<std::string> args;
      std::istringstream strbuf;
      int ns;


      fp.open(filename.c_str());
      if (!fp)
	aborterror("Error: Could not find file " + filename + ". Exiting.");

      int nlines = 0;
      while (true){
	utils::get_line(fp, line);
	ns = utils::get_substrings( line, args, "\t :,()[]=");
	if (ns==0 && fp) continue;

	if (ns>=2){
	  if (nlines==0){
	    strbuf.str(args[0]);
	    strbuf >> rmin;
	    strbuf.clear();
	  }
	  nlines++;
	}

	if (! fp) break;


      }
      fp.close();
      fp.clear();




      pot_Reppot[ivec].elemname1 = elem.idx2name(i1);
      pot_Reppot[ivec].elemname2 = elem.idx2name(i2);


      /* ***************************************************************************
	 Make sure the minimum r value rmin is larger than zero. Else V(rmin) will
	 be infinite, which will cause problems.
	 *************************************************************************** */

      int Nr_rep = nlines;
      bool skip_first_value = false;
      if (fp_is_small(rmin)){
	skip_first_value = true;
	Nr_rep--; // first r value is 0, discard this
      }
      // Allocate memory:
      pot_Reppot[ivec].Nr_rep = Nr_rep;
      pot_Reppot[ivec].r_rep.resize(Nr_rep);
      pot_Reppot[ivec].V_rep.resize(Nr_rep);
      pot_Reppot[ivec].d2_V_rep.resize(Nr_rep);


      k=0;
      if (skip_first_value) k=-1;

      fp.open(filename.c_str());

      while (true){
	utils::get_line(fp, line);
	ns = utils::get_substrings( line, args, "\t :,()[]=");
	if (ns==0 && fp) continue;

	if (ns>=2){
	  if (k==-1 && skip_first_value) { k=0; continue; }

	  strbuf.str(args[0]);
	  strbuf >>   pot_Reppot[ivec].r_rep[k];
	  strbuf.clear();

	  strbuf.str(args[1]);
	  strbuf >>   pot_Reppot[ivec].V_rep[k];
	  strbuf.clear();

	  k++;
	}


	if (!fp) break; 
      }
      fp.close();
      fp.clear();


#if 0
      if (! fp_is_small(rmin)){ // first r value is not 0
	pot_Reppot[ivec].r_rep[0] = 0.0;

	double a, b;
	a = (pot_Reppot[ivec].V_rep[2] - pot_Reppot[ivec].V_rep[1])/
	  (pot_Reppot[ivec].r_rep[2] - pot_Reppot[ivec].r_rep[1]);
	b = 0.5 * ( pot_Reppot[ivec].V_rep[1] + pot_Reppot[ivec].V_rep[2]
		    - a * (pot_Reppot[ivec].r_rep[1] + pot_Reppot[ivec].r_rep[2]) );
	pot_Reppot[ivec].V_rep[0] = a * pot_Reppot[ivec].r_rep[0] + b;
      }
#endif

              
      /* Get the spline table: */
      spline(pot_Reppot[ivec].r_rep,
	     pot_Reppot[ivec].V_rep, 
	     pot_Reppot[ivec].d2_V_rep);


      std::cout << "Largest r value: " << pot_Reppot[ivec].r_rep[ pot_Reppot[ivec].Nr_rep-1 ] << std::endl;


      /* Print out data for manual checking: */

      filename = "reppot." + s1 + "." + s2 + ".test";
      fpo.open(filename.c_str());

      dr = pot_Reppot[ivec].r_rep[ pot_Reppot[ivec].Nr_rep-1 ] /
	(3 * pot_Reppot[ivec].Nr_rep);

      r = rmin;
      while (r <= pot_Reppot[ivec].r_rep[ pot_Reppot[ivec].Nr_rep-1 ]){
	fpo << r << " "
	    << splint(pot_Reppot[ivec].r_rep,
		      pot_Reppot[ivec].V_rep, 
		      pot_Reppot[ivec].d2_V_rep,
		      r) << " "
	    << splint_dy(pot_Reppot[ivec].r_rep,
			 pot_Reppot[ivec].V_rep, 
			 pot_Reppot[ivec].d2_V_rep,
			 r)
	    << std::endl;
	r += dr;

      }
      fpo.close();
      fpo.clear();



      std::cout << "Repulsive potential core for " << s1 << "-" << s2 << " interaction read in. "
	   << "rmin = " << pot_Reppot[ivec].r_rep[ 0 ] << " "
	   << "rmax = " << pot_Reppot[ivec].r_rep[ pot_Reppot[ivec].Nr_rep-1 ]
	   << std::endl;

    }
  }



  return;
}




