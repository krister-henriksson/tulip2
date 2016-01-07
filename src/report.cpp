
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include <cmath>

#include "chisq-basics.hpp"
#include "funcfit-basics.hpp"
#include "funcfit-conjgrad.hpp"
#include "funcfit-errors.hpp"
#include "funcfit-gauss-newton.hpp"
#include "funcfit-leve-marq.hpp"
#include "funcfit-powelldogleg.hpp"
#include "funcfit-simplexfit.hpp"
#include "constants.hpp"
#include "nr-f1dim.hpp"
#include "nr-golden.hpp"
#include "nr-linemethod.hpp"
#include "param.hpp"
#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-matrix3.hpp"
#include "utils-matrix-Choleskydecomp.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix-LUdecomp.hpp"
#include "utils-matrix-QRdecomp.hpp"
#include "utils-string.hpp"
#include "utils-vector.hpp"
#include "utils-errors.hpp"

#include "atomsystem.hpp"
#include "compound.hpp"
#include "elem-iacs.hpp"
#include "helpfuns.hpp"
#include "mdsystem.hpp"
#include "mtwister.hpp"
#include "param-pot.hpp"
#include "physconst.hpp"
#include "potclasses.hpp"
#include "potinfo.hpp"
#include "propfun.hpp"
#include "specs-fit-prop-pot.hpp"
#include "get-comp-prop.hpp"
#include "report.hpp"


using namespace utils;
using namespace constants;
using namespace physconst;
using namespace funcfit;
using boost::format;



// This signature is required. Otherwise the function cannot be
// set as a reporting function for ChiSq objects.
void report_pot_prop(ParamPot & param,
		     Vector<CompoundStructureFit> & DX,
		     Vector<double> & DY,
		     Vector<double> & MDY
		     ){
  report_pot_prop_ext(param, DX, DY, MDY, std::cout, std::cout);
}


void report_pot_prop_ext(ParamPot & param,
			 Vector<CompoundStructureFit> & DX,
			 Vector<double> & DY,
			 Vector<double> & MDY,
			 std::ostream & fout_pot,
			 std::ostream & fout_prop
			 ){

  std::cout << "-------------------------------------------------------------------------------" << std::endl;


  param.update_pot();


  // Report the current settings of the potentials:
  report_pot( param.p_potinfo, true, false, fout_pot );

  // Report the current properties of the compounds (assumes they have been
  // calculated using the current settings of the potentials):
  report_prop( DX, fout_prop );


  //  if (param.p_potinfo->specs_prop.mds_specs_common.quick_mode)
  //  aborterror("Done with quick mode, performing dirty exit.");
}



/*

  potentials: - non-fittable                        X
              - fittable
	         * fixed parameters                 X
		 * free with limits parameters            X
		 * free parameters                        X

 */


void report_pot(PotentialInformationFit * p_potinfo,
		bool fittable_ones,
		bool fixed_ones,
		std::ostream & fout
		){
  int i1, i2, i3, j, k, nel = p_potinfo->elem.nelem();
  std::string s1, s2, s3, partypestring;
  bool o1, o2, oprint;
  std::string formatf = "%20.10f";
  std::string formate = "%20.10e";
  const double llim=1.0e-4;
  const double ulim=1.0e+4;
  double td;

  if (! fittable_ones && ! fixed_ones) return;



  // ##########################################################################
  // Report on potential
  // ##########################################################################



  for (i1=0; i1<nel; ++i1){
    s1 = p_potinfo->elem.idx2name(i1);
    
    for (i2=i1; i2<nel; ++i2){
      s2 = p_potinfo->elem.idx2name(i2);
	
      j = p_potinfo->basepot_vecidx(s1,s2);
      if (j<0) continue;


      if (p_potinfo->basepot(s1,s2) == "ABOP"){
	/*
	  fout << "ABOP: " << p_potinfo->pot_ABOP[j].elemname1 << "-"
	  << p_potinfo->pot_ABOP[j].elemname2 << ": is fittable?: "
	  << p_potinfo->is_fittable(s1,s2) << std::endl;
	*/

	for (k=0; k<p_potinfo->pot_ABOP[j].parname.size(); ++k){

	  o1 = fixed_ones
	    && 
	    ( 
	     ((p_potinfo->is_fittable(s1,s2)) && (p_potinfo->pot_ABOP[j].partype[k] == PARAM_FIXED))
	     ||
	     (!(p_potinfo->is_fittable(s1,s2)))
	      );
	  o2 = fittable_ones
	    && (p_potinfo->is_fittable(s1,s2)) 
	    && (p_potinfo->pot_ABOP[j].partype[k] != PARAM_FIXED);
	  if (!o1 && !o2) continue;

	  if (o1 || o2){
	    fout << "ABOP "
		 << format("%2s-%2s: ") % s1 % s2
		 << format("%20s : ") % p_potinfo->pot_ABOP[j].parname[k];

	    td = p_potinfo->pot_ABOP[j].parval[k];
	    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
	    else                              fout << format(formatf) % td;
	  }

	  if (o2){
	    partypestring="unknown";
	    if      (p_potinfo->pot_ABOP[j].partype[k] == PARAM_FIXED) partypestring="FIXED parameter";
	    else if (p_potinfo->pot_ABOP[j].partype[k] == PARAM_FREE_WITH_LIMITS) partypestring="CONSTRAINED parameter";
	    else if (p_potinfo->pot_ABOP[j].partype[k] == PARAM_FREE) partypestring="FREE parameter";
	    fout << "   min: " << format("%20.10e") % p_potinfo->pot_ABOP[j].parmin[k]
		 << "   max: " << format("%20.10e") % p_potinfo->pot_ABOP[j].parmax[k]
		 << "   " <<  partypestring;
	  }
	  if (o1 || o2) fout << std::endl;
	} // loop: k

	for (k=0; k<p_potinfo->pot_ABOP[j].rcs.parname.size(); ++k){

	  o1 = fixed_ones
	    && 
	    ( 
	     ((p_potinfo->is_fittable(s1,s2)) && (p_potinfo->pot_ABOP[j].rcs.partype[k] == PARAM_FIXED))
	     ||
	     (!(p_potinfo->is_fittable(s1,s2)))
	      );
	  o2 = fittable_ones
	    && (p_potinfo->is_fittable(s1,s2)) 
	    && (p_potinfo->pot_ABOP[j].rcs.partype[k] != PARAM_FIXED);
	  if (!o1 && !o2) continue;

	  if (o1 || o2){
	    fout << "ABOP "
		 << format("%2s-%2s: ") % s1 % s2
		 << format("%20s : ") % p_potinfo->pot_ABOP[j].rcs.parname[k];

	    td = p_potinfo->pot_ABOP[j].rcs.parval[k];
	    if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
	    else                              fout << format(formatf) % td;
	  }

	  if (o2){
	    partypestring="unknown";
	    if      (p_potinfo->pot_ABOP[j].rcs.partype[k] == PARAM_FIXED) partypestring="FIXED parameter";
	    else if (p_potinfo->pot_ABOP[j].rcs.partype[k] == PARAM_FREE_WITH_LIMITS) partypestring="CONSTRAINED parameter";
	    else if (p_potinfo->pot_ABOP[j].rcs.partype[k] == PARAM_FREE) partypestring="FREE parameter";
	    fout << "   min: " << format("%20.10e") % p_potinfo->pot_ABOP[j].rcs.parmin[k]
		 << "   max: " << format("%20.10e") % p_potinfo->pot_ABOP[j].rcs.parmax[k]
		 << "   " <<  partypestring;
	  }
	  if (o1 || o2) fout << std::endl;
	} // loop: k

      } // ABOP
    } // loop: i2
  }


  //if (counter==0) fout << "(no items detected)" << std::endl;

  for (i1=0; i1<nel; ++i1){
    s1 = p_potinfo->elem.idx2name(i1);
    for (i2=0; i2<nel; ++i2){
      s2 = p_potinfo->elem.idx2name(i2);
      for (i3=0; i3<nel; ++i3){
	s3 = p_potinfo->elem.idx2name(i3);


	if (! (*p_potinfo).use_abop_alpha.elem(i1,i2,i3)) continue;

	o1 = fixed_ones    && (*p_potinfo).abop_alpha_partype.elem(i1,i2,i3) == PARAM_FIXED;
	o2 = fittable_ones && (*p_potinfo).abop_alpha_partype.elem(i1,i2,i3) != PARAM_FIXED;
	if (!o1 && !o2) continue;

	if (o1 || o2){
	  fout << "ABOP "
	       << format("%2s-%2s-%2s: ") % s1 % s2 % s3
	       << format("%17s : ") % "alpha";

	  td = p_potinfo->abop_alpha.elem(i1,i2,i3);
	  if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
	  else                              fout << format(formatf) % td;
	}

	if (o2){
	  partypestring="unknown";
	  if      ((*p_potinfo).abop_alpha_partype.elem(i1,i2,i3) == PARAM_FIXED) partypestring="FIXED parameter";
	  else if ((*p_potinfo).abop_alpha_partype.elem(i1,i2,i3) == PARAM_FREE_WITH_LIMITS) partypestring="CONSTRAINED parameter";
	  else if ((*p_potinfo).abop_alpha_partype.elem(i1,i2,i3) == PARAM_FREE) partypestring="FREE parameter";
	  fout << "   min: " << format("%20.10e") % p_potinfo->abop_alpha_parmin.elem(i1,i2,i3)
	       << "   max: " << format("%20.10e") % p_potinfo->abop_alpha_parmax.elem(i1,i2,i3)
	       << "   " <<  partypestring;
	}
	if (o1 || o2) fout << std::endl;

      }
    }
  }


  for (i1=0; i1<nel; ++i1){
    s1 = p_potinfo->elem.idx2name(i1);
    for (i2=0; i2<nel; ++i2){
      s2 = p_potinfo->elem.idx2name(i2);
      for (i3=0; i3<nel; ++i3){
	s3 = p_potinfo->elem.idx2name(i3);

	if (! (*p_potinfo).use_abop_omega.elem(i1,i2,i3)) continue;

	o1 = fixed_ones    && (*p_potinfo).abop_omega_partype.elem(i1,i2,i3) == PARAM_FIXED;
	o2 = fittable_ones && (*p_potinfo).abop_omega_partype.elem(i1,i2,i3) != PARAM_FIXED;
	if (!o1 && !o2) continue;

	if (o1 || o2){
	  fout << "ABOP "
	       << format("%2s-%2s-%2s: ") % s1 % s2 % s3
	       << format("%17s : ") % "omega";

	  td = p_potinfo->get_abop_omega(s1,s2,s3);
	  if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
	  else                              fout << format(formatf) % td;
	}

	if (o2){
	  partypestring="unknown";
	  if      ((*p_potinfo).abop_omega_partype.elem(i1,i2,i3) == PARAM_FIXED) partypestring="FIXED parameter";
	  else if ((*p_potinfo).abop_omega_partype.elem(i1,i2,i3) == PARAM_FREE_WITH_LIMITS) partypestring="CONSTRAINED parameter";
	  else if ((*p_potinfo).abop_omega_partype.elem(i1,i2,i3) == PARAM_FREE) partypestring="FREE parameter";
	  fout << "   min: " << format("%20.10e") % p_potinfo->abop_omega_parmin.elem(i1,i2,i3)
	       << "   max: " << format("%20.10e") % p_potinfo->abop_omega_parmax.elem(i1,i2,i3)
	       << "   " <<  partypestring;
	}
	if (o1 || o2) fout << std::endl;

      }
    }
  }


  for (i1=0; i1<nel; ++i1){
    s1 = p_potinfo->elem.idx2name(i1);
    for (i2=0; i2<nel; ++i2){
      s2 = p_potinfo->elem.idx2name(i2);

	if (! (*p_potinfo).use_abop_2mu.elem(i1,i2)) continue;

	o1 = fixed_ones    && (*p_potinfo).abop_2mu_partype.elem(i1,i2) == PARAM_FIXED;
	o2 = fittable_ones && (*p_potinfo).abop_2mu_partype.elem(i1,i2) != PARAM_FIXED;
	if (!o1 && !o2) continue;

	if (o1 || o2){
	  fout << "ABOP "
	       << format("%2s-%2s: ") % s1 % s2
	       << format("%20s : ") % "2mu";

	  td = p_potinfo->abop_2mu.elem(i1,i2);
	  if (abs(td)<llim || abs(td)>ulim) fout << format(formate) % td;
	  else                              fout << format(formatf) % td;
	}

	if (o2){
	  partypestring="unknown";
	  if      ((*p_potinfo).abop_2mu_partype.elem(i1,i2) == PARAM_FIXED) partypestring="FIXED parameter";
	  else if ((*p_potinfo).abop_2mu_partype.elem(i1,i2) == PARAM_FREE_WITH_LIMITS) partypestring="CONSTRAINED parameter";
	  else if ((*p_potinfo).abop_2mu_partype.elem(i1,i2) == PARAM_FREE) partypestring="FREE parameter";
	  fout << "   min: " << format("%20.10e") % p_potinfo->abop_2mu_parmin.elem(i1,i2)
	       << "   max: " << format("%20.10e") % p_potinfo->abop_2mu_parmax.elem(i1,i2)
	       << "   " <<  partypestring;
	}
	if (o1 || o2) fout << std::endl;

    }
  }




  fout << "-------------------------------------------------------------------------------" << std::endl;

}





void report_prop(Vector<CompoundStructureFit> & DX,
		 std::ostream & fout,
		 bool firsttime
		 ){

  int i;
  std::string s1, s2, s3;
  double td1, td2, td3, td4;
  int k,p;
  bool tb1, tb2;
  std::string propstr="empty";
  std::string compstr;


  // ##########################################################################
  // Report on compounds
  // ##########################################################################


  //  std::string compstr = stream.str();


  for (i=0; i<DX.size(); ++i){
    CompoundStructureFit cmpfit = DX[i];

    //fout << "Compound: " << cmpfit.name << std::endl;

    if (cmpfit.prop_use.a){
      td1 = cmpfit.prop_pred.a;
      td2 = cmpfit.prop_readin.a;
      tb1 = cmpfit.use_u.a;
      tb2 = cmpfit.use_w.a;
      td3 = cmpfit.prop_u.a;
      td4 = cmpfit.prop_w.a;

      //printf("Lattice parameter a                   : read-in  %15.10f", td2);
      //printf("  predicted  %15.10f  rel. change  %10.3e", td1, fp_divide(td1-td2, td2));
      //printf ("  uncertainty  %15.10f", td3);
      //printf ("  weight       %15.10f", td4);
      //printf("\n");

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Lattice parameter a                   : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.b){
      td1 = cmpfit.prop_pred.b;
      td2 = cmpfit.prop_readin.b;
      tb1 = cmpfit.use_u.b;
      tb2 = cmpfit.use_w.b;
      td3 = cmpfit.prop_u.b;
      td4 = cmpfit.prop_w.b;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Lattice parameter b                   : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.c){
      td1 = cmpfit.prop_pred.c;
      td2 = cmpfit.prop_readin.c;
      tb1 = cmpfit.use_u.c;
      tb2 = cmpfit.use_w.c;
      td3 = cmpfit.prop_u.c;
      td4 = cmpfit.prop_w.c;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Lattice parameter c                   : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.bpa){
      td1 = cmpfit.prop_pred.bpa;
      td2 = cmpfit.prop_readin.bpa;
      tb1 = cmpfit.use_u.bpa;
      tb2 = cmpfit.use_w.bpa;
      td3 = cmpfit.prop_u.bpa;
      td4 = cmpfit.prop_w.bpa;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Lattice parameter ratio b/a           : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.cpa){
      td1 = cmpfit.prop_pred.cpa;
      td2 = cmpfit.prop_readin.cpa;
      tb1 = cmpfit.use_u.cpa;
      tb2 = cmpfit.use_w.cpa;
      td3 = cmpfit.prop_u.cpa;
      td4 = cmpfit.prop_w.cpa;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Lattice parameter ratio c/a           : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.r0){
      td1 = cmpfit.prop_pred.r0;
      td2 = cmpfit.prop_readin.r0;
      tb1 = cmpfit.use_u.r0;
      tb2 = cmpfit.use_w.r0;
      td3 = cmpfit.prop_u.r0;
      td4 = cmpfit.prop_w.r0;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Dimer bond length r0                  : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.angle_ab){
      td1 = cmpfit.prop_pred.angle_ab;
      td2 = cmpfit.prop_readin.angle_ab;
      tb1 = cmpfit.use_u.angle_ab;
      tb2 = cmpfit.use_w.angle_ab;
      td3 = cmpfit.prop_u.angle_ab;
      td4 = cmpfit.prop_w.angle_ab;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Angle btw lat param a and b           : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.angle_ac){
      td1 = cmpfit.prop_pred.angle_ac;
      td2 = cmpfit.prop_readin.angle_ac;
      tb1 = cmpfit.use_u.angle_ac;
      tb2 = cmpfit.use_w.angle_ac;
      td3 = cmpfit.prop_u.angle_ac;
      td4 = cmpfit.prop_w.angle_ac;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Angle btw lat param a and c           : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.angle_bc){
      td1 = cmpfit.prop_pred.angle_bc;
      td2 = cmpfit.prop_readin.angle_bc;
      tb1 = cmpfit.use_u.angle_bc;
      tb2 = cmpfit.use_w.angle_bc;
      td3 = cmpfit.prop_u.angle_bc;
      td4 = cmpfit.prop_w.angle_bc;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Angle btw lat param b and c           : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.Vatom){
      td1 = cmpfit.prop_pred.Vatom;
      td2 = cmpfit.prop_readin.Vatom;
      tb1 = cmpfit.use_u.Vatom;
      tb2 = cmpfit.use_w.Vatom;
      td3 = cmpfit.prop_u.Vatom;
      td4 = cmpfit.prop_w.Vatom;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Atomic volume Vatom                   : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.Ecoh){
      td1 = cmpfit.prop_pred.Ecoh;
      td2 = cmpfit.prop_readin.Ecoh;
      tb1 = cmpfit.use_u.Ecoh;
      tb2 = cmpfit.use_w.Ecoh;
      td3 = cmpfit.prop_u.Ecoh;
      td4 = cmpfit.prop_w.Ecoh;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Cohesive energy Ecoh                  : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.Ecoh_delta){
      td1 = cmpfit.prop_pred.Ecoh_delta;
      td2 = cmpfit.prop_readin.Ecoh_delta;
      tb1 = cmpfit.use_u.Ecoh_delta;
      tb2 = cmpfit.use_w.Ecoh_delta;
      td3 = cmpfit.prop_u.Ecoh_delta;
      td4 = cmpfit.prop_w.Ecoh_delta;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Change in cohesive energy Ecoh_delta  : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.Emix){
      td1 = cmpfit.prop_pred.Emix;
      td2 = cmpfit.prop_readin.Emix;
      tb1 = cmpfit.use_u.Emix;
      tb2 = cmpfit.use_w.Emix;
      td3 = cmpfit.prop_u.Emix;
      td4 = cmpfit.prop_w.Emix;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Mixing energy Emix                    : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.B){
      td1 = cmpfit.prop_pred.B;
      td2 = cmpfit.prop_readin.B;
      tb1 = cmpfit.use_u.B;
      tb2 = cmpfit.use_w.B;
      td3 = cmpfit.prop_u.B;
      td4 = cmpfit.prop_w.B;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Bulk modulus B                        : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    if (cmpfit.prop_use.Bp){
      td1 = cmpfit.prop_pred.Bp;
      td2 = cmpfit.prop_readin.Bp;
      tb1 = cmpfit.use_u.Bp;
      tb2 = cmpfit.use_w.Bp;
      td3 = cmpfit.prop_u.Bp;
      td4 = cmpfit.prop_w.Bp;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Pressure derivative of bulk modulus B': ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }

    for (k=0; k<6; ++k)
      for (p=0; p<6; ++p)
	if (cmpfit.prop_use.C.elem(k,p)){
	  td1 = cmpfit.prop_pred.C.elem(k,p);
	  td2 = cmpfit.prop_readin.C.elem(k,p);
	  tb1 = cmpfit.use_u.C.elem(k,p);
	  tb2 = cmpfit.use_w.C.elem(k,p);
	  td3 = cmpfit.prop_u.C.elem(k,p);
	  td4 = cmpfit.prop_w.C.elem(k,p);

	  std::ostringstream sstream;
	  sstream << format("%15s") % cmpfit.name;
	  compstr = "Compound: " + sstream.str() + " : ";

	  propstr = compstr + "Elastic constant C" + tostring(k+1) + tostring(p+1) + "                  : ";
	  print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
	}

    if (cmpfit.prop_use.Fmax){
      td1 = cmpfit.prop_pred.Fmax;
      td2 = cmpfit.prop_readin.Fmax;
      tb1 = cmpfit.use_u.Fmax;
      tb2 = cmpfit.use_w.Fmax;
      td3 = cmpfit.prop_u.Fmax;
      td4 = cmpfit.prop_w.Fmax;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Maximum force Fmax                    : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }
    else if (! cmpfit.prop_use.Fmax){
      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "INFO: Maximum force Fmax              : ";
      print_prop_pred(fout, propstr, cmpfit.prop_pred.Fmax);
    }


    if (cmpfit.prop_use.Pmax){
      td1 = cmpfit.prop_pred.Pmax;
      td2 = cmpfit.prop_readin.Pmax;
      tb1 = cmpfit.use_u.Pmax;
      tb2 = cmpfit.use_w.Pmax;
      td3 = cmpfit.prop_u.Pmax;
      td4 = cmpfit.prop_w.Pmax;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Maximum pressure Pmax                 : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }
    else if (! cmpfit.prop_use.Pmax){
      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "INFO: Maximum pressure Pmax           : ";
      print_prop_pred(fout, propstr, cmpfit.prop_pred.Pmax);
    }


    if (cmpfit.prop_use.displmax){
      td1 = cmpfit.prop_pred.displmax;
      td2 = cmpfit.prop_readin.displmax;
      tb1 = cmpfit.use_u.displmax;
      tb2 = cmpfit.use_w.displmax;
      td3 = cmpfit.prop_u.displmax;
      td4 = cmpfit.prop_w.displmax;

      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "Maximum displacement displmax         : ";
      print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
    }
    else if (! cmpfit.prop_use.displmax){
      std::ostringstream sstream;
      sstream << format("%15s") % cmpfit.name;
      compstr = "Compound: " + sstream.str() + " : ";

      propstr = compstr + "INFO: Maximum displacement displmax   : ";
      print_prop_pred(fout, propstr, cmpfit.prop_pred.displmax);
    }




    // Forces !!!
    if (cmpfit.prop_use.frc){
      tb1 = cmpfit.use_u.frc;
      tb2 = cmpfit.use_w.frc;
      td3 = -1.0;
      td4 = -1.0;

      int nb = cmpfit.basis_elems.size();
      for (int i=0; i<nb; ++i){
	for (int k=0; k<3; ++k){
	  td1 = cmpfit.prop_pred.frc[i][k];
	  td2 = cmpfit.prop_readin.frc[i][k];
	  if (tb1) td3 = cmpfit.prop_u.frc[i][k];
	  else     td4 = cmpfit.prop_w.frc[i][k];

	  std::ostringstream sstream;
	  sstream << format("%15s") % cmpfit.name;
	  compstr = "Compound: " + sstream.str() + " : ";

	  propstr = compstr + "Force component  iat " + tostring(nb) + " direction " + tostring(k+1) + ")       : ";
	  print_prop_readin_pred_comp(fout, firsttime, tb1, tb2, propstr, td1, td2, td3, td4);
	}
      }
    }







    fout << "..............................................................................." << std::endl;

  }

  if (DX.size()>0)
    fout << "-------------------------------------------------------------------------------" << std::endl;

}




void print_prop_readin_pred_comp(std::ostream & fout,
				 bool firsttime,
				 bool tb1, bool tb2,
				 std::string propstr,
				 double td1, double td2,
				 double td3, double td4){

  const double llim=1.0e-4;
  const double ulim=1.0e+4;
  std::string formatf = "%20.10f";
  std::string formate = "%20.10e";

  fout << propstr;
  fout << "read-in  ";
  if (abs(td2)<llim || abs(td2)>ulim) fout << format(formate) % td2;
  else                                fout << format(formatf) % td2;

  // -----------------------------------------------------------------------
  if (!firsttime){
    fout << "  predicted  ";
    if (abs(td1)<llim || abs(td1)>ulim) fout << format(formate) % td1;
    else                                fout << format(formatf) % td1;

    fout << "  rel. change  " << format("%10.3e")  % fp_divide(td1-td2, td2);
  }
  // -----------------------------------------------------------------------

  if (tb1){
    fout << "  uncertainty  ";
    if (abs(td3)<llim || abs(td3)>ulim) fout << format(formate) % td3;
    else                                fout << format(formatf) % td3;
  }
  else {
    fout << "  weight       ";
    if (abs(td4)<llim || abs(td4)>ulim) fout << format(formate) % td4;
    else                                fout << format(formatf) % td4;
  }
  fout << std::endl;
}


void print_prop_pred(std::ostream & fout,
		     std::string propstr,
		     double td1){

  const double llim=1.0e-4;
  const double ulim=1.0e+4;
  std::string formatf = "%20.10f";
  std::string formate = "%20.10e";

  fout << propstr;
  if (abs(td1)<llim || abs(td1)>ulim) fout << format(formate) % td1;
  else                                fout << format(formatf) % td1;
  fout << std::endl;
}




