
#ifndef REPORT_HPP
#define REPORT_HPP

#include <fstream>

#include "chisq-basics.hpp"
#include "funcfit-basics.hpp"
#include "funcfit-conjgrad.hpp"
#include "funcfit-errors.hpp"
#include "funcfit-ls-gauss-newton.hpp"
#include "funcfit-ls-leve-marq.hpp"
#include "funcfit-ls-powelldogleg.hpp"
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
#include "compoundfit.hpp"
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


using namespace std;
using namespace utils;


void report_pot_prop(ParamPot & param,
		     Vector<CompoundStructureFit> & DX,
		     Vector<double> & DY,
		     Vector<double> & MDY
		     );

void report_pot_prop_ext(ParamPot & param,
			 Vector<CompoundStructureFit> & DX,
			 Vector<double> & DY,
			 Vector<double> & MDY,
			 ostream & fout_pot=cout,
			 ostream & fout_prop=cout
			 );

void report_pot(PotentialInformationFit * p_potinfo,
		bool fittable_ones=true,
		bool fixed_ones=false,
		ostream & fout_pot=cout
		);

void report_prop(Vector<CompoundStructureFit> & DX,
		 ostream & fout_prop=cout,
		 bool firsttime=false
		 );


void print_prop_readin_pred_comp(ostream & fout,
				 bool firsttime,
				 bool tb1, bool tb2,
				 string propstr,
				 double td1, double td2,
				 double td3, double td4);



#endif

