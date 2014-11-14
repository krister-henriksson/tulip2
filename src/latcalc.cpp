
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include <cmath>

#include "chisq-basics.hpp"
#include "exiterrors.hpp"
#include "funcfit-basics.hpp"
#include "funcfit-conjgrad.hpp"
#include "funcfit-exceptions.hpp"
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
#include "latcalc.hpp"


using namespace std;
using namespace utils;
using namespace exiterrors;
using namespace constants;
using namespace physconst;
using namespace funcfit;
using boost::format;


// ##############################################################################
// ##############################################################################
//
// Update potentials from fitting parameters, and then calculate properties
// of the used compounds.
//
// ##############################################################################
// ##############################################################################

Vector<double> latcalc(ParamPot & param, Vector<CompoundStructureFit> & DX){

  Vector<double> MDY;
  CompoundStructureFit cmpfit;
  int iDX,k,p;

  // Make sure potentials have been updated from the fitting parameters:
  param.update_pot();

  MDY.resize(0);
  for (iDX=0; iDX<DX.size(); ++iDX){
    cmpfit = DX[iDX];

    cout << "Getting properties of compound " << cmpfit.name << " ..." << endl;
    cmpfit.getprop(param);


    // ###################################################################
    // Make sure we get data back to calling function:
    // ###################################################################
    DX[iDX] = cmpfit;


    // ###################################################################
    // Fill the ModelDataY vector
    // ###################################################################
    if (cmpfit.prop_use.a) MDY.push_back(cmpfit.prop_pred.a);
    if (cmpfit.prop_use.b) MDY.push_back(cmpfit.prop_pred.b);
    if (cmpfit.prop_use.c) MDY.push_back(cmpfit.prop_pred.c);
    if (cmpfit.prop_use.bpa) MDY.push_back(cmpfit.prop_pred.bpa);
    if (cmpfit.prop_use.cpa) MDY.push_back(cmpfit.prop_pred.cpa);
    if (cmpfit.prop_use.r0) MDY.push_back(cmpfit.prop_pred.r0);
    if (cmpfit.prop_use.angle_ab) MDY.push_back(cmpfit.prop_pred.angle_ab);
    if (cmpfit.prop_use.angle_ac) MDY.push_back(cmpfit.prop_pred.angle_ac);
    if (cmpfit.prop_use.angle_bc) MDY.push_back(cmpfit.prop_pred.angle_bc);
    if (cmpfit.prop_use.Vatom) MDY.push_back(cmpfit.prop_pred.Vatom);
    if (cmpfit.prop_use.Ecoh) MDY.push_back(cmpfit.prop_pred.Ecoh);
    if (cmpfit.prop_use.Emix) MDY.push_back(cmpfit.prop_pred.Emix);
    if (cmpfit.prop_use.B) MDY.push_back(cmpfit.prop_pred.B);
    if (cmpfit.prop_use.Bp) MDY.push_back(cmpfit.prop_pred.Bp);
    for (k=0; k<6; ++k)
      for (p=0; p<6; ++p)
	if (cmpfit.prop_use.C.elem(k,p)) MDY.push_back(cmpfit.prop_pred.C.elem(k,p));
    if (cmpfit.prop_use.Fmax) MDY.push_back(cmpfit.prop_pred.Fmax);
    if (cmpfit.prop_use.Pmax) MDY.push_back(cmpfit.prop_pred.Pmax);
    if (cmpfit.prop_use.displmax) MDY.push_back(cmpfit.prop_pred.displmax);



  }


  // If e.g. last compound checked is not a reference compound, then we have
  // debugged forces of read-in compound, and we can quit the program.
  if ( ! DX[ DX.size()-1 ].mds_specs.is_ref_comp
       &&
       param.p_potinfo->specs_prop.mds_specs_common.debug_forces)
    aborterror("Warning: Debugging of forces done, performing quick and dirty exit.");
  

  return MDY;
}




