
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/format.hpp>

#include <cmath>

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
#include "errors.hpp"




#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;





using namespace utils;
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



Vector<double> get_comp_prop(ParamPot & param, Vector<CompoundStructureFit> & DX){

  Vector<double> MDY;
  CompoundStructureFit cmpfit;
  int iDX,sizeDX,k,p;
  bad_point err_bad_point;

  // Make sure potentials have been updated from the fitting parameters:
  param.update_pot();

  sizeDX = DX.size();

  // ###################################################################
  // Get properties of fittable compounds
  // ###################################################################
  double Ecoh_delta_ref = 0.0;
  for (iDX=0; iDX<sizeDX; ++iDX){
    // -------------------------------------------------------------------
    // Make a local copy which is easier to work with:
    // -------------------------------------------------------------------
    cmpfit = DX[iDX];

    std::cout << "Getting properties of compound " << cmpfit.name << " ..." << std::endl;
    try {
      cmpfit.getprop(param);
    }
    catch (bad_mds & err_bad_mds){
      std::cout << "ERROR: Bad MDS run detected!" << std::endl;
      throw err_bad_point;
    }

    if (cmpfit.Ecoh_delta_refcomp){
      Ecoh_delta_ref = cmpfit.prop_pred.Ecoh;
      // std::cout << "Ecoh_delta_ref = " << Ecoh_delta_ref << " from compound " << cmpfit.name << std::endl;
    }

    // -------------------------------------------------------------------
    // Make sure we get local updates back to calling function:
    // -------------------------------------------------------------------
    DX[iDX] = cmpfit;
  }



  // ###################################################################
  // Fill the ModelDataY vector
  // ###################################################################
  MDY.resize(0);
  for (iDX=0; iDX<sizeDX; ++iDX){
    // -------------------------------------------------------------------
    // Make a local copy which is easier to work with:
    // -------------------------------------------------------------------
    cmpfit = DX[iDX];

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

    if (cmpfit.prop_use.Ecoh_delta){
      // std::cout << "Using Ecoh_delta with preloaded value " << cmpfit.prop_pred.Ecoh_delta;
      cmpfit.prop_pred.Ecoh_delta -= Ecoh_delta_ref;
      MDY.push_back(cmpfit.prop_pred.Ecoh_delta);
      /*
	std::cout << " and subtracting a term " << Ecoh_delta_ref
	<< " so that predicted value is " << cmpfit.prop_pred.Ecoh_delta;
      */
    }


    if (cmpfit.prop_use.Emix) MDY.push_back(cmpfit.prop_pred.Emix);
    if (cmpfit.prop_use.B) MDY.push_back(cmpfit.prop_pred.B);
    if (cmpfit.prop_use.Bp) MDY.push_back(cmpfit.prop_pred.Bp);

    for (k=0; k<6; ++k)
      for (p=0; p<6; ++p)
	if (cmpfit.prop_use.C.elem(k,p)) MDY.push_back(cmpfit.prop_pred.C.elem(k,p));

    if (cmpfit.prop_use.Fmax) MDY.push_back(cmpfit.prop_pred.Fmax);
    if (cmpfit.prop_use.Pmax) MDY.push_back(cmpfit.prop_pred.Pmax);
    if (cmpfit.prop_use.displmax) MDY.push_back(cmpfit.prop_pred.displmax);

    if (cmpfit.prop_use.frc){
      int nb = cmpfit.basis_elems.size();
      for (int iat=0; iat<nb; ++iat)
	for (int k=0; k<3; ++k)
	  MDY.push_back(cmpfit.prop_pred.frc[iat][k]);
    }

    // -------------------------------------------------------------------
    // Make sure we get local updates back to calling function:
    // -------------------------------------------------------------------
    DX[iDX] = cmpfit;
  }




  // If e.g. last compound checked is not a reference compound, then we have
  // debugged forces of read-in compound, and we can quit the program.
  if ( ! DX[ DX.size()-1 ].mds_specs.is_ref_comp
       &&
       param.p_potinfo->specs_prop.mds_specs_common.debug_forces)
    aborterror("Warning: Debugging of forces done, performing quick and dirty exit.");
  

  return MDY;
}




