


#ifndef COMPOUND_FIT_HPP
#define COMPOUND_FIT_HPP



#include <iostream>
#include <string>

#include <cstdio>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"

#include "elem-iacs.hpp"
#include "specs-fit-prop-pot.hpp"
#include "mdsettings.hpp"
//#include "mdsystem.hpp"
//#include "param-pot.hpp"


#include "compound.hpp"


using namespace std;
using namespace utils;




class CompoundStructureFit
  :
  public CompoundStructure
{
public:
  CompoundPropertiesUse              prop_use;
  CompoundProperties                 prop_readin;
  CompoundProperties                 prop_pred;

  MDSettings                         mds_specs;

  CompoundPropertiesUseUncertainties use_u;
  CompoundPropertiesUseWeights       use_w;
  CompoundPropertiesUncertainties    prop_u;
  CompoundPropertiesWeights          prop_w;


  CompoundStructureFit();
  // CompoundStructureFit(const CompoundStructureFit & sv);
  // CompoundStructureFit & operator=(const CompoundStructureFit & sv);

  void read_forces(void);

  void check_and_fix_uses();
  int NData();

  void getprop(ParamPot                  & param);
  void get_B_Bp(MDSystem                 & mds,
		ParamPot                 & param,
		Vector< Vector<double> > & pos_bak,
		Matrix<double>           & boxdir_bak,
		Vector<double>           & boxlen_bak,
		double                   & E0,
		double                   & V0);
  void get_Cij(MDSystem                  & mds,
	       ParamPot                  & param,
	       Vector< Vector<double> >  & pos_bak,
	       Matrix<double>            & boxdir_bak,
	       Vector<double>            & boxlen_bak,
	       double                    & E0,
	       double                    & V0);


} ;
  



// #####################################################################
// #####################################################################




class CompoundListFit {
public:
  Elements elem;
  int ncompounds;
  Vector<CompoundStructureFit> compounds;

  // Read structures:
  CompoundListFit(const Elements & el,
		  MDSettings     & mds_specs,
		  string         filename);
  int NData();
} ;


#endif
