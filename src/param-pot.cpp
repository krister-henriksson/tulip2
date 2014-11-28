


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"
#include "elem-iacs.hpp"
#include "potinfo.hpp"
#include "param.hpp"
#include "param-pot.hpp"
#include "exiterrors.hpp"


using namespace exiterrors;


// construct
ParamPot::ParamPot(){
  p_potinfo = 0;
}


// destroy
ParamPot::~ParamPot(){
  p_potinfo = 0;
}


// copy
ParamPot::ParamPot(const ParamPot & sv)
  : Param(sv)
{
  p_potinfo = sv.p_potinfo;
  update_par();
}


// assign
ParamPot & ParamPot::operator= (const ParamPot & sv){
  if (this==&sv) return *this;
  
  // 1. Call base class assignment operator:
  Param::operator= (sv);

  // 2. Handle members in derived class:
  p_potinfo = sv.p_potinfo;
  update_par();

  return *this;
}


// alt. construct
ParamPot::ParamPot(PotentialInformationFit * p_pi){
  p_potinfo = p_pi;
  update_par();
}




// ######################################################################
// Insert fitting parameters into the free potential parameters
// ######################################################################
  
void ParamPot::update_pot(){
  if (p_potinfo==0){
    aborterror("Error: update_pot(): Potentialinformation is not set. Exiting.");
  }

  cout << "Updating potential from parameters ..." << endl;


  int ipar=0,j,k, nel = (*p_potinfo).elem.nelem();

  for (int i1=0; i1<nel; ++i1){
    for (int i2=i1; i2<nel; ++i2){

      j = (*p_potinfo).basepot_vecidx(i1,i2);
      /*
      string s1 = (*p_potinfo).elem.idx2name(i1);
      string s2 = (*p_potinfo).elem.idx2name(i2);

      j = (*p_potinfo).basepot_vecidx(s1,s2);
      */
      if (j<0) continue;

      if (! (*p_potinfo).is_fittable(i1,i2)) continue;
      //if (! (*p_potinfo).is_fittable(s1,s2)) continue;

      if ((*p_potinfo).basepot(i1,i2) == "ABOP"){
	//if ((*p_potinfo).basepot(s1,s2) == "ABOP"){
	for (k=0; k<(*p_potinfo).pot_ABOP[j].partype.size(); ++k)
	  if ((*p_potinfo).pot_ABOP[j].partype[k] != PARAM_FIXED) (*p_potinfo).pot_ABOP[j].parval[k] = X(ipar++);
      }
    }
  }

  for (int i1=0; i1<nel; ++i1){
    for (int i2=0; i2<nel; ++i2){
      for (int i3=0; i3<nel; ++i3){

	if ((*p_potinfo).use_abop_alpha.elem(i1,i2,i3)){
	  if ((*p_potinfo).abop_alpha_partype.elem(i1,i2,i3) != PARAM_FIXED)
	    (*p_potinfo).abop_alpha.elem(i1,i2,i3) = X(ipar++);
	}

      }
    }
  }

  for (int i1=0; i1<(*p_potinfo).elem.nelem(); ++i1){
    for (int i2=0; i2<(*p_potinfo).elem.nelem(); ++i2){
      for (int i3=0; i3<(*p_potinfo).elem.nelem(); ++i3){

	if ((*p_potinfo).use_abop_omega.elem(i1,i2,i3)){
	  if ((*p_potinfo).abop_omega_partype.elem(i1,i2,i3) != PARAM_FIXED)
	    (*p_potinfo).set_abop_omega(i1,i2,i3, X(ipar++));
	}
	
      }
    }
  }

  for (int i1=0; i1<(*p_potinfo).elem.nelem(); ++i1){
    for (int i2=0; i2<(*p_potinfo).elem.nelem(); ++i2){

      if ((*p_potinfo).use_abop_2mu.elem(i1,i2)){

	if ((*p_potinfo).abop_2mu_partype.elem(i1,i2) != PARAM_FIXED)
	  (*p_potinfo).abop_2mu.elem(i1,i2) = X(ipar++);
      }
      
    }
  }


  
}



// ######################################################################
// Insert free potential parameters into fitting parameters
// ######################################################################

void ParamPot::update_par(){
  if (p_potinfo==0){
    aborterror("Error: update_par(): Potentialinformation is not set. Exiting.");
  }

  int nel = (*p_potinfo).elem.nelem();

  int ipar=0,j,k;
  Vector<double> xi, ximin, ximax;
  Vector<parametertype> xitype;

  xi.resize(0);
  ximin.resize(0);
  ximax.resize(0);
  xitype.resize(0);


  cout << "Updating parameters from potential ..." << endl;

  /*
  cout << "PARAM_FIXED " << PARAM_FIXED << endl;
  cout << "PARAM_FREE " << PARAM_FREE << endl;
  cout << "PARAM_FREE_WITH_LIMITS " << PARAM_FREE_WITH_LIMITS << endl;
  */

  for (int i1=0; i1<nel; ++i1){
    string s1 = (*p_potinfo).elem.idx2name(i1);
    for (int i2=i1; i2<nel; ++i2){
      string s2 = (*p_potinfo).elem.idx2name(i2);

      j = (*p_potinfo).basepot_vecidx(s1,s2);
      if (j<0) continue;

      if (! (*p_potinfo).is_fittable(s1,s2)) continue;


      if ((*p_potinfo).basepot(s1,s2) == "ABOP"){

	for (k=0; k<(*p_potinfo).pot_ABOP[j].partype.size(); ++k)
	  if ((*p_potinfo).pot_ABOP[j].partype[k] != PARAM_FIXED){
	    ++ipar;
	    xi.push_back(   (*p_potinfo).pot_ABOP[j].parval[k]);
	    xitype.push_back( (*p_potinfo).pot_ABOP[j].partype[k] );
	    ximin.push_back((*p_potinfo).pot_ABOP[j].parmin[k]);
	    ximax.push_back((*p_potinfo).pot_ABOP[j].parmax[k]);
	  }

      }

    }
  }


  for (int i1=0; i1<nel; ++i1){
    string s1 = (*p_potinfo).elem.idx2name(i1);
    for (int i2=0; i2<nel; ++i2){
      string s2 = (*p_potinfo).elem.idx2name(i2);
      for (int i3=0; i3<nel; ++i3){
	string s3 = (*p_potinfo).elem.idx2name(i3);

	if ((*p_potinfo).use_abop_alpha.elem(i1,i2,i3)){

	  if ((*p_potinfo).abop_alpha_partype.elem(i1,i2,i3) != PARAM_FIXED){
	    ++ipar;
	    xi.push_back(   (*p_potinfo).abop_alpha.elem(i1,i2,i3));
	    ximin.push_back((*p_potinfo).abop_alpha_parmin.elem(i1,i2,i3));
	    ximax.push_back((*p_potinfo).abop_alpha_parmax.elem(i1,i2,i3));
	    xitype.push_back( (*p_potinfo).abop_alpha_partype.elem(i1,i2,i3) );
	  }

	}

      }
    }
  }


  for (int i1=0; i1<nel; ++i1){
    string s1 = (*p_potinfo).elem.idx2name(i1);
    for (int i2=0; i2<nel; ++i2){
      string s2 = (*p_potinfo).elem.idx2name(i2);
      for (int i3=0; i3<nel; ++i3){
	string s3 = (*p_potinfo).elem.idx2name(i3);

	if ((*p_potinfo).use_abop_omega.elem(i1,i2,i3)){

	  if ((*p_potinfo).abop_omega_partype.elem(i1,i2,i3) != PARAM_FIXED){
	    ++ipar;
	    xi.push_back(   (*p_potinfo).get_abop_omega(s1,s2,s3));
	    ximin.push_back((*p_potinfo).abop_omega_parmin.elem(i1,i2,i3));
	    ximax.push_back((*p_potinfo).abop_omega_parmax.elem(i1,i2,i3));
	    xitype.push_back( (*p_potinfo).abop_omega_partype.elem(i1,i2,i3) );
	  }
	}

      }
    }
  }


  for (int i1=0; i1<nel; ++i1){
    string s1 = (*p_potinfo).elem.idx2name(i1);
    for (int i2=0; i2<nel; ++i2){
      string s2 = (*p_potinfo).elem.idx2name(i2);

      if ((*p_potinfo).use_abop_2mu.elem(i1,i2)){

	if ((*p_potinfo).abop_2mu_partype.elem(i1,i2) != PARAM_FIXED){
	  ++ipar;
	  xi.push_back(   (*p_potinfo).abop_2mu.elem(i1,i2));
	  ximin.push_back((*p_potinfo).abop_2mu_parmin.elem(i1,i2));
	  ximax.push_back((*p_potinfo).abop_2mu_parmax.elem(i1,i2));
	  xitype.push_back( (*p_potinfo).abop_2mu_partype.elem(i1,i2) );
	}
	
      }
      
    }
  }
  



  X()    = xi;
  Xmin() = ximin;
  Xmax() = ximax;
  Xtype() = xitype;


}





// Given any vector (of correct type of course), update all the parameters
// (the full parameter vector).
utils::Vector<double> ParamPot::Xupdate(const utils::Vector<double> & xi){
  Param::Xupdate(xi);
  update_pot();
  return X();
}


