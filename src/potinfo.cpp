


#include <iostream>
#include <string>

#include <boost/format.hpp>

#include <cstdio>
#include <cstdlib>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "utils-errors.hpp"

#include "utils-string.hpp"

#include "elem-iacs.hpp"
#include "specs-fit-prop-pot.hpp"
#include "potinfo.hpp"
#include "potclasses.hpp"
#include "compound.hpp"

#include "omp-basics.hpp"


using namespace utils;
using boost::format;



PotentialInformation::PotentialInformation() { }
// call read_info(...) manually


PotentialInformation::PotentialInformation(std::string filename_info){
  read_info(filename_info);
}




// ###############################################################################

// For any pair of elements, report the base potential used:
std::string PotentialInformation::basepot(std::string s1, std::string s2){
  int i1 = elem.name2idx(s1);
  int i2 = elem.name2idx(s2);
  if (i1<=i2) return mbasepot.elem(i1, i2);
  else        return mbasepot.elem(i2, i1);
}
std::string PotentialInformation::basepot(int i1, int i2){
  if (i1<=i2) return mbasepot.elem(i1, i2);
  else        return mbasepot.elem(i2, i1);
}


// For any pair of elements, report the index of the potential vector
// where the parametrization is stored. Base potential type is
// also needed, to determine which potential vector to select.
int PotentialInformation::basepot_vecidx(std::string s1, std::string s2){
  int i1 = elem.name2idx(s1);
  int i2 = elem.name2idx(s2);
  if (i1<=i2) return mbasepot_vecidx.elem(i1, i2);
  else        return mbasepot_vecidx.elem(i2, i1);
}
int PotentialInformation::basepot_vecidx(int i1, int i2){
  if (i1<=i2) return mbasepot_vecidx.elem(i1, i2);
  else        return mbasepot_vecidx.elem(i2, i1);
}


// ###############################################################################

// Record if reppot is used or not:
bool & PotentialInformation::use_reppot(std::string s1, std::string s2){
  int i1 = elem.name2idx(s1);
  int i2 = elem.name2idx(s2);
  return PotentialInformation::use_reppot(i1, i2);
}

bool & PotentialInformation::use_reppot(int i1, int i2){
  int k, res=-1;
  bool done=false;
  if (i1>i2){
    k = i2;
    i2 = i1;
    i1 = k;
  }
  k=0;
  int nel = elem.nelem();
  for (int i=0; i<nel; ++i){
    for (int j=i; j<nel; ++j){
      if (i==i1 && j==i2){
	res = k; done = true;
      }
      if (done) break;
      k++;
    }
    if (done) break;
  }
  if (res==-1)
    aborterror("Error: Did not find reppot use indication for pair "
	       + elem.idx2name(i1) + "-"
	       + elem.idx2name(i2) + ". Exiting.");

  return muse_reppot[res];
}

void PotentialInformation::reppot_finalize(void){
  int ipair=0, k=0, nel = elem.nelem();
  for (int i=0; i<nel; ++i){
    for (int j=i; j<nel; ++j){
      if (muse_reppot[k]){
	  mreppot_vecidx.elem(i,j) = ipair++;
	  pot_Reppot.resize( ipair );
      }
      k++;
    }
  }
}

int PotentialInformation::reppot_vecidx(std::string s1, std::string s2){
  if (!use_reppot(s1,s2))
    return -1;

  int i1 = elem.name2idx(s1);
  int i2 = elem.name2idx(s2);
  if (i1<=i2) return mreppot_vecidx.elem(i1, i2);
  else        return mreppot_vecidx.elem(i2, i1);
}

int PotentialInformation::reppot_vecidx(int i1, int i2){
  if (!use_reppot(i1,i2))
    return -1;

  if (i1<=i2) return mreppot_vecidx.elem(i1, i2);
  else        return mreppot_vecidx.elem(i2, i1);
}

// ###############################################################################
// ###############################################################################
// ###############################################################################
// ###############################################################################


double PotentialInformation::get_abop_omega(std::string s1, std::string s2, std::string s3){
  int i = elem.name2idx(s1);
  int j = elem.name2idx(s2);
  int k = elem.name2idx(s3);
  if (! use_abop_omega.elem(i,j,k)){
    aborterror("ERROR: Tried to get ABOP omega(" + s1 + " " + s2 + " " + s3
	       + ", but omega is not specified as being used! Exiting.");
  }
  return mabop_omega.elem(i,j,k);
}

double PotentialInformation::get_abop_omega(int i, int j, int k){
  if (! use_abop_omega.elem(i,j,k)){
    aborterror("ERROR: Tried to get ABOP omega(" + tostring(i) + " " + tostring(j) + " " + tostring(k)
	       + ", but omega is not specified as being used! Exiting.");
  }
  return mabop_omega.elem(i,j,k);
}



void PotentialInformation::set_abop_omega(std::string s1, std::string s2, std::string s3, double val){
  int i = elem.name2idx(s1);
  int j = elem.name2idx(s2);
  int k = elem.name2idx(s3);
  if (! use_abop_omega.elem(i,j,k)){
    aborterror("ERROR: Tried to set ABOP omega(" + s1 + " " + s2 + " " + s3
	       + ", but omega is not specified as being used! Exiting.");
  }
  mabop_omega.elem(i,j,k) = val;
}

void PotentialInformation::set_abop_omega(int i, int j, int k, double val){
  if (! use_abop_omega.elem(i,j,k)){
    aborterror("ERROR: Tried to set ABOP omega(" + tostring(i) + " " + tostring(j) + " " + tostring(k)
	       + ", but omega is not specified as being used! Exiting.");
  }
  mabop_omega.elem(i,j,k) = val;
}


// ###############################################################################
// ###############################################################################
// ###############################################################################
// ###############################################################################




double PotentialInformation::get_rcut_max(void){
  int i,j,ipair,p, nel = elem.nelem();
  double rc=0.0, rcut_max=0.0;
  
  p=0;
  for (i=0; i<nel; ++i){
    for (j=i; j<nel; ++j){
      /*
      std::string s1 = elem.idx2name(i);
      std::string s2 = elem.idx2name(j);

      if (basepot(s1,s2) == "EAM"){
	ipair = basepot_vecidx(s1,s2);
	rc    = pot_EAM[ipair].rcut();
      }
      else if (basepot(s1,s2) == "ABOP"){
	ipair = basepot_vecidx(s1,s2);
	rc    = pot_ABOP[ipair].rcut();
      }
      */

      if (basepot(i,j) == "EAM"){
	ipair = basepot_vecidx(i,j);
	rc    = pot_EAM[ipair].rcut();
      }
      else if (basepot(i,j) == "ABOP"){
	ipair = basepot_vecidx(i,j);
	rc    = pot_ABOP[ipair].rcut();
      }


      if (p==0 || (p>0 && rc>rcut_max))
	rcut_max=rc;
      p++;
    }
  }

  return rcut_max;
}


double PotentialInformation::get_rcut_max(Vector<std::string> elemnames){
  int i,j,ipair,p;
  double rc=0.0, rcut_max=0.0;

  p=0;
  for (i=0; i<elemnames.size(); ++i){
    for (j=i; j<elemnames.size(); ++j){

      std::string s1 = elemnames[i];
      std::string s2 = elemnames[j];

      if (basepot(s1,s2) == "EAM"){
	ipair = basepot_vecidx(s1,s2);
	rc    = pot_EAM[ipair].rcut();
      }
      else if (basepot(s1,s2) == "ABOP"){
	ipair = basepot_vecidx(s1,s2);
	rc    = pot_ABOP[ipair].rcut();
      }

      if (p==0 || (p>0 && rc>rcut_max))
	rcut_max=rc;
      p++;
    }
  }

  return rcut_max;
}






// ###############################################################################
// ###############################################################################


PotentialInformationFit::PotentialInformationFit(){ }


PotentialInformationFit::PotentialInformationFit(std::string filename_info,
						 std::string filename_specs
						 )
  : PotentialInformation::PotentialInformation(filename_info)
{
  read_info_fit(filename_info);
  read_specs(filename_specs);
}


// ###############################################################################

bool & PotentialInformationFit::is_fittable(std::string s1, std::string s2){
  int i1 = elem.name2idx(s1);
  int i2 = elem.name2idx(s2);
  int k, res=-1;
  bool done=false;
  if (i1>i2){
    k = i2;
    i2 = i1;
    i1 = k;
  }
  k=0;
  for (int i=0; i<elem.nelem(); ++i){
    for (int j=i; j<elem.nelem(); ++j){
      if (i==i1 && j==i2){
	res = k; done = true;
	return mfit[k];
      }
      if (done) break;
      k++;
    }
    if (done) break;
  }
  if (res==-1)
    aborterror("Error: Did not find fit use indication for pair " + s1 + "-" + s2 + ". Exiting.");

  return mfit[res];
}


bool & PotentialInformationFit::is_fittable(int i1, int i2){
  int k, res=-1;
  bool done=false;
  if (i1>i2){
    k = i2;
    i2 = i1;
    i1 = k;
  }
  k=0;
  int nel = elem.nelem();
  for (int i=0; i<nel; ++i){
    for (int j=i; j<nel; ++j){
      if (i==i1 && j==i2){
	res = k; done = true;
	return mfit[k];
      }
      if (done) break;
      k++;
    }
    if (done) break;
  }
  if (res==-1)
    aborterror("Error: Did not find fit use indication for pair "
	       + elem.idx2name(i1) + "-"
	       + elem.idx2name(i2) + ". Exiting.");

  return mfit[res];
}

