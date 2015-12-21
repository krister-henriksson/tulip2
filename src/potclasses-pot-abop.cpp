


#include <string>
#include <new>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "potclasses.hpp"
#include "param.hpp"


using namespace utils;



Potential_ABOP::Potential_ABOP()
  :
  elemname1("none"), elemname2("none"),
  maxindex(17)
{

  rcut_fun = "none";

  parval.resize(maxindex+1);
  parname.resize(maxindex+1);

  parname[0]="D0";
  parname[1]="r0";
  parname[2]="beta";
  parname[3]="S";
  parname[4]="gamma";
  parname[5]="c";
  parname[6]="d";
  parname[7]="h";
  parname[8]="bfermi";
  parname[9]="rfermi";
  parname[10]="p";

  parname[11]="D";  // Tersoff cutoff parameter
  parname[12]="R";  // Tersoff cutoff parameter

  parname[13]="pn";  // Perriot cutoff parameter
  parname[14]="pm";  // Perriot cutoff parameter
  parname[15]="prcut"; // Perriot cutoff parameter
  parname[16]="prmin"; // Perriot cutoff parameter
  parname[17]="prmax"; // Perriot cutoff parameter

  for (int i=0; i<=maxindex; ++i) parval[i]=0.0;
  parval[8] = 10.0;
  parval[9] =  1.0;
  parval[10] = 0.5;

  parval[11] = 0.1;
  parval[12] = 5.0;

  parval[13] =  5;
  parval[14] = 48;
  parval[15] = 5.0;
  parval[16] = 4.8;
  parval[17] = 5.0;

 
}



int Potential_ABOP::parname2idx(std::string name){
  for (int i=0; i<parname.size(); ++i){
    if (name == parname[i]) return i;
  }
  return -1;
}

std::string Potential_ABOP::paridx2name(int idx){
  if (idx<0 || idx>maxindex) return "none";
  else return parname[idx];
}

double & Potential_ABOP::parname2val(std::string name){
  int i = parname2idx(name);
  return parval[i];
}

double Potential_ABOP::rcut(void){
  if      (rcut_fun=="tersoff")
    return parname2val("D") + parname2val("R");
  else if (rcut_fun=="perriot")
    return parname2val("prcut");
  else
    return 0.0;
}



void Potential_ABOP::init_lims(void){
  int n = parname.size();
  partype.resize(n);
  parmin.resize(n);
  parmax.resize(n);
  for (int i=0; i<n; ++i){
    partype[i] = PARAM_FIXED;
    parmin[i] = 1;
    parmax[i] = 1;
  }
}



