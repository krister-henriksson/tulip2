


#include <string>
#include <new>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-errors.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "potclasses.hpp"
#include "param.hpp"


using namespace utils;




// #########################################################################################################
// #########################################################################################################



CutoffScreeningPair::CutoffScreeningPair(){
  tersoff    =false;
  perriot_cut=false;
  perriot_scr=false;
  R=D=0;
  pn=pm=prcut=prmin=prmax=0;
  rcut=0.0;
}



// #########################################################################################################
// #########################################################################################################


CutoffScreening::CutoffScreening(){
  mode = "none";
  name = "none";
  use = "false";

  parval.resize(0);
  parname.resize(0);
}


void CutoffScreening::set_name_mode(std::string na, std::string mo){
  mode = "none";
  name = "none";
  use = "false";
  parval.resize(0);
  parname.resize(0);

  if      (na=="tersoff"){
    name="tersoff";
    mode="cut";
    use=true;
    parname.resize(2);
    parval.resize(2);
    parname[0]="R";
    parname[1]="D";
    parval[0]=parval[1]=0.0;
  }
  else if (na=="perriot" && mo=="cut"){
    name="perriot";
    mode="cut";
    use=true;
    parname.resize(2);
    parval.resize(2);
    parname[0]="prmin";
    parname[1]="prmax";
    parval[0]=parval[1]=0.0;
  }
  else if (na=="perriot" && mo=="scr"){
    name="perriot";
    mode="scr";
    use=true;
    parname.resize(5);
    parval.resize(5);
    parname[0]="pn";
    parname[1]="pm";
    parname[2]="prcut";
    parname[3]="prmin";
    parname[4]="prmax";
    parval[0]=parval[1]=parval[2]=parval[3]=parval[4]=0.0;
  }
  else {
    aborterror("ERROR: Unknown combination of cutoff/screening name and mode: " + na + " " + mode + ". Exiting.");
  }
}


void   CutoffScreening::set_parval(std::string name, double pv){
  for (int i=0; i<parname.size(); ++i)
    if (parname[i]==name){ parval[i]=pv; return; }
}
double CutoffScreening::get_parval(std::string name){
  for (int i=0; i<parname.size(); ++i)
    if (parname[i]==name) return parval[i];
  //aborterror("ERROR: Parameter name " + name + " not found in the CutoffScreening sub-object. Exiting.");
}


double CutoffScreening::rcut(void){
  if      (name=="tersoff") return get_parval("R")+get_parval("D");
  else if (name=="perriot" && mode=="cut") return get_parval("prmax");
  else if (name=="perriot" && mode=="scr") return get_parval("prcut");
  else return 0.0;
}


void CutoffScreening::init_lims(void){
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
void   CutoffScreening::set_par_extr(std::string name, double pv, std::string minmax){
  for (int i=0; i<parname.size(); ++i){
    if (parname[i]==name){
      if      (minmax=="min"){ parmin[i]=pv; return; }
      else if (minmax=="max"){ parmax[i]=pv; return; }
    }
  }
  //aborterror("ERROR: Parameter name " + name + " not found in the CutoffScreening sub-object. Exiting.");
}
double CutoffScreening::get_par_extr(std::string name, std::string minmax){
  for (int i=0; i<parname.size(); ++i){
    if (parname[i]==name){
      if      (minmax=="min") return parmin[i];
      else if (minmax=="max") return parmax[i];
    }
  }
  //aborterror("ERROR: Parameter name " + name + " not found in the CutoffScreening sub-object. Exiting.");
}

int           CutoffScreening::npar(){
  return parname.size();
}
void          CutoffScreening::set_par_types(){
  for (int i=0; i<parname.size(); ++i)
    set_param_type(parmin[i], parmax[i], partype[i]);
}
void          CutoffScreening::set_par_type(std::string name){
  for (int i=0; i<parname.size(); ++i)
    if (parname[i]==name){
      set_param_type(parmin[i], parmax[i], partype[i]);
      return;
    }
  //aborterror("ERROR: Parameter name " + name + " not found in the CutoffScreening sub-object. Exiting.");
}
parametertype CutoffScreening::get_par_type(std::string name){
  for (int i=0; i<parname.size(); ++i)
    if (parname[i]==name) return partype[i];
  //aborterror("ERROR: Parameter name " + name + " not found in the CutoffScreening sub-object. Exiting.");
}



// #########################################################################################################
// #########################################################################################################

Potential_ABOPPair::Potential_ABOPPair(){
  elemname1="none";
  elemname2="none";
  D0=r0=beta=S=p=gamma=c=d=h=rfermi=bfermi=0;
}


Potential_ABOP::Potential_ABOP()
  :
  elemname1("none"), elemname2("none"),
  maxindex(10)
{
  /*
  use_cutoff_only = false;
  use_screening   = false;
  rcut_fun = "none";
  rcut_scr = "none";
  */

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

  for (int i=0; i<=maxindex; ++i) parval[i]=0.0;

  // Defaults:
  parval[10] = 0.5;

  /*
  parname[11]="D";  // Tersoff cutoff parameter
  parname[12]="R";  // Tersoff cutoff parameter

  parname[13]="pn";  // Perriot cutoff parameter
  parname[14]="pm";  // Perriot cutoff parameter
  parname[15]="prcut"; // Perriot cutoff parameter
  parname[16]="prmin"; // Perriot cutoff parameter
  parname[17]="prmax"; // Perriot cutoff parameter

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
  */
}


void   Potential_ABOP::set_parval(std::string name, double pv){
  int n=parname.size();
  for (int i=0; i<n; ++i){
    if (parname[i]==name){
      parval[i]=pv; return;
    }
  }
  // Try the Cutoff/Screening sub-object:
  rcs.set_parval(name, pv);
}
double Potential_ABOP::get_parval(std::string name){
  int n=parname.size();
  for (int i=0; i<n; ++i){
    if (parname[i]==name){
      return parval[i];
    }
  }
  // Try the Cutoff/Screening sub-object:
  return rcs.get_parval(name);
}



double Potential_ABOP::rcut(void){
  return rcs.rcut();
}
/*
  if      (rcut_fun=="tersoff")
    return parname2val("D") + parname2val("R");
  else if (rcut_fun=="perriot")
    return parname2val("prmax");
  else if (rcut_scr=="perriot")
    return parname2val("prcut");
  else
    return 0.0;
*/




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
  rcs.init_lims();
}
void   Potential_ABOP::set_par_extr(std::string name, double pv, std::string minmax){
  for (int i=0; i<parname.size(); ++i){
    if (parname[i]==name){
      if      (minmax=="min"){ parmin[i]=pv; return; }
      else if (minmax=="max"){ parmax[i]=pv; return; }
    }
  }
  rcs.set_par_extr(name, pv, minmax);
}
double Potential_ABOP::get_par_extr(std::string name, std::string minmax){
  for (int i=0; i<parname.size(); ++i){
    if (parname[i]==name){
      if      (minmax=="min") return parmin[i];
      else if (minmax=="max") return parmax[i];
    }
  }
  return rcs.get_par_extr(name, minmax);
}

int           Potential_ABOP::npar(){
  int n1 = parname.size();
  int n2 = rcs.npar();
  return n1+n2;
}
void          Potential_ABOP::set_par_types(){
  for (int i=0; i<parname.size(); ++i)
    set_param_type(parmin[i], parmax[i], partype[i]);
  rcs.set_par_types();
}
void          Potential_ABOP::set_par_type(std::string name){
  for (int i=0; i<parname.size(); ++i)
    if (parname[i]==name){
      set_param_type(parmin[i], parmax[i], partype[i]);
      return;
    }
  rcs.set_par_type(name);
}
parametertype Potential_ABOP::get_par_type(std::string name){
  for (int i=0; i<parname.size(); ++i)
    if (parname[i]==name) return partype[i];
  return rcs.get_par_type(name);
}

