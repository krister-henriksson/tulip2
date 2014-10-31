


#include <string>

#include <boost/format.hpp>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "elem-iacs.hpp"

using namespace std;
using namespace utils;
using boost::format;



void Elements::add_elem(string name){
  int ie = -1;
  for (int i=0; i<melemnames.size(); i++){
    if (melemnames[i]==name){ ie=i; break;}
  }
  if (ie==-1){// element not already stored
    cout << "trying to push back " << name << " at end of vector " << melemnames << endl; 
    melemnames.push_back(name);
    int m = melemnames.size();
    // Initialize other properties to defaults:
    matomtype.resize(m);   matomtype[m-1] = m;
    mmass.resize(m);       mmass[m-1]     = 1.0;
    mreflat.resize(m);     mreflat[m-1]   = "none";
    mreflat_a.resize(m);   mreflat_a[m-1] = -1.0;
    mreflat_b.resize(m);   mreflat_b[m-1] = -1.0;
    mreflat_c.resize(m);   mreflat_c[m-1] = -1.0;
    mreflat_bpa.resize(m); mreflat_bpa[m-1] = -1.0;
    mreflat_cpa.resize(m); mreflat_cpa[m-1] = -1.0;
    
  }
  // else element already stored, do nothing
}


int Elements::nelem(void){
  return melemnames.size();
}


int Elements::name2idx(string name){
  for (int i=0; i<melemnames.size(); ++i)
    if (melemnames[i]==name) return i;
  // If we make it here then element name has not been stored. Store it now.
  add_elem(name);
  return nelem()-1;
}

string Elements::idx2name(int i){
  return melemnames[i];
}




int & Elements::atomtype(string s){
  return matomtype[name2idx(s)];
}

double & Elements::mass(string s){
  return mmass[name2idx(s)];
}

string & Elements::reflat(string s){
  return mreflat[name2idx(s)];
}

double & Elements::reflat_a(string s){
  return mreflat_a[name2idx(s)];
}
double & Elements::reflat_b(string s){
  return mreflat_b[name2idx(s)];
}
double & Elements::reflat_c(string s){
  return mreflat_c[name2idx(s)];
}
double & Elements::reflat_bpa(string s){
  return mreflat_bpa[name2idx(s)];
}
double & Elements::reflat_cpa(string s){
  return mreflat_cpa[name2idx(s)];
}





int & Elements::atomtype(int i){
  return matomtype[i];
}

double & Elements::mass(int i){
  return mmass[i];
}

string & Elements::reflat(int i){
  return mreflat[i];
}

double & Elements::reflat_a(int i){
  return mreflat_a[i];
}
double & Elements::reflat_b(int i){
  return mreflat_b[i];
}
double & Elements::reflat_c(int i){
  return mreflat_c[i];
}
double & Elements::reflat_bpa(int i){
  return mreflat_bpa[i];
}
double & Elements::reflat_cpa(int i){
  return mreflat_cpa[i];
}


Vector<double> Elements::masses(){
  return mmass;
}


// ******************************************************************************
// ******************************************************************************

void Interactions::init(Elements el){
  mel = el;
  miacnames.resize(mel.nelem(), mel.nelem());
  for (int i=0; i<mel.nelem(); ++i)
    for (int j=i; j<mel.nelem(); ++j)
      miacnames.elem(i,j) = "none";
}

string & Interactions::name(string n1, string n2){
  int i = mel.name2idx(n1);
  int j = mel.name2idx(n2);
  if (i<=j) return miacnames.elem(i,j);
  else      return miacnames.elem(j,i);
}

string & Interactions::name(int i, int j){
  if (i<=j) return miacnames.elem(i,j);
  else      return miacnames.elem(j,i);
}

