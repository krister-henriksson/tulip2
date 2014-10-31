
#ifndef LATCALC_HPP
#define LATCALC_HPP




#include "utils.hpp"
#include "utils-vector.hpp"

#include "compound.hpp"
#include "param-pot.hpp"


using namespace std;
using namespace utils;


Vector<double> latcalc(ParamPot & param, Vector<CompoundStructureFit> & DX);



#endif

