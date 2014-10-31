

#ifndef PROPFUN_HPP
#define PROPFUN_HPP



#include "utils-vector.hpp"
#include "param.hpp"

using namespace utils;



// Birch-Murnaghan EOS
Vector<double> fun_bmeos(Param & P, Vector<double> & DX);

Vector<double> fun_poly2(Param & P, Vector<double> & DX);
Vector<double> fun_poly2f(Param & P, Vector<double> & DX);

Vector<double> fun_poly3(Param & P, Vector<double> & DX);
Vector<double> fun_poly3f(Param & P, Vector<double> & DX);




#endif


