

#include <cmath>

#include "utils-vector.hpp"
#include "param.hpp"

using namespace utils;



// Birch-Murnaghan EOS
Vector<double> fun_bmeos(Param & P, Vector<double> & DX){
  Vector<double> DY(DX);
  double E0, V0, B0, Bp0;
    
  E0  = P.X(0);
  V0  = P.X(1);
  B0  = P.X(2);
  Bp0 = P.X(3);
    
  for (int i=0; i<DX.size(); ++i)
    DY[i] = E0 + 9.0*V0*B0/16.0 * pow( pow(V0/DX[i], 2.0/3.0) - 1.0, 2 )
      * (
	 ( pow(V0/DX[i], 2.0/3.0) - 1.0 ) * Bp0
	 + 6.0 - 4.0 * pow(V0/DX[i], 2.0/3.0)
	 );
  return DY;
}



Vector<double> fun_poly2(Param & P, Vector<double> & DX){
  Vector<double> DY(DX);
  double t;

  for (int i=0; i<DX.size(); ++i){
    t = DX[i] - P.X(0);
    DY[i] = P.X(1) + P.X(2)*t*t;
  }
  return DY;
}

Vector<double> fun_poly2f(Param & P, Vector<double> & DX){
  Vector<double> DY(DX);
  double t;

  for (int i=0; i<DX.size(); ++i){
    t = DX[i] - P.X(0);
    DY[i] = P.X(1) + P.X(2)*t + P.X(3)*t*t;
  }
  return DY;
}


Vector<double> fun_poly3(Param & P, Vector<double> & DX){
  Vector<double> DY(DX);
  double t;

  for (int i=0; i<DX.size(); ++i){
    t = DX[i] - P.X(0);
    DY[i] = P.X(1) + P.X(2)*t*t + P.X(3)*t*t*t;
  }
  return DY;
}

Vector<double> fun_poly3f(Param & P, Vector<double> & DX){
  Vector<double> DY(DX);
  double t;

  for (int i=0; i<DX.size(); ++i){
    t = DX[i] - P.X(0);
    DY[i] = P.X(1) + P.X(2)*t + P.X(3)*t*t + P.X(4)*t*t*t;
  }
  return DY;
}





