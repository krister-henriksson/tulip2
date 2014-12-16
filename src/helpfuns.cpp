


#include <limits>
#include <string>

#include <boost/format.hpp>

#include "utils-vector.hpp"
#include "utils-math.hpp"
#include "utils-string.hpp"
#include "utils-errors.hpp"


using namespace utils;
using boost::format;
using std::cout;
using std::endl;


/*
  For an arbitrary x value between x[1] and x[N], end points inclusive,
  get the corresponding y value by doing a spline interpolation
  between the closest data values x[j] and x[j+1].

  For x values outside the interval the function returns an error.

  Codes adapted from Press et al., Numerical Recipes in C.
*/

void spline(Vector<double> & x,
	    Vector<double> & y,
	    Vector<double> & d2y
	    ){

  Vector<double> u(x.size(), 0);
  double s, p, qn, un;
  int i;
  
  d2y[0] = u[0] = 0.0;

  for (i=1; i<x.size()-1; i++){
    s = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
    p = s * d2y[i-1] + 2.0;
    d2y[i] = (s - 1.0) / p;
    u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
    u[i] = ( 6.0 * u[i] / (x[i+1]-x[i-1]) - s * u[i-1] ) / p;
  }

  qn = un = 0.0;
  i = x.size()-1;
  d2y[i] = (un - qn*u[i-1]) / (qn*d2y[i-1] + 1.0);

  for (i=x.size()-2; i>=0; i--){
    d2y[i] = d2y[i] * d2y[i+1] + u[i];
  }
}



double splint(Vector<double> & x,
	      Vector<double> & y,
	      Vector<double> & d2y,
	      double xp
	      ){
  int i, j, N = x.size();
  double A, B, C, D, u, usq, yp;

  if (fp_is_small(xp - x[0]))   return y[0];
  if (fp_is_small(xp - x[N-1])) return y[N-1];

  j = -1;
  for (i=0; i<N-1; i++){
    if (fp_is_small(xp - x[i])) return y[i];
    
    if (xp>=x[i] && xp<=x[i+1]){
      j = i; break;
    }
  }
  if (j <= -1)
    aborterror("ERROR: Supplied point " + tostring(xp) + " is outside the interval [" +
	       tostring(x[0]) + ", " + tostring(x[N-1]) + "]. Exiting.");
  if (j+1 >= N)
    aborterror("ERROR: Supplied point " + tostring(xp) + " is outside the interval [" +
	       tostring(x[0]) + ", " + tostring(x[N-1]) + "]. Exiting.");

  A = (x[j+1] - xp  ) / (x[j+1] - x[j]);
  B = (xp     - x[j]) / (x[j+1] - x[j]);
  u   = (x[j+1] - x[j]);
  usq = u*u;
  C = 1.0/6.0 * A * (A*A - 1.0) * usq;
  D = 1.0/6.0 * B * (B*B - 1.0) * usq;

  yp = A * y[j] + B * y[j+1] + C * d2y[j] + D * d2y[j+1];

  return yp;
}





/* Return derivative at point xp. */

double splint_dy(Vector<double> & x,
		 Vector<double> & y,
		 Vector<double> & d2y,
		 double xp
		 ){
  int i, j, N = x.size();
  double A, B, dyp;

  if (fp_is_small(xp-x[0]))   return 0.0;
  if (fp_is_small(xp-x[N-1])) return 0.0;

  j = -1;
  for (i=0; i<N-1; i++){
    if (xp >= x[i]  && xp < x[i+1]){
      j = i; break;
    }
  }
  if (j <= -1)
    aborterror("ERROR: Supplied point " + tostring(xp) + " is outside the interval [" +
	       tostring(x[0]) + ", " + tostring(x[N-1]) + "]. Exiting.");
  if (j+1 >=N)
    aborterror("ERROR: Supplied point " + tostring(xp) + " is outside the interval [" +
	       tostring(x[0]) + ", " + tostring(x[N-1]) + "]. Exiting.");
  
  A = (x[j+1] - xp  ) / (x[j+1] - x[j]);
  B = (xp     - x[j]) / (x[j+1] - x[j]);
    
  dyp = (y[j+1]-y[j]) / (x[j+1]-x[j])
    - (3.0*A*A-1.0)/6.0 * (x[j+1]-x[j]) * d2y[j]
    + (3.0*B*B-1.0)/6.0 * (x[j+1]-x[j]) * d2y[j+1];

  return dyp;
}










double calc_rel_change(double pred, double readin){
  if (fp_is_small(readin))
    return std::numeric_limits<double>::epsilon();
  else
    return 100.0 * (pred / readin - 1.0);
}


Vector<double> convert_unc_to_wei(Vector<double> & x){
  Vector<double> y(x);
  double small = sqrt( std::numeric_limits<double>::epsilon() );

  for (int i=0; i<x.size(); ++i){
    x[i] = (x[i]<0) ? -x[i] : x[i];
    if (x[i] < small) y[i] = 1.0/small;
    else              y[i] = 1.0/x[i];
  }
  return y;
}



bool get_boolean_choice(std::string ts){
  int nc = ts.size();
  int c1=' ', c2=' ';

  if (nc>=1) c1=ts[0];
  if (nc>=2) c2=ts[1];

  if (c1=='y' || c1=='Y' || c1=='t' || c1=='T')
    return true;
  else if ( (c1=='o' && c2=='n') || (c1=='O' && c2=='N') )
    return true;
  else if ( c1=='1' )
    return true;
  else
    return false;

}



void get_parabolic_fit_from_triplet(Vector<double> x,
				    Vector<double> y,
				    double & a0,
				    double & a2,
				    double & x0,
				    bool debug){
  
  double td1, td2;
  double x1 = x[0], x2 = x[1], x3 = x[2];
  double y1 = y[0], y2 = y[1], y3 = y[2];

  if (debug){
    cout << "x1 x2 x3 : " << x1 << " " << x2 << " " << x3 <<endl;
    cout << "y1 y2 y3 : " << y1 << " " << y2 << " " << y3 <<endl;
  }

  td1 = (y1-y2)*(x1*x1 - x3*x3) - (y1-y3)*(x1*x1 - x2*x2);
  td2 = (y1-y2)*(x1 - x3) - (y1-y3)*(x1 - x2);
  x0 = td1/td2;


  a2 = (y1-y2)/( (x1*x1 - x2*x2) - 2*x0*(x1-x2) );

  a0 = y3 - a2 * (x3-x0)*(x3-x0);

  if (debug){
    cout << "parabolic fit: xmin: " << x0 << endl;
    cout << "parabolic fit: ymin: " << a0 << endl;
    cout << "parabolic fit: curvature: " << a2 << endl;
  }

}

