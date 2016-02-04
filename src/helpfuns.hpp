

#ifndef HELPFUNS_HPP
#define HELPFUNS_HPP


#include <string>

#include "utils-vector.hpp"
#include "utils-vector3.hpp"


using utils::Vector;


void spline(utils::Vector<double> & x,
	    utils::Vector<double> & y,
	    utils::Vector<double> & d2y
	    );

double splint(utils::Vector<double> & x,
	      utils::Vector<double> & y,
	      utils::Vector<double> & d2y,
	      double xp
	      );

double splint_dy(utils::Vector<double> & x,
		 utils::Vector<double> & y,
		 utils::Vector<double> & d2y,
		 double xp
		 );

double calc_rel_change(double pred, double readin);

Vector<double> convert_unc_to_wei(Vector<double> & x);


void get_parabolic_fit_from_triplet(utils::Vector3<double> x,
				    utils::Vector3<double> y,
				    double & a0,
				    double & a2,
				    double & x0,
				    bool debug);


#endif
