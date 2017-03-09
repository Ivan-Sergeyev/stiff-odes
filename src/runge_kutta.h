#ifndef RUNGE_KUTTA_H_
#define RUNGE_KUTTA_H_

#include "variables.h"


template <class VarType>
using rhs_f_homog_with_param = VarType (*) (VarType, double);

template <class VarType>
VarType rk1_iter(double tau, VarType arg, double epsilon,
				 rhs_f_homog_with_param<VarType> rhs) {
	VarType k1 = rhs(arg, epsilon);
	return arg + k1 * tau;
}

template <class VarType>
VarType rk4_iter(double tau, VarType arg, double epsilon,
				 rhs_f_homog_with_param<VarType> rhs) {
	VarType k1 = rhs(arg, epsilon) * tau;
	VarType k2 = rhs(arg + k1 / 2, epsilon) * tau;
	VarType k3 = rhs(arg + k2 / 2, epsilon) * tau;
	VarType k4 = rhs(arg + k3, epsilon) * tau;
	return (arg * 6 + k1 + k2 * 2 + k3 * 2 + k4) / 6;
}

// todo : implement the third method

#endif  // RUNGE_KUTTA_H_
