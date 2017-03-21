#ifndef RUNGE_KUTTA_H_
#define RUNGE_KUTTA_H_

#include "variables.h"


template <class VarType>
using rhs_f_homog = VarType (*)(VarType, double);

template <class VarType>
using rk_iter_func = VarType (*)(double, VarType, rhs_f_homog<VarType>, double);

template <class VarType>
VarType rk1_expl_iter(double tau, VarType arg,
                      rhs_f_homog<VarType> rhs, double epsilon) {
    VarType k1 = rhs(arg, epsilon) * tau;
    return arg + k1;
}

template <class VarType>
VarType rk4_expl_iter(double tau, VarType arg,
                      rhs_f_homog<VarType> rhs, double epsilon) {
    VarType k1 = rhs(arg, epsilon) * tau;
    VarType k2 = rhs(arg + k1 * 0.5, epsilon) * tau;
    VarType k3 = rhs(arg + k2 * 0.5, epsilon) * tau;
    VarType k4 = rhs(arg + k3, epsilon) * tau;
    return arg + (k1 + k2 * 2 + k3 * 2 + k4) / 6;
}

template <class VarType>
VarType rk4_impl_iter(double tau, VarType arg,
                      rhs_f_homog<VarType> rhs, double epsilon) {
// todo : implement the third method
    return arg;
}

const int RK_METHOD_NUM = 3;

template <class VarType>
const rk_iter_func<VarType> RK_METHOD[RK_METHOD_NUM] = {
    rk1_expl_iter, rk4_expl_iter, rk4_impl_iter
};

const int RK_METHOD_APPROX[RK_METHOD_NUM] = {
    1, 4, 4
};

#endif  // RUNGE_KUTTA_H_
