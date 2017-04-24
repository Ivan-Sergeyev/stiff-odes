#ifndef RUNGE_KUTTA_H_
#define RUNGE_KUTTA_H_

#include "linear_algebra.h"


template <class VarType>
using rhs_f_homog = VarType (*)(VarType, double);

template <class VarType>
using rhs_newton_solve = VarType (*)(VarType, double);

template <class VarType>
using rk_iter_func = VarType (*)(double, VarType, rhs_f_homog<VarType>, double,
                                 double);

template <class VarType>
VarType rk1_expl_iter(double tau, VarType arg,
                      rhs_f_homog<VarType> rhs, double epsilon,
                      double allow_error) {
    VarType k1(rhs(arg, epsilon));
    return arg + k1 * tau;
}

template <class VarType>
VarType rk4_expl_iter(double tau, VarType arg,
                      rhs_f_homog<VarType> rhs, double epsilon,
                      double allow_error) {
    VarType k1(rhs(arg, epsilon));
    VarType k2(rhs(arg + k1 * tau / 2, epsilon));
    VarType k3(rhs(arg + k2 * tau / 2, epsilon));
    VarType k4(rhs(arg + k3 * tau, epsilon));
    return VarType(arg + (k1 + k2 * 2 + k3 * 2 + k4) * tau / 6);
}

template <class VarType>
VarType rk4_impl_iter(double tau, VarType arg,
                      rhs_f_homog<VarType> rhs,
                      double epsilon,
                      double allow_error) {
    // VarType l_curr = arg;
    VarType l_next = arg;
    // do {
    //     l_curr = l_next;
    //     l_next = rhs_newton_solve(tau, arg, epsilon) + l_curr;
    // } while ((l_curr - l_next).norm_max() > allow_error);
    return l_next;
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
