#ifndef RHS_H_
#define RHS_H_

#include "variables.h"

template <int N>
using rhs_func = Variables<N> (*)(Variables<N>, double);

Variables<3> rhs_1(Variables<3> arg, double epsilon) {
    Variables<3> res;
    double x = arg[0], y = arg[1], a = arg[2];
    res[0] = (1 - x / 2 - 2 * y / 7 / a / a) * x;
    res[1] = (2 * a - 7 * a * a * x / 2 - y / 2) * y;
    res[2] = (2 - 7 * a * x) * epsilon;
    return res;
}

Variables<4> rhs_2(Variables<4> arg, double epsilon) {
    Variables<4> res;
    double x = arg[0], y = arg[1], a1 = arg[2], a2 = arg[3];
    res[0] = (2 * a1 - x / 2 - a1 * a1 * y / a2 / a2) * x;
    res[1] = (2 * a2 - x * a2 * a2 / a1 / a1 - y / 2) * y;
    res[2] = (2 - 2 * a1 * y / a2 / a2) * epsilon;
    res[3] = (2 - 2 * a2 * x / a1 / a1) * epsilon;
    return res;
}

Variables<4> rhs_3(Variables<4> arg, double epsilon) {
    Variables<4> res;
    double x = arg[0], y = arg[1], a1 = arg[2], a2 = arg[3];
    res[0] = (2 * a1 - x / 2 - a1 * a1 * a1 * y / a2 / a2 / a2) * x;
    res[1] = (2 * a2 - x * a2 * a2 * a2 / a1 / a1 / a1 - y / 2) * y;
    res[2] = (2 - 3 * a1 * a1 * y / a2 / a2 / a2) * epsilon;
    res[3] = (2 - 3 * a2 * a2 * x / a1 / a1 / a1) * epsilon;
    return res;
}

#endif  // RHS_H_
