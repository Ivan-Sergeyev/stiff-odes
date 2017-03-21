#ifndef RHS_H_
#define RHS_H_

#include "variables.h"

template <int N>
using rhs_func = Variables<N> (*)(Variables<N>, double);

Variables<4> rhs_1(Variables<4> arg, double epsilon) {
    Variables<4> res;
    double x = arg[0], y = arg[1], a = arg[2];
    res[0] = (1 - x / 2 - 2 * y / 7 / a / a) * x;
    res[1] = (2 * a - 7 * a * a * x / 2 - y / 2) * y;
    res[2] = (2 - 7 * a * x) * epsilon;
    res[3] = 0;  // not used
    return res;
}

Variables<4> rhs_2(Variables<4> arg, double epsilon) {
    Variables<4> res;
    double x = arg[0], y = arg[1], a1 = arg[2], a2 = arg[3];
    res[0] = (2 * a1 - y * a1 * a1 / a2 / a2 - x / 2) * x;
    res[1] = (2 * a2 - x * a2 * a2 / a1 / a1 - y / 2) * y;
    res[2] = (2 - 2 * y * a1 / a2 / a2) * epsilon;
    res[3] = (2 - 2 * x * a2 / a1 / a1) * epsilon;
    return res;
}

Variables<4> rhs_3(Variables<4> arg, double epsilon) {
    Variables<4> res;
    double x = arg[0], y = arg[1], a1 = arg[2], a2 = arg[3];
    res[0] = (2 * a1 - y * a1 * a1 * a1 / a2 / a2 / a2 - x / 2) * x;
    res[1] = (2 * a2 - x * a2 * a2 * a2 / a1 / a1 / a1 - y / 2) * y;
    res[2] = (2 - 3 * y * a1 * a1 / a2 / a2 / a2) * epsilon;
    res[3] = (2 - 3 * x * a2 * a2 / a1 / a1 / a1) * epsilon;
    return res;
}

const int RHS_FUNC_NUM = 3;

const rhs_func<4> RHS_FUNC[RHS_FUNC_NUM] = {rhs_1, rhs_2, rhs_3};

#endif  // RHS_H_
