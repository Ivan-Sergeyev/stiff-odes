#ifndef RHS_H_
#define RHS_H_

#include "linear_algebra.h"


Vector<4> rhs_1(Vector<4> arg, double epsilon) {
    Vector<4> res;
    double x = arg[0], y = arg[1], a = arg[2];
    res[0] = (1 - x / 2 - 2 * y / 7 / a / a) * x;
    res[1] = (2 * a - 7 * a * a * x / 2 - y / 2) * y;
    res[2] = (2 - 7 * a * x) * epsilon;
    res[3] = 0;  // not used
    return res;
}

Vector<4> rhs_2(Vector<4> arg, double epsilon) {
    Vector<4> res;
    double x = arg[0], y = arg[1], a1 = arg[2], a2 = arg[3];
    res[0] = (2 * a1 - y * a1 * a1 / a2 / a2 - x / 2) * x;
    res[1] = (2 * a2 - x * a2 * a2 / a1 / a1 - y / 2) * y;
    res[2] = (2 - 2 * y * a1 / a2 / a2) * epsilon;
    res[3] = (2 - 2 * x * a2 / a1 / a1) * epsilon;
    return res;
}

Vector<4> rhs_3(Vector<4> arg, double epsilon) {
    Vector<4> res;
    double x = arg[0], y = arg[1], a1 = arg[2], a2 = arg[3];
    res[0] = (2 * a1 - y * a1 * a1 * a1 / a2 / a2 / a2 - x / 2) * x;
    res[1] = (2 * a2 - x * a2 * a2 * a2 / a1 / a1 / a1 - y / 2) * y;
    res[2] = (2 - 3 * y * a1 * a1 / a2 / a2 / a2) * epsilon;
    res[3] = (2 - 3 * x * a2 * a2 / a1 / a1 / a1) * epsilon;
    return res;
}


Matrix<4, 4> rhs_1_jacobi(Vector<4> arg, double epsilon) {
    Matrix<4, 4> j;
    double x = arg[0], y = arg[1], a = arg[2];
    j(0, 0) = 1 - x - 2 * y / 7 / a / a;
    j(0, 1) = -2 * x / 7 / a / a;
    j(0, 2) = 4 * x * y / 7 / a / a / a;
    j(1, 0) = -7 * a * a * y / 2;
    j(1, 1) = 2 * a - 7 * a * a * x / 2 - y;
    j(1, 2) = 2 * y - 7 * a * x * y;
    j(2, 0) = -7 * epsilon * a;
    j(2, 1) = 0;
    j(2, 2) = -7 * epsilon * x;
    j(0, 3) = j(1, 3) = j(2, 3) = 0;            // not used
    j(3, 0) = j(3, 1) = j(3, 2) = j(3, 3) = 0;  // not used
    return j;
}

Matrix<4, 4> rhs_2_jacobi(Vector<4> arg, double epsilon) {
    Matrix<4, 4> j;
    double x = arg[0], y = arg[1], a1 = arg[2], a2 = arg[3];
    j(0, 0) = 2 * a1 - x - a1 * a1 * y / a2 / a2;
    j(0, 1) = -a1 * a1 * x / a2 / a2;
    j(0, 2) = 2 * x - 2 * a1 * x * y / a2 / a2;
    j(0, 3) = 2 * a1 * a1 * x * y / a2 / a2 / a2;
    j(1, 0) = -a2 * a2 * y / a1 / a1;
    j(1, 1) = 2 * a2 - a2 * a2 * x / a1 / a1 - y;
    j(1, 2) = 2 * a2 * a2 * x * y / a1 / a1 / a1;
    j(1, 3) = 2 * y - 2 * a2 * x * y / a1 / a1;
    j(2, 0) = 0;
    j(2, 1) = -2 * epsilon * a1 / a2 / a2;
    j(2, 2) = -2 * epsilon * y / a2 / a2;
    j(2, 3) = 4 * epsilon * a1 * y / a2 / a2 / a2;
    j(3, 0) = -2 * epsilon * a2 / a1 / a1;
    j(3, 1) = 0;
    j(3, 2) = 4 * epsilon * a2 * x / a1 / a1 / a1;
    j(3, 3) = -2 * epsilon * x / a1 / a1;
    return j;
}

Matrix<4, 4> rhs_3_jacobi(Vector<4> arg, double epsilon) {
    Matrix<4, 4> j;
    double x = arg[0], y = arg[1], a1 = arg[2], a2 = arg[3];
    j(0, 0) = 2 * a1 - x - a1 * a1 * a1 * y / a2 / a2 / a2;
    j(0, 1) = -a1 * a1 * a1 * x / a2 / a2 / a2;
    j(0, 2) = 2 * x - 3 * a1 * a1 * x * y / a2 / a2 / a2;
    j(0, 3) = 3 * a1 * a1 * a1 * x * y / a2 / a2 / a2 / a2;
    j(1, 0) = -a2 * a2 * a2 * y / a1 / a1 / a1;
    j(1, 1) = 2 * a2 - a2 * a2 * a2 * x / a1 / a1 / a1 - y;
    j(1, 2) = 3 * a2 * a2 * a2 * x * y / a1 / a1 / a1 / a1;
    j(1, 3) = 2 * y - 3 * a2 * a2 * x * y / a1 / a1 / a1;
    j(2, 0) = 0;
    j(2, 1) = -2 * epsilon * a1 * a1 / a2 / a2 / a2;
    j(2, 2) = -4 * epsilon * a1 * y / a2 / a2 / a2;
    j(2, 3) = 6 * epsilon * a1 * a1 * y / a2 / a2 / a2 / a2;
    j(3, 0) = -2 * epsilon * a2 * a2 / a1 / a1 / a1;
    j(3, 1) = 0;
    j(3, 2) = 6 * epsilon * a2 * a2 * x / a1 / a1 / a1 / a1;
    j(3, 3) = -4 * epsilon * a2 * x / a1 / a1 / a1;
    return j;
}

Vector<4> rhs_1_newton_solve(double tau, Vector<4> u, Vector<4> l, double epsilon) {
    double c1 = (3 + 2 * sqrt(3)) / 12;
    double c2 = (2 * sqrt(3) - 3) * tau / 3;
    double c3 = 0.25;
    Vector<4> arg_1 = (3 * l - rhs_1(u + l, epsilon)) * c2 + u;
    Vector<4> arg_2 = u + l;
    Matrix<4, 4> jac1 = rhs_1_jacobi(arg_1, epsilon);
    Matrix<4, 4> jac2 = rhs_1_jacobi(arg_2, epsilon);
    Matrix<4, 4> j = c1 * (-c2) * jac1 * jac2 + c3 * jac2;
    Matrix<4, 4> inv;
    double det = j(0, 0) * j(1, 1) * j(2, 2) + j(0, 1) * j(1, 2) * j(2, 0) - \
               - j(0, 1) * j(1, 0) * j(2, 2) - j(0, 2) * j(1, 1) * j(2, 0);
    inv(0, 0) = j(1, 1) * j(2, 2) / det;
    inv(0, 1) = -j(0, 1) * j(2, 2) / det;
    inv(0, 2) = (j(0, 1) * j(1, 2) - j(0, 2) * j(1, 1)) / det;
    inv(1, 0) = (j(1, 2) * j(2, 0) - j(1, 0) * j(2, 2)) / det;
    inv(1, 1) = (j(0, 0) * j(2, 2) - j(0, 2) * j(2, 0)) / det;
    inv(1, 2) = (j(0, 2) * j(1, 0) - j(0, 0) * j(1, 2)) / det;
    inv(2, 0) = -j(1, 1) * j(2, 0) / det;
    inv(2, 1) = j(0, 1) * j(2, 0) / det;
    inv(2, 2) = (j(0, 0) * j(1, 1) - j(0, 1) * j(1, 0)) / det;

    Vector<4> res;
    return res;
}

template <int N>
using rhs_func = Vector<N> (*)(Vector<N>, double);

template <int N>
using rhs_jacobi = Matrix<N, N> (*)(Vector<N>, double);

template <int N>
using rhs_newton = Vector<N> (*)(Vector<N>, double);

const int RHS_FUNC_NUM = 3;

const rhs_func<4> RHS_FUNC[RHS_FUNC_NUM] = {
    rhs_1, rhs_2, rhs_3
};

const rhs_jacobi<4> RHS_JACOBI[RHS_FUNC_NUM] = {
    rhs_1_jacobi, rhs_2_jacobi, rhs_3_jacobi
};

#endif  // RHS_H_
