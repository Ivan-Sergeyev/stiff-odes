#ifndef RHS_H_
#define RHS_H_

#include "linear_algebra.h"


template <int N>
using matr_inv_f = Matrix<N, N> (*)(Matrix<N, N>);


template <int N>
class RHS_Function {
 public:
    matr_inv_f<N> m_inv;

    RHS_Function(matr_inv_f<N> matr_inv) :
        m_inv(matr_inv)
    {}

    virtual Vector<N> f(Vector<N> arg) const = 0;
    virtual Matrix<N, N> j(Vector<N> arg) const = 0;
};


class RHS_1 : public RHS_Function<4> {
 private:
    long double epsilon;

 public:
    RHS_1(long double eps) :
        RHS_Function<4>(matr_inv_3), epsilon(eps)
    {}

    Vector<4> f(Vector<4> arg) const {
        Vector<4> res;
        long double x = arg(0), y = arg(1), a = arg(2);
        res(0) = (1 - x / 2 - 2 * y / 7 / a / a) * x;
        res(1) = (2 * a - 7 * a * a * x / 2 - y / 2) * y;
        res(2) = (2 - 7 * a * x) * epsilon;
        res(3) = 0;  // not used
        return res;
    }

    Matrix<4, 4> j(Vector<4> arg) const {
        Matrix<4, 4> j;
        long double x = arg(0), y = arg(1), a = arg(2);
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
};


class RHS_2 : public RHS_Function<4> {
 private:
    long double epsilon;

 public:
    RHS_2(long double eps) :
        RHS_Function<4>(matr_inv_4), epsilon(eps)
    {}

    Vector<4> f(Vector<4> arg) const {
        Vector<4> res;
        long double x = arg(0), y = arg(1), a1 = arg(2), a2 = arg(3);
        res(0) = (2 * a1 - y * a1 * a1 / a2 / a2 - x / 2) * x;
        res(1) = (2 * a2 - x * a2 * a2 / a1 / a1 - y / 2) * y;
        res(2) = (2 - 2 * y * a1 / a2 / a2) * epsilon;
        res(3) = (2 - 2 * x * a2 / a1 / a1) * epsilon;
        return res;
    }

    Matrix<4, 4> j(Vector<4> arg) const {
        Matrix<4, 4> j;
        long double x = arg(0), y = arg(1), a1 = arg(2), a2 = arg(3);
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
};


class RHS_3 : public RHS_Function<4> {
 private:
    long double epsilon;

 public:
    RHS_3(long double eps) :
        RHS_Function<4>(matr_inv_4), epsilon(eps)
    {}

    Vector<4> f(Vector<4> arg) const {
        Vector<4> res;
        long double x = arg(0), y = arg(1), a1 = arg(2), a2 = arg(3);
        res(0) = (2 * a1 - y * a1 * a1 * a1 / a2 / a2 / a2 - x / 2) * x;
        res(1) = (2 * a2 - x * a2 * a2 * a2 / a1 / a1 / a1 - y / 2) * y;
        res(2) = (2 - 3 * y * a1 * a1 / a2 / a2 / a2) * epsilon;
        res(3) = (2 - 3 * x * a2 * a2 / a1 / a1 / a1) * epsilon;
        return res;
    }

    Matrix<4, 4> j(Vector<4> arg) const {
        Matrix<4, 4> j;
        long double x = arg(0), y = arg(1), a1 = arg(2), a2 = arg(3);
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
};

#endif  // RHS_H_
