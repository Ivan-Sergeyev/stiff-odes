#ifndef RUNGE_KUTTA_H_
#define RUNGE_KUTTA_H_

#include <fstream>

#include "linear_algebra.h"
#include "rhs.h"


using std::ofstream;


template <int N>
class RK_Method {
 public:
    int order;

    RK_Method(int method_order) :
        order(method_order)
    {}

    virtual Vector<N> iter(const RHS_Function<N> *rhs,
                           long double tau, Vector<N> arg) const = 0;
};


template <int N>
class RK1_Explicit : public RK_Method<N> {
 public:
    RK1_Explicit() :
        RK_Method<N>(1)
    {}

    Vector<N> iter(const RHS_Function<N> *rhs,
                   long double tau, Vector<N> arg) const {
        Vector<N> k1 = rhs->f(arg);
        return arg + k1 * tau;
    }
};


template <int N>
class RK4_Explicit : public RK_Method<N> {
 public:
    RK4_Explicit() :
        RK_Method<N>(1)
    {}

    Vector<N> iter(const RHS_Function<N> *rhs,
                   long double tau, Vector<N> arg) const {
        Vector<N> k1 = rhs->f(arg);
        Vector<N> k2 = rhs->f(arg + k1 * tau / 2);
        Vector<N> k3 = rhs->f(arg + k2 * tau / 2);
        Vector<N> k4 = rhs->f(arg + k3 * tau);
        return arg + (k1 + k2 * 2 + k3 * 2 + k4) * tau / 6;
    }
};


template <int N>
class RK4_Implicit : public RK_Method<N> {
 private:
    long double allow_error;

    Vector<N> newton_solve(const RHS_Function<N> *rhs, long double tau,
                           Vector<N> u, Vector<N> l) const {
        long double c1 = (3 + 2 * sqrt(3)) / 12;
        long double c2 = (2 * sqrt(3) - 3) * tau / 3;

        Vector<N> arg1 = l * tau + u;
        Vector<N> f1 = rhs->f(arg1);
        Matrix<N, N> jac1 = rhs->j(arg1);

        Vector<N> arg2 = (l * 3 - f1) * c2 + u;
        Vector<N> f2 = rhs->f(arg2);
        Matrix<N, N> jac2 = rhs->j(arg2);

        Matrix<N, N> e = unit_matrix<N>();
        Matrix<N, N> h = jac2 * (e * 3 - jac1 * tau) * tau / 12 + \
                         jac1 * tau / 4 - e;
        Matrix<N, N> inv = rhs->m_inv(h);

        Vector<N> g = f2 * c1 - l + f1 / 4;
        return l - inv * g;
    }

 public:
    RK4_Implicit(long double error) :
        RK_Method<N>(1), allow_error(error)
    {}

    Vector<N> iter(const RHS_Function<N> *rhs,
                   long double tau, Vector<N> arg) const {
        Vector<N> l_curr;
        Vector<N> l_next = Vector<N>();
        for (int i = 0; i < 3; ++i) {
            l_curr = l_next;
            l_next = newton_solve(rhs, tau, arg, l_curr);
        }

        Vector<N> l1 = l_next;
        Vector<N> k2 = rhs->f(arg + l1 * tau);
        Vector<N> k1 = (l1 * 4 - k2) * (2 * sqrt(3) - 3);
        return arg + (k1 + k2) * tau / 2;
    }
};


template <int N>
void simulate_problem(const RHS_Function<N> *rhs, const RK_Method<N> *rkm,
                      long double t, long double T, long double dt0,
                      Vector<N> vars,
                      long double allow_error, ofstream &out) {
    out << "# parameters:\n"
        << "#   time: from " << t << " to " << T << '\n'
        << "#   starting vars: " << vars << '\n'
        << "#   accuracy: " << allow_error << '\n';

    Vector<N> res1, res2;
    long double dt = dt0;
    long double step_error;
    bool step_inc = false;
    bool step_dec = false;
    const long double error_denom = pow(2, rkm->order) - 1;

    const int STEPS_UNTIL_FLUSH = 1000;
    int steps = 0;

    out << "# data:\n"
        << t << '\t' << vars << '\n';

    while (t < T) {
        res1 = rkm->iter(rhs, dt,     vars);
        res2 = rkm->iter(rhs, dt / 2, vars);
        res2 = rkm->iter(rhs, dt / 2, res2);
        step_error = (res1 - res2).norm_max() / error_denom / res1.norm_max();

        if (step_error >= allow_error && dt / 4 > 0 &&
            (step_dec || (!step_dec && !step_inc))) {
            dt /= 2;
            step_dec = true;
        } else if (step_error < allow_error && dt < dt0 &&
                   (step_inc || (!step_dec && !step_inc))) {
            dt *= 2;
            step_inc = true;
        } else {
            step_inc = step_dec = false;
            t += dt;
            vars = res1;
            out << t << '\t' << vars << '\n';

            ++steps;
            if (!(steps % STEPS_UNTIL_FLUSH)) {
                cerr << "t = " << t << "\tdt = " << dt << '\n';
                steps = 0;
                out.flush();
            }
        }
    }
}

#endif  // RUNGE_KUTTA_H_
