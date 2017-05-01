#ifndef RUNGE_KUTTA_TEST_H_
#define RUNGE_KUTTA_TEST_H_

#include "rhs.h"
#include "runge_kutta.h"


class RHS_Test_Function : public RHS_Function<4> {
 public:
    RHS_Test_Function() :
        RHS_Function<4>(matr_inv_4)
    {}

    Vector<4> f(Vector<4> arg) const {
        return arg;
    }

    Matrix<4, 4> j(Vector<4> arg) const {
        return unit_matrix<4>();
    }
};


void runge_kutta_test() {
    long double allow_error = 0.00390625;

    const int RK_METHOD_NUM = 3;
    RK1_Explicit<4> rk1_expl;
    RK4_Explicit<4> rk4_expl;
    RK4_Implicit<4> rk4_impl(allow_error);
    const RK_Method<4> *RK_METHOD[RK_METHOD_NUM] = {
        &rk1_expl, &rk4_expl, &rk4_impl
    };

    RHS_Test_Function rtf;
    ofstream out;

    for (int m = 0; m < RK_METHOD_NUM; ++m) {
        cerr << "==== running method #" << m << " ====\n";

        Vector<4> vars({0.25, 0.5, 1, 2});
        string data_file_name = \
            "results_test/data_" + to_string(m + 1) + ".txt";

        out.open(data_file_name, std::ios::out);
        simulate_problem<4>(&rtf, RK_METHOD[m],
                            0, 10, 0.00390625, vars,
                            allow_error, out);
        out.close();
    }
}

#endif  // RUNGE_KUTTA_TEST_H_
