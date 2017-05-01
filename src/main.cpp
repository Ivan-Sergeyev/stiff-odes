#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

#include "rhs.h"
#include "runge_kutta.h"


using std::cerr;
using std::ofstream;
using std::string;
using std::to_string;


int main() {
    long double epsilon[3] = {pow(2, -7), pow(2, -7), pow(2, -7)};
    long double t[3] = {0, 0, 0};
    long double T[3] = {4000, 4000, 2000};
    long double dt[3] = {pow(2, -8), pow(2, -8), pow(2, -8)};
    array<long double, 4> prob_vars_start[3] = {
        {2, 1, 1, 0},
        {2, 1, 1, 2},
        {2, 1, 1, 2}
    };
    long double allow_error = pow(2, -4);
    ofstream out;

    const int RHS_FUNC_NUM = 3;
    RHS_1 rhs_1(epsilon[0]);
    RHS_2 rhs_2(epsilon[1]);
    RHS_3 rhs_3(epsilon[2]);
    const RHS_Function<4> *RHS_FUNC[RHS_FUNC_NUM] = {
        &rhs_1, &rhs_2, &rhs_3
    };

    const int RK_METHOD_NUM = 3;
    RK1_Explicit<4> rk1_expl;
    RK4_Explicit<4> rk4_expl;
    RK4_Implicit<4> rk4_impl(allow_error);
    const RK_Method<4> *RK_METHOD[RK_METHOD_NUM] = {
        &rk1_expl, &rk4_expl, &rk4_impl
    };

    for (int prob_num = 0; prob_num < RHS_FUNC_NUM; ++prob_num) {
        cerr << "==== solving problem #" << prob_num << " ====\n";

        Vector<4> vars(prob_vars_start[prob_num]);

        for (int m = 0; m < RK_METHOD_NUM; ++m) {
            cerr << "==== running method #" << m << " ====\n";

            string data_file_name = \
                "results/data_" + to_string(prob_num + 1) + \
                '_' + to_string(m + 1) + ".txt";

            out.open(data_file_name, std::ios::out);
            simulate_problem<4>(RHS_FUNC[prob_num], RK_METHOD[m],
                                t[prob_num], T[prob_num], dt[prob_num], vars,
                                allow_error, out);
            out.close();
        }
    }

    return 0;
}
