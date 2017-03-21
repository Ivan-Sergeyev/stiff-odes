#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

#include "rhs.h"
#include "runge_kutta.h"


using std::ofstream;
using std::string;
using std::to_string;


template <int N>
void simulate_problem(rhs_func<N> rhs_f, double epsilon,
                      double t, double T, double dt0, Variables<N> vars,
                      rk_iter_func<Variables<N>> rk_iter, int rk_approx,
                      double allow_error, string data_file_name) {
    ofstream out;
    out.open(data_file_name, std::ios::out);

    out << "# parameters:\n"
        << "#   time: from " << t << " to " << T << '\n'
        << "#   starting vars: " << vars << '\n'
        << "#   accuracy: " << allow_error << '\n';

    Variables<N> res1, res2;
    double dt = dt0;
    double step_error;

    out << "# data:\n"
        << t << '\t' << vars << '\n';

    while (t < T) {
        // std::cerr << rk_approx << "  t=" << t << '\n';
        res1 = rk_iter(dt, vars, rhs_f, epsilon);
        res2 = rk_iter(dt / 2, vars, rhs_f, epsilon);
        res2 = rk_iter(dt / 2, res2, rhs_f, epsilon);
        step_error = (res1 - res2).norm_max() / (pow(2, rk_approx) - 1);

        if (step_error > allow_error) {
            dt /= 2;
        } else if (step_error * 100 < allow_error && dt < dt0) {
            dt *= 2;
        } else {
            t += dt;
            vars = res1;
            out << t << '\t' << vars << '\n';
        }
    }

    out.close();
}


template <int N>
void run_methods(rhs_func<N> rhs_f, double epsilon,
                 double t, double T, double dt0, Variables<N> vars,
                 double allow_error, string data_file_name_prefix) {
    for (int m = 0; m < RK_METHOD_NUM; ++m) {
        string data_file_name(data_file_name_prefix +
                              '_' + to_string(m + 1) + ".txt");
        rk_iter_func<Variables<N>> rk_iter_f = RK_METHOD<Variables<N>>[m];
        int approx = RK_METHOD_APPROX[m];
        simulate_problem<N>(rhs_f, epsilon, t, T, dt0, vars,
                            rk_iter_f, approx, allow_error,
                            data_file_name);
    }
}

int main() {
    double prob_epsilon[3] = {pow(2, -7), pow(2, -7), pow(2, -7)};
    double prob_t_start[3] = {0, 0, 0};
    double prob_t_end[3] = {50, 50, 50};
    double prob_dt_start[3] = {pow(2, -10), pow(2, -10), pow(2, -10)};
    array<double, 4> prob_vars_start[3] = {
        {2, 1, 1, 0},
        {2, 1, 1, 5},
        {2, 1, 1, 5}
    };

    for (int prob_num = 0; prob_num < 3; ++prob_num) {
        double epsilon = prob_epsilon[prob_num];
        double t = prob_t_start[prob_num];
        double T = prob_t_end[prob_num];
        double dt = prob_dt_start[prob_num];
        double allow_error = pow(2, -4);
        string data_file_name_prefix = \
            "results/data_" + to_string(prob_num + 1);

        rhs_func<4> rhs_f = RHS_FUNC[prob_num];
        Variables<4> vars(prob_vars_start[prob_num]);
        run_methods<4>(rhs_f, epsilon, t, T, dt, vars, allow_error,
                       data_file_name_prefix);
    }

    return 0;
}
