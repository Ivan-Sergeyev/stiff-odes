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


template <int N>
void simulate_problem(rhs_func<N> rhs_f, double epsilon,
                      double t, double T, double dt0, double dt_min,
                      Vector<N> vars,
                      rk_iter_func<Vector<N>> rk_iter, int rk_approx,
                      double allow_error, string data_file_name) {
    ofstream out;
    out.open(data_file_name, std::ios::out);

    out << "# parameters:\n"
        << "#   time: from " << t << " to " << T << '\n'
        << "#   starting vars: " << vars << '\n'
        << "#   accuracy: " << allow_error << '\n';

    Vector<N> res1, res2;
    double dt = dt0;
    double step_error;

    const int STEPS_UNTIL_FLUSH = 1000;
    int steps = 0;

    out << "# data:\n"
        << t << '\t' << vars << '\n';

    while (t < T) {
        res1 = rk_iter(dt, vars, rhs_f, epsilon, allow_error);
        res2 = rk_iter(dt / 2, vars, rhs_f, epsilon, allow_error);
        res2 = rk_iter(dt / 2, res2, rhs_f, epsilon, allow_error);
        step_error = (res1 - res2).norm_max() / (pow(2, rk_approx) - 1);

        if (dt / 2 >= dt_min && step_error > allow_error) {
            dt /= 2;
        } else if (step_error * 64 < allow_error && dt < dt0) {
            dt *= 2;
        } else {
            t += dt;
            vars = res1;
            out << t << '\t' << vars << '\n';

            ++steps;
            if (!(steps % STEPS_UNTIL_FLUSH)) {
                cerr << t << ' ' << dt << '\n';
                steps = 0;
                out.flush();
            }
        }
    }

    out.close();
}


template <int N>
void run_methods(rhs_func<N> rhs_f, double epsilon,
                 double t, double T, double dt0, Vector<N> vars,
                 double allow_error, string data_file_name_prefix) {
    for (int m = 0; m < 2; ++m) {
        cerr << "==== running method #" << m << " ====\n";
        string data_file_name(data_file_name_prefix +
                              '_' + to_string(m + 1) + ".txt");
        rk_iter_func<Vector<N>> rk_iter_f = RK_METHOD<Vector<N>>[m];
        int approx = RK_METHOD_APPROX[m];
        simulate_problem<N>(rhs_f, epsilon, t, T, dt0, dt0 / 256, vars,
                            rk_iter_f, approx, allow_error,
                            data_file_name);
    }
}

int main() {
    // todo : use actual initial conditions for the problem
    // todo : read initial conditions from file
    double prob_epsilon[3] = {pow(2, -7), pow(2, -7), pow(2, -7)};
    double prob_t_start[3] = {0, 0, 0};
    double prob_t_end[3] = {4000, 4000, 2000};
    double prob_dt_start[3] = {pow(2, -8), pow(2, -8), pow(2, -8)};
    array<double, 4> prob_vars_start[3] = {
        {2, 1, 1, 0},
        {2, 1, 1, 3},
        {2, 1, 1, 3}
    };
    double allow_error = pow(2, -4);

    for (int prob_num = 0; prob_num < 3; ++prob_num) {
        cerr << "==== solving problem #" << prob_num << " ====\n";
        double epsilon = prob_epsilon[prob_num];
        double t = prob_t_start[prob_num];
        double T = prob_t_end[prob_num];
        double dt = prob_dt_start[prob_num];
        string data_file_name_prefix = \
            "results/data_" + to_string(prob_num + 1);

        rhs_func<4> rhs_f = RHS_FUNC[prob_num];
        Vector<4> vars(prob_vars_start[prob_num]);
        run_methods<4>(rhs_f, epsilon, t, T, dt, vars, allow_error,
                       data_file_name_prefix);
    }

    return 0;
}
