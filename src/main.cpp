#include <iomanip>
#include <fstream>
#include <string>

#include "rhs.h"
#include "runge_kutta.h"


using std::ofstream;
using std::string;
using std::to_string;


template <int N>
void simulate_problem(rhs_func<N> rhs_f, double epsilon,
                      double t, double T, double dt0, Variables<N> vars,
                      rk_iter_func<Variables<N>> rk_iter, double accuracy,
                      string data_file_name) {
    ofstream out;
    out.open(data_file_name, std::ios::out);

    out << "# parameters:\n"
        << "#   epsilon = " << epsilon << '\n'
        << "# initial values:\n"
        << "#   t = " << t << '\n'
        << "#   T = " << T << '\n'
        << "#   vars = " << vars << '\n'
        << "#   accuracy = " << accuracy << '\n';

    Variables<N> res1, res2;
    double dt = dt0;

    out << "# data:\n"
        << t << '\t' << vars << '\n';
    while (t < T) {
        res1 = rk_iter(dt, vars, epsilon, rhs_f);
        res2 = rk_iter(dt / 2, vars, epsilon, rhs_f);
        res2 = rk_iter(dt / 2, res2, epsilon, rhs_f);
        double step_acc = (res1 - res2).norm_eucl();
        // todo : add denominator 2^(p-1)-1 (see Petrov, p. 200)

        if (step_acc > accuracy) {
            dt /= 2;
        } else if (step_acc * 100 < accuracy && dt < dt0) {
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
                 double accuracy, string data_file_name_prefix) {
    for (int m = 0; m < NUM_METHODS; ++m) {
        string data_file_name(data_file_name_prefix +
                              '_' + to_string(m + 1) + ".txt");
        simulate_problem<N>(rhs_f, epsilon, t, T, dt0, vars,
                            RK_METHOD<Variables<N>>[m], accuracy,
                            data_file_name);
    }
}

int main() {
    rhs_func<4> prob_rhs_func[3] = {nullptr, rhs_2, rhs_3};
    double prob_epsilon[3] = {0.0078125, 0.0078125, 0.0009765625};

    double prob_t_start[3] = {0, 0, 0};
    double prob_t_end[3] = {20, 20, 20};
    double prob_dt_start[3] = {0.0625, 0.0625, 0.0625};
    array<double, 4> prob_vars_start[3] = {
        {2, 1, 1, 0}, {2, 1, 1, 1}, {2, 1, 1, 1}
    };

    for (int prob_num = 0; prob_num < 3; ++prob_num) {
        double t = prob_t_start[prob_num];
        double T = prob_t_end[prob_num];
        double dt = prob_dt_start[prob_num];
        double epsilon = prob_epsilon[prob_num];
        double accuracy = 0.0078125;
        string data_file_name_prefix = \
            "results/data_" + to_string(prob_num + 1);

        if (prob_num == 0) {
            rhs_func<3> rhs_f = rhs_1;
            array<double, 3> vars_start = {
                prob_vars_start[0][0],
                prob_vars_start[0][1],
                prob_vars_start[0][2]
            };
            Variables<3> vars(vars_start);

            run_methods<3>(rhs_f, epsilon, t, T, dt, vars, accuracy,
                           data_file_name_prefix);
        } else {
            rhs_func<4> rhs_f = prob_rhs_func[prob_num];
            Variables<4> vars(prob_vars_start[prob_num]);

            run_methods<4>(rhs_f, epsilon, t, T, dt, vars, accuracy,
                           data_file_name_prefix);
        }
    }

    return 0;
}
