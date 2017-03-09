#include <iomanip>
#include <fstream>
#include <string>

#include "runge_kutta.h"
#include "variables.h"


using std::ofstream;
using std::string;
using std::to_string;


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


template <int N>
void simulate_problem(double t, double T, double dt, double epsilon,
					  rhs_func<N> rhs_f, Variables<N> vars, string data_file_name) {
	ofstream out;
	out.open(data_file_name, std::ios::out);

	// todo : write parameters and initial conditions
	//        as comments in the beginning of the data file

	while (t < T) {
		out << t << '\t' << vars << '\n';
		vars = rk1_iter(dt, vars, epsilon, rhs_f);
		// todo : add other methods as function pointer parameter
		t += dt;
	}
	out << t << '\t' << vars << '\n';

	out.close();
}


int main() {
	double prob_t_start[3] = {0, 0, 0};
	double prob_t_end[3] = {20, 20, 20};
	double prob_dt_start[3] = {0.0625, 0.0625, 0.0625};
	double prob_epsilon[3] = {0.0078125, 0.0078125, 0.0009765625};
	rhs_func<4> prob_rhs_func[3] = {nullptr, rhs_2, rhs_3};
	array<double, 4> prob_vars_start[3] = {
		{2, 1, 1, 0}, {2, 1, 1, 1}, {2, 1, 1, 1}
	};

	for (int prob_num = 0; prob_num < 3; ++prob_num) {
		double t = prob_t_start[prob_num];
		double T = prob_t_end[prob_num];
		double dt = prob_dt_start[prob_num];
		double epsilon = prob_epsilon[prob_num];
		string data_file_name = \
			"results/data_" + to_string(prob_num + 1) + ".txt";
		// todo : add method identifier to data file name

		if (prob_num == 0) {
			rhs_func<3> rhs_f = rhs_1;
			array<double, 3> vars_start = {
				prob_vars_start[0][0],
				prob_vars_start[0][1],
				prob_vars_start[0][2]
			};
			Variables<3> vars(vars_start);

			simulate_problem<3> (t, T, dt, epsilon, rhs_f, vars, data_file_name);
			// todo : run other methods
		} else {
			rhs_func<4> rhs_f = prob_rhs_func[prob_num];
			Variables<4> vars(prob_vars_start[prob_num]);

			simulate_problem<4> (t, T, dt, epsilon, rhs_f, vars, data_file_name);
			// todo : run other methods
		}
	}

	return 0;
}
