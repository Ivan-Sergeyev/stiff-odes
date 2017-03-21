#ifndef VARIABLES_H_
#define VARIABLES_H_

#include <algorithm>
#include <array>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <string>


using std::array;
using std::copy;
using std::fill_n;
using std::ostream;
using std::string;
using std::to_string;


template <int N>
class Variables {
 public:
    array<double, N> var;

    Variables() {
        var.fill(0);
    }

    Variables(array<double, N> values) {
        copy(values.begin(), values.end(), var.begin());
    }

    Variables(const Variables<N>& other) {
        copy(other.var.begin(), other.var.end(), var.begin());
    }

    double& operator [] (int idx) {
        return var[idx];
    }

    const double operator [] (int idx) const {
        return var[idx];
    }

    Variables<N> operator + (const Variables<N>& other) const {
        Variables v;
        for (int i = 0; i < N; ++i) {
            v[i] = var[i] + other[i];
        }
        return v;
    }

    Variables<N> operator - (const Variables<N>& other) const {
        Variables v;
        for (int i = 0; i < N; ++i) {
            v[i] = var[i] - other[i];
        }
        return v;
    }

    Variables<N> operator * (double n) const {
        Variables v;
        for (int i = 0; i < N; ++i) {
            v[i] = var[i] * n;
        }
        return v;
    }

    double norm_max() const {
        double n = var[0];
        for (int i = 1; i < N; ++i) {
            if (var[i] > n) {
                n = var[i];
            }
        }
        return n;
    }

    double norm_abs() const {
        double n = 0;
        for (int i = 0; i < N; ++i) {
            n += abs(var[i]);
        }
        return n;
    }

    double norm_eucl() const {
        double n = 0;
        for (int i = 0; i < N; ++i) {
            n += var[i] * var[i];
        }
        return sqrt(n);
    }
};

template <int N>
inline Variables<N> operator * (double n, const Variables<N>& v) {
    return v * n;
}

template <int N>
inline Variables<N> operator / (const Variables<N>& v, double n) {
    assert(n);
    return v * (1 / n);
}

template <int N>
ostream& operator << (ostream& out, const Variables<N> v) {
    for (int i = 0; i < N - 1; ++i) {
        out << v[i] << '\t';
    }
    return out << v[N - 1];
}

#endif  // VARIABLES_H_
