#ifndef LINEAR_ALGEBRA_H_
#define LINEAR_ALGEBRA_H_

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
class Vector {
 public:
    array<double, N> var;

    Vector() {
        var.fill(0);
    }

    Vector(array<double, N> values) {
        copy(values.begin(), values.end(), var.begin());
    }

    Vector(const Vector<N>& other) {
        copy(other.var.begin(), other.var.end(), var.begin());
    }

    double& operator [] (int idx) {
        return var[idx];
    }

    const double operator [] (int idx) const {
        return var[idx];
    }

    Vector<N> operator + (const Vector<N>& other) const {
        Vector<N> v;
        for (int i = 0; i < N; ++i) {
            v[i] = var[i] + other[i];
        }
        return v;
    }

    Vector<N> operator - (const Vector<N>& other) const {
        Vector<N> v;
        for (int i = 0; i < N; ++i) {
            v[i] = var[i] - other[i];
        }
        return v;
    }

    Vector<N> operator * (double n) const {
        Vector<N> v;
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

    double norm_eucl_sq() const {
        double n = 0;
        for (int i = 0; i < N; ++i) {
            n += var[i] * var[i];
        }
        return n;
    }
};

template <int N>
inline Vector<N> operator * (double n, const Vector<N>& v) {
    return v * n;
}

template <int N>
inline Vector<N> operator / (const Vector<N>& v, double n) {
    assert(n);
    return v * (1 / n);
}

template <int N>
ostream& operator << (ostream& out, const Vector<N>& v) {
    for (int i = 0; i < N - 1; ++i) {
        out << v[i] << '\t';
    }
    return out << v[N - 1];
}


template <int N, int M>
class Matrix {
 public:
    array<double, N * M> var;

    Matrix() {
        var.fill(0);
    }

    Matrix(array<double, N * M> values) {
        copy(values.begin(), values.end(), var.begin());
    }

    Matrix(const Matrix<N, M>& other) {
        copy(other.var.begin(), other.var.end(), var.begin());
    }

    double& operator () (double row, double col) {
        return var[row * N + col];
    }

    const double operator () (double row, double col) const {
        return var[row * N + col];
    }

    Matrix<N, M> operator + (const Matrix<N, M>& m) const {
        Matrix<N, M> res;
        for (int r = 0; r < M; ++r) {
            for (int c = 0; c < N; ++c) {
                res(r, c) = (*this)(r, c) + m(r, c);
            }
        }
        return res;
    }

    template <int K>
    Matrix<K, M> operator * (const Matrix<K, N>& m) const {
        Matrix<K, M> res;
        for (int c = 0; c < K; ++c) {
            for (int r = 0; r < M; ++r) {
                for (int t = 0; t < N; ++t) {
                    res(r, c) += (*this)(r, t) * m(t, c);
                }
            }
        }
        return res;
    }

    Matrix<N, M> operator * (double n) const {
        Matrix<N, M> m;
        for (int r = 0; r < M; ++r) {
            for (int c = 0; c < N; ++c) {
                m(r, c) = (*this)(r, c) * n;
            }
        }
        return m;
    }
};

template <int N, int M>
inline Matrix<N, M> operator * (double n, const Matrix<N, M>& m) {
    return m * n;
}

#endif  // LINEAR_ALGEBRA_H_
