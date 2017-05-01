#ifndef LINEAR_ALGEBRA_H_
#define LINEAR_ALGEBRA_H_

#include <algorithm>
#include <array>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <string>


using std::abs;
using std::array;
using std::cerr;
using std::copy;
using std::fill_n;
using std::ostream;
using std::string;
using std::to_string;


template <int N>
class Vector {
 public:
    array<long double, N> var;

    Vector() {
        var.fill(0);
    }

    Vector(array<long double, N> values) {
        copy(values.begin(), values.end(), var.begin());
    }

    Vector(const Vector<N>& other) {
        copy(other.var.begin(), other.var.end(), var.begin());
    }

    long double& operator () (int idx) {
        return var[idx];
    }

    const long double operator () (int idx) const {
        return var[idx];
    }

    Vector<N> operator + (const Vector<N>& other) const {
        Vector<N> v;
        for (int i = 0; i < N; ++i) {
            v(i) = var[i] + other(i);
        }
        return v;
    }

    Vector<N> operator - (const Vector<N>& other) const {
        Vector<N> v;
        for (int i = 0; i < N; ++i) {
            v(i) = var[i] - other(i);
        }
        return v;
    }

    Vector<N> operator * (long double n) const {
        Vector<N> v;
        for (int i = 0; i < N; ++i) {
            v(i) = var[i] * n;
        }
        return v;
    }

    Vector<N> operator / (long double n) const {
        assert(n);
        Vector<N> v;
        for (int i = 0; i < N; ++i) {
            v(i) = var[i] / n;
        }
        return v;
    }

    long double norm_max() const {
        long double n = var[0];
        for (int i = 1; i < N; ++i) {
            if (abs(var[i]) > n) {
                n = var[i];
            }
        }
        return n;
    }

    long double norm_abs() const {
        long double n = 0;
        for (int i = 0; i < N; ++i) {
            n += abs(var[i]);
        }
        return n;
    }

    long double norm_eucl_sq() const {
        long double n = 0;
        for (int i = 0; i < N; ++i) {
            n += var[i] * var[i];
        }
        return n;
    }
};

template <int N>
inline Vector<N> operator * (long double n, const Vector<N>& v) {
    return v * n;
}

template <int N>
ostream& operator << (ostream& out, const Vector<N>& v) {
    for (int i = 0; i < N - 1; ++i) {
        out << v(i) << '\t';
    }
    return out << v(N - 1);
}


template <int N, int M>
class Matrix {
 public:
    array<long double, N * M> var;

    Matrix() {
        var.fill(0);
    }

    Matrix(array<long double, N * M> values) {
        copy(values.begin(), values.end(), var.begin());
    }

    Matrix(const Matrix<N, M>& other) {
        copy(other.var.begin(), other.var.end(), var.begin());
    }

    long double& operator () (long double row, long double col) {
        return var[row * N + col];
    }

    const long double operator () (long double row, long double col) const {
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

    Matrix<N, M> operator - (const Matrix<N, M>& m) const {
        Matrix<N, M> res;
        for (int r = 0; r < M; ++r) {
            for (int c = 0; c < N; ++c) {
                res(r, c) = (*this)(r, c) - m(r, c);
            }
        }
        return res;
    }

    Vector<M> operator * (const Vector<N>& v) const {
        Vector<M> res;
        for (int r = 0; r < M; ++r) {
            for (int t = 0; t < N; ++t) {
                res(r) += (*this)(r, t) * v(t);
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

    Matrix<N, M> operator * (long double n) const {
        Matrix<N, M> m;
        for (int r = 0; r < M; ++r) {
            for (int c = 0; c < N; ++c) {
                m(r, c) = (*this)(r, c) * n;
            }
        }
        return m;
    }

    Matrix<N, M> operator / (long double n) const {
        assert(n);
        Matrix<N, M> m;
        for (int r = 0; r < M; ++r) {
            for (int c = 0; c < N; ++c) {
                m(r, c) = (*this)(r, c) / n;
            }
        }
        return m;
    }

    template <int K, int L>
    bool operator == (const Matrix<K, L>& m) const {
        if (K != N || L != M) {
            return false;
        }

        for (int r = 0; r < M; ++r) {
            for (int c = 0; c < N; ++c) {
                if ((*this)(r, c) != m(r, c)) {
                    return false;
                }
            }
        }

        return true;
    }

    template <int K, int L>
    bool operator != (const Matrix<K, L>& m) const {
        return !(*this == m);
    }
};

template <int N, int M>
inline Matrix<N, M> operator * (long double n, const Matrix<N, M>& m) {
    return m * n;
}

template <int N, int M>
ostream& operator << (ostream& out, const Matrix<N, M>& m) {
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N - 1; ++c) {
            out << m(r, c) << '\t';
        }
        out << m(r, N - 1) << '\n';
    }
    return out;
}


template <int N>
Matrix<N, N> unit_matrix() {
    Matrix<N, N> unit;
    for (int i = 0; i < N; ++i) {
        unit(i, i) = 1;
    }
    return unit;
}


Matrix<4, 4> matr_inv_3(Matrix<4, 4> j) {
    long double a11 = j(0, 0), a12 = j(0, 1), a13 = j(0, 2),
                a21 = j(1, 0), a22 = j(1, 1), a23 = j(1, 2),
                a31 = j(2, 0), a32 = j(2, 1), a33 = j(2, 2);

    long double det = + (a11 * a22 * a33 - a11 * a23 * a32) \
                      + (a12 * a23 * a31 - a12 * a21 * a33) \
                      + (a13 * a21 * a32 - a13 * a22 * a31);

    Matrix<4, 4> inv;

    inv(0, 0) = a22 * a33 - a23 * a32;
    inv(0, 1) = a13 * a32 - a12 * a33;
    inv(0, 2) = a12 * a23 - a13 * a22;

    inv(1, 0) = a23 * a31 - a21 * a33;
    inv(1, 1) = a11 * a33 - a13 * a31;
    inv(1, 2) = a13 * a21 - a11 * a23;

    inv(2, 0) = a21 * a32 - a22 * a31;
    inv(2, 1) = a12 * a31 - a11 * a32;
    inv(2, 2) = a11 * a22 - a12 * a21;

    inv(3, 3) = det;

    inv(0, 3) = inv(1, 3) = inv(2, 3) = 0;  // not used
    inv(3, 0) = inv(3, 1) = inv(3, 2) = 0;  // not used

    return inv / det;
}

Matrix<4, 4> matr_inv_4(Matrix <4, 4> j) {
    long double a11 = j(0, 0), a12 = j(0, 1), a13 = j(0, 2), a14 = j(0, 3),
                a21 = j(1, 0), a22 = j(1, 1), a23 = j(1, 2), a24 = j(1, 3),
                a31 = j(2, 0), a32 = j(2, 1), a33 = j(2, 2), a34 = j(2, 3),
                a41 = j(3, 0), a42 = j(3, 1), a43 = j(3, 2), a44 = j(3, 3);

    long double det = + (a11 * a22 * a33 * a44 - a11 * a22 * a34 * a43) \
                      + (a11 * a23 * a34 * a42 - a11 * a23 * a32 * a44) \
                      + (a11 * a24 * a32 * a43 - a11 * a24 * a33 * a42) \
                      + (a12 * a21 * a34 * a43 - a12 * a21 * a33 * a44) \
                      + (a12 * a23 * a31 * a44 - a12 * a23 * a34 * a41) \
                      + (a12 * a24 * a33 * a41 - a12 * a24 * a31 * a43) \
                      + (a13 * a21 * a32 * a44 - a13 * a21 * a34 * a42) \
                      + (a13 * a22 * a34 * a41 - a13 * a22 * a31 * a44) \
                      + (a13 * a24 * a31 * a42 - a13 * a24 * a32 * a41) \
                      + (a14 * a21 * a33 * a42 - a14 * a21 * a32 * a43) \
                      + (a14 * a22 * a31 * a43 - a14 * a22 * a33 * a41) \
                      + (a14 * a23 * a32 * a41 - a14 * a23 * a31 * a42);

    Matrix<4, 4> inv;

    inv(0, 0) = + (a22 * a33 * a44 - a22 * a34 * a43) \
                + (a23 * a34 * a42 - a23 * a32 * a44) \
                + (a24 * a32 * a43 - a24 * a33 * a42);
    inv(0, 1) = + (a12 * a34 * a43 - a12 * a33 * a44) \
                + (a13 * a32 * a44 - a13 * a34 * a42) \
                + (a14 * a33 * a42 - a14 * a32 * a43);
    inv(0, 2) = + (a12 * a23 * a44 - a12 * a24 * a43) \
                + (a13 * a24 * a42 - a13 * a22 * a44) \
                + (a14 * a22 * a43 - a14 * a23 * a42);
    inv(0, 3) = + (a12 * a24 * a33 - a12 * a23 * a34) \
                + (a13 * a22 * a34 - a13 * a24 * a32) \
                + (a14 * a23 * a32 - a14 * a22 * a33);
    inv(1, 0) = + (a21 * a34 * a43 - a21 * a33 * a44) \
                + (a23 * a31 * a44 - a23 * a34 * a41) \
                + (a24 * a33 * a41 - a24 * a31 * a43);
    inv(1, 1) = + (a11 * a33 * a44 - a11 * a34 * a43) \
                + (a13 * a34 * a41 - a13 * a31 * a44) \
                + (a14 * a31 * a43 - a14 * a33 * a41);
    inv(1, 2) = + (a11 * a24 * a43 - a11 * a23 * a44) \
                + (a13 * a21 * a44 - a13 * a24 * a41) \
                + (a14 * a23 * a41 - a14 * a21 * a43);
    inv(1, 3) = + (a11 * a23 * a34 - a11 * a24 * a33) \
                + (a13 * a24 * a31 - a13 * a21 * a34) \
                + (a14 * a21 * a33 - a14 * a23 * a31);
    inv(2, 0) = + (a21 * a32 * a44 - a21 * a34 * a42) \
                + (a22 * a34 * a41 - a22 * a31 * a44) \
                + (a24 * a31 * a42 - a24 * a32 * a41);
    inv(2, 1) = + (a11 * a34 * a42 - a11 * a32 * a44) \
                + (a12 * a31 * a44 - a12 * a34 * a41) \
                + (a14 * a32 * a41 - a14 * a31 * a42 );
    inv(2, 2) = + (a11 * a22 * a44 - a11 * a24 * a42) \
                + (a12 * a24 * a41 - a12 * a21 * a44) \
                + (a14 * a21 * a42 - a14 * a22 * a41);
    inv(2, 3) = + (a11 * a24 * a32 - a11 * a22 * a34) \
                + (a12 * a21 * a34 - a12 * a24 * a31) \
                + (a14 * a22 * a31 - a14 * a21 * a32 );
    inv(3, 0) = + (a21 * a33 * a42 - a21 * a32 * a43) \
                + (a22 * a31 * a43 - a22 * a33 * a41) \
                + (a23 * a32 * a41 - a23 * a31 * a42);
    inv(3, 1) = + (a11 * a32 * a43 - a11 * a33 * a42) \
                + (a12 * a33 * a41 - a12 * a31 * a43) \
                + (a13 * a31 * a42 - a13 * a32 * a41);
    inv(3, 2) = + (a11 * a23 * a42 - a11 * a22 * a43) \
                + (a12 * a21 * a43 - a12 * a23 * a41) \
                + (a13 * a22 * a41 - a13 * a21 * a42);
    inv(3, 3) = + (a11 * a22 * a33 - a11 * a23 * a32) \
                + (a12 * a23 * a31 - a12 * a21 * a33) \
                + (a13 * a21 * a32 - a13 * a22 * a31);

    return inv / det;
}

#endif  // LINEAR_ALGEBRA_H_
