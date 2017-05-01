#ifndef LINEAR_ALGEBRA_TEST_H_
#define LINEAR_ALGEBRA_TEST_H_


#include "linear_algebra.h"


void test_inverse_3() {
    const int num_tests = 6;
    const Matrix<4, 4> test_matr[num_tests] = {
        Matrix<4, 4>({1, 0,  0, 0,  0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}),
        Matrix<4, 4>({2, 0,  0, 0,  0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1}),
        Matrix<4, 4>({3, 0,  0, 0,  0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1}),
        Matrix<4, 4>({1, 1,  1, 0, -1, 2, 1, 0, 1, 4, 1, 0, 0, 0, 0, 1}),
        Matrix<4, 4>({1, 1,  1, 0,  0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1}),
        Matrix<4, 4>({1, 2,  3, 0,  0, 1, 2, 0, 0, 0, 1, 0, 0, 0, 0, 1})
    };

    Matrix<4, 4> m;
    Matrix<4, 4> inv;
    Matrix<4, 4> e = unit_matrix<4>();

    cerr << "testing matr_inv_3...\n";
    for(int i = 0; i < num_tests; ++i) {
        m = test_matr[i];
        inv = matr_inv_3(m);
        if (inv * m != e) {
            cerr << "test " << i + 1 << ": incorrect left inverse\n"
                 << "original matrix:\n" << m << '\n'
                 << "calculated inverse:\n" << inv << '\n'
                 << "inv * m  is:\n" << inv * m << '\n'
                 << "=============================================\n";
        }
        if (m * inv != e) {
            cerr << "test " << i + 1 << ": incorrect right inverse\n"
                 << "original matrix:\n" << m << '\n'
                 << "calculated inverse:\n" << inv << '\n'
                 << "m * inv  is:\n" << m * inv << '\n'
                 << "=============================================\n";
        }
    }
}


void test_inverse_4() {
    const int num_tests = 7;
    const Matrix<4, 4> test_matr[num_tests] = {
        Matrix<4, 4>({1, 0,  0, 0,  0, 1, 0, 0, 0, 0, 1, 0, 0,  0, 0, 1}),
        Matrix<4, 4>({2, 0,  0, 0,  0, 2, 0, 0, 0, 0, 2, 0, 0,  0, 0, 2}),
        Matrix<4, 4>({3, 0,  0, 0,  0, 3, 0, 0, 0, 0, 3, 0, 0,  0, 0, 3}),
        Matrix<4, 4>({0, 0, -1, 0,  0, 0, 0, 1, 1, 0, 0, 0, 0, -1, 0, 0}),
        Matrix<4, 4>({1, 1,  1, 0, -1, 2, 1, 0, 1, 4, 1, 0, 0,  0, 0, 3}),
        Matrix<4, 4>({1, 1,  1, 1,  0, 1, 1, 1, 0, 0, 1, 1, 0,  0, 0, 1}),
        Matrix<4, 4>({1, 2,  3, 4,  0, 1, 2, 3, 0, 0, 1, 2, 0,  0, 0, 1})
    };

    Matrix<4, 4> m;
    Matrix<4, 4> inv;
    Matrix<4, 4> e = unit_matrix<4>();

    cerr << "testing matr_inv_4...\n";
    for(int i = 0; i < num_tests; ++i) {
        m = test_matr[i];
        inv = matr_inv_4(m);
        if (inv * m != e) {
            cerr << "test " << i + 1 << ": incorrect left inverse\n"
                 << "original matrix:\n" << m << '\n'
                 << "calculated inverse:\n" << inv << '\n'
                 << "inv * m  is:\n" << inv * m << '\n'
                 << "=============================================\n";
        }
        if (m * inv != e) {
            cerr << "test " << i + 1 << ": incorrect right inverse\n"
                 << "original matrix:\n" << m << '\n'
                 << "calculated inverse:\n" << inv << '\n'
                 << "m * inv  is:\n" << m * inv << '\n'
                 << "=============================================\n";
        }
    }
}

void linear_algebra_test() {
    test_inverse_3();
    test_inverse_4();
}

#endif  // LINEAR_ALGEBRA_TEST_H_
