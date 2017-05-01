#include <stdio.h>

#include "linear_algebra_test.h"
#include "runge_kutta_test.h"


using std::cerr;


int main() {
    linear_algebra_test();
    runge_kutta_test();
    return 0;
}
