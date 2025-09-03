#include <cassert>
#include <cmath>

#include "test_math_utils.hpp"

double cyl_bessel_j(const double nu, const double x)
{
    assert(x < 30);
    double sum = 0;
    for (int k = 0; k < 100; ++k) {
        sum += std::pow(-1, k) * std::pow(x / 2, 2 * k + nu)
               / (std::tgamma(k + 1) * std::tgamma(nu + k + 1));
    }
    return sum;
}
