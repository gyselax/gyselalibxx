#include <cassert>
#include <cmath>

#include "test_math_utils.hpp"

double cyl_bessel_j(const double nu, const double x)
{
    assert(nu < 30);
    double sum = std::pow(x / 2, nu) / std::tgamma(nu + 1);
    for (int k = 1; k < 100; ++k) {
        sum += std::pow(-1, k) * std::pow(x / 2, 2 * k + nu)
               / (std::tgamma(k) * std::tgamma(nu + k + 1));
    }
    return sum;
}
