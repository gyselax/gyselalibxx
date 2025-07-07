#include "test_math_utils.hpp"

KOKKOS_FUNCTION double cyl_bessel_j(const double nu, const double x)
{
    double fct = 1;
    double sum = 0;
    for (int k = 0; k < 10; fct *= ++k) {
        sum += Kokkos::pow(-1, k) * Kokkos::pow(x / 2, 2 * k + nu)
               / (fct * Kokkos::tgamma(nu + k + 1));
    }
    return sum;
}
