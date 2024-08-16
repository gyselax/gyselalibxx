// SPDX-License-Identifier: MIT
#include "test_cases.hpp"

template <>
double ManufacturedPoissonTest<CartesianSolution<CircularToCartesian<X, Y, R, Theta>>>::
        solution_at_pole(Coord<R, Theta> const& coord) const
{
    return 0.0;
}

template <>
double ManufacturedPoissonTest<CartesianSolution<CircularToCartesian<X, Y, R, Theta>>>::
        non_singular_solution(Coord<R, Theta> const& coord) const
{
    const double r = ddc::get<R>(coord);
    const double theta = ddc::get<Theta>(coord);

    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);

    const double tanh_term = std::tanh(20.0 * r - 14.0);
    const double coeff_alpha = std::exp(tanh_term);

    const double plus_pow_6 = ipow(r + 1.0, 6);
    const double minus_pow_6 = ipow(1.0 - r, 6);

    const double cos_sin = std::cos(2.0 * M_PI * r * sin_theta);
    const double cos_cos = std::cos(2.0 * M_PI * r * cos_theta);
    const double sin_sin = std::sin(2.0 * M_PI * r * sin_theta);
    const double sin_cos = std::sin(2.0 * M_PI * r * cos_theta);

    return 0.4096 * minus_pow_6 * plus_pow_6 * coeff_alpha * sin_sin * cos_cos
           - (r * (20.0 * ipow(tanh_term, 2) - 20.0)
                      * (0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta * cos_sin * cos_cos
                         - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_sin * sin_cos * cos_theta
                         + 2.4576 * minus_pow_6 * ipow((r + 1.0), 5) * sin_sin * cos_cos
                         - 2.4576 * ipow((1.0 - r), 5) * plus_pow_6 * sin_sin * cos_cos)
                      * std::exp(-tanh_term)
              + r
                        * ((-1.6384) * (M_PI * M_PI) * minus_pow_6 * plus_pow_6 * ipow(sin_theta, 2)
                                   * sin_sin * cos_cos
                           - 3.2768 * (M_PI * M_PI) * minus_pow_6 * plus_pow_6 * sin_theta * sin_cos
                                     * cos_theta * cos_sin
                           - 1.6384 * (M_PI * M_PI) * minus_pow_6 * plus_pow_6 * sin_sin
                                     * ipow(cos_theta, 2) * cos_cos
                           + 9.8304 * M_PI * minus_pow_6 * ipow(r + 1.0, 5) * sin_theta * cos_sin
                                     * cos_cos
                           - 9.8304 * M_PI * minus_pow_6 * ipow(r + 1.0, 5) * sin_sin * sin_cos
                                     * cos_theta
                           + 12.288 * minus_pow_6 * ipow(r + 1.0, 4) * sin_sin * cos_cos
                           - 9.8304 * M_PI * ipow(1.0 - r, 5) * plus_pow_6 * sin_theta * cos_sin
                                     * cos_cos
                           + 9.8304 * M_PI * ipow(1.0 - r, 5) * plus_pow_6 * sin_sin * sin_cos
                                     * cos_theta
                           - 29.4912 * ipow(1.0 - r, 5) * ipow(r + 1.0, 5) * sin_sin * cos_cos
                           + 12.288 * ipow(1.0 - r, 4) * plus_pow_6 * sin_sin * cos_cos)
                        * std::exp(-tanh_term)
              + (0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta * cos_sin * cos_cos
                 - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_sin * sin_cos * cos_theta
                 + 2.4576 * minus_pow_6 * ipow(r + 1.0, 5) * sin_sin * cos_cos
                 - 2.4576 * ipow(1.0 - r, 5) * plus_pow_6 * sin_sin * cos_cos)
                        * std::exp(-tanh_term)
              + ((-1.6384) * (M_PI * M_PI) * r * minus_pow_6 * plus_pow_6 * ipow(sin_theta, 2)
                         * sin_sin * cos_cos
                 + 3.2768 * (M_PI * M_PI) * r * minus_pow_6 * plus_pow_6 * sin_theta * sin_cos
                           * cos_theta * cos_sin
                 - 1.6384 * (M_PI * M_PI) * r * minus_pow_6 * plus_pow_6 * sin_sin
                           * ipow(cos_theta, 2) * cos_cos
                 - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta * cos_sin * cos_cos
                 + 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_sin * sin_cos * cos_theta)
                        * std::exp(-tanh_term))
                     / r;
}

//---------------------------------------------------------------------

template <>
double ManufacturedPoissonTest<CartesianSolution<CzarnyToCartesian<X, Y, R, Theta>>>::
        solution_at_pole(Coord<R, Theta> const& coord) const
{
    return 0.0;
}

template <>
double ManufacturedPoissonTest<CartesianSolution<CzarnyToCartesian<X, Y, R, Theta>>>::
        non_singular_solution(Coord<R, Theta> const& coord) const
{
    const double r = ddc::get<R>(coord);
    const double theta = ddc::get<Theta>(coord);
    const double epsilon = m_coordinate_converter.epsilon();
    const double e = m_coordinate_converter.e();

    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);

    const double tanh_term = std::tanh(20.0 * r - 14.0);
    const double coeff_alpha = std::exp(tanh_term);

    const double plus_pow_6 = ipow(r + 1.0, 6);
    const double minus_pow_6 = ipow(1.0 - r, 6);

    const double xi = 1. / std::sqrt(1. - epsilon * epsilon * 0.25);
    const double tmp1 = std::sqrt(epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0);

    return 0.4096 * minus_pow_6 * plus_pow_6 * coeff_alpha
                   * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
           - ((-r)
                      * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                          + e * cos_theta * xi / ((2.0 - tmp1)))
                                 * (e * epsilon * r * sin_theta * cos_theta * xi
                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                    + e * sin_theta * xi / ((2.0 - tmp1)))
                         - sin_theta * cos_theta
                                   / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                      * (0.4096 * minus_pow_6 * plus_pow_6
                                 * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                    + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                 * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                 * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                         - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta
                                   * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                   / tmp1)
                      * (20.0 * ipow(tanh_term, 2) - 20.0) * std::exp(-tanh_term)
                      / std::sqrt(
                              (-ipow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                   + 1.0)),
                                     2))
                              + (ipow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow((2.0 - tmp1), 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                      2)
                                 + ipow(sin_theta, 2)
                                           / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                        * (ipow((e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow((2.0 - tmp1), 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1))),
                                                2)
                                           + ipow(cos_theta, 2)
                                                     / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0)))
              - r
                        * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                    / (ipow((2.0 - tmp1), 2) * tmp1)
                            + e * cos_theta * xi / ((2.0 - tmp1)))
                                   * (e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow((2.0 - tmp1), 2) * tmp1)
                                      + e * sin_theta * xi / ((2.0 - tmp1)))
                           - sin_theta * cos_theta
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (0.4096 * minus_pow_6 * plus_pow_6
                                   * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                              / (ipow((2.0 - tmp1), 2) * tmp1)
                                      + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1)
                        * (1.0 / 2.0
                                   * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow((2.0 - tmp1), 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow((2.0 - tmp1), 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                   * (4.0 * epsilon * sin_theta * ipow(cos_theta, 2)
                                              / ipow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                      + 1.0),
                                                     2)
                                      + 2.0
                                                * ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                           / (ipow((2.0 - tmp1), 2) * tmp1)
                                                   + e * cos_theta * xi / ((2.0 - tmp1)))
                                                * ((-e) * (epsilon * epsilon) * r * sin_theta
                                                           * ipow(cos_theta, 2) * xi
                                                           / (ipow(2.0 - tmp1, 2)
                                                              * pow((epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0),
                                                                    (3.0 / 2.0)))
                                                   + 2.0 * e * (epsilon * epsilon) * r * sin_theta
                                                             * ipow(cos_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 3)
                                                                * (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))
                                                   + 2.0 * e * epsilon * sin_theta * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1))
                                      + 2.0
                                                * (e * epsilon * r * sin_theta * cos_theta * xi
                                                           / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * sin_theta * xi / ((2.0 - tmp1)))
                                                * (e * (epsilon * epsilon) * r * ipow(sin_theta, 2)
                                                           * cos_theta * xi
                                                           / (ipow(2.0 - tmp1, 2)
                                                              * pow((epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0),
                                                                    (3.0 / 2.0)))
                                                   - 2.0 * e * (epsilon * epsilon) * r
                                                             * ipow(sin_theta, 2) * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 3)
                                                                * (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))
                                                   - e * epsilon * ipow(sin_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * epsilon * ipow(cos_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)))
                           - 1.0 / 2.0
                                     * ((-2.0) * epsilon * ipow(cos_theta, 3)
                                                / ipow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0),
                                                       2)
                                        + (e * epsilon * r * sin_theta * cos_theta * xi
                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                  * ((-2.0) * e * (epsilon * epsilon) * r
                                                             * sin_theta * ipow(cos_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2)
                                                                * pow((epsilon
                                                                               * (epsilon
                                                                                  + 2.0 * r * cos_theta)
                                                                       + 1.0),
                                                                      (3.0 / 2.0)))
                                                     + 4.0 * e * (epsilon * epsilon) * r * sin_theta
                                                               * ipow(cos_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 3)
                                                                  * (epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0))
                                                     + 4.0 * e * epsilon * sin_theta * cos_theta
                                                               * xi / (ipow(2.0 - tmp1, 2) * tmp1)))
                                     * (ipow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                      / (ipow(2.0 - tmp1, 2) * tmp1)
                                              + e * cos_theta * xi / ((2.0 - tmp1))),
                                             2)
                                        + ipow(sin_theta, 2)
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0))
                           - 1.0 / 2.0
                                     * ((-2.0) * epsilon * ipow(sin_theta, 2) * cos_theta
                                                / ipow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0),
                                                       2)
                                        + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                           + e * cos_theta * xi / ((2.0 - tmp1)))
                                                  * (2.0 * e * (epsilon * epsilon) * r
                                                             * ipow(sin_theta, 2) * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2)
                                                                * pow((epsilon
                                                                               * (epsilon
                                                                                  + 2.0 * r * cos_theta)
                                                                       + 1.0),
                                                                      (3.0 / 2.0)))
                                                     - 4.0 * e * (epsilon * epsilon) * r
                                                               * ipow(sin_theta, 2) * cos_theta * xi
                                                               / (ipow(2.0 - tmp1, 3)
                                                                  * (epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0))
                                                     - 2.0 * e * epsilon * ipow(sin_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + 2.0 * e * epsilon * ipow(cos_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)))
                                     * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                             + e * sin_theta * xi / ((2.0 - tmp1))),
                                            2.0)
                                        + ipow(cos_theta, 2)
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0)))
                        * std::exp(-tanh_term)
                        / pow(((-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                   + 1.0)),
                                     2.0))
                               + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                      2.0)
                                  + ipow(sin_theta, 2)
                                            / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                         * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1))),
                                                2.0)
                                            + ipow(cos_theta, 2)
                                                      / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                         + 1.0))),
                              (3.0 / 2.0))
              - r
                        * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + e * cos_theta * xi / ((2.0 - tmp1)))
                                   * (e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + e * sin_theta * xi / ((2.0 - tmp1)))
                           - sin_theta * cos_theta
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (0.8192 * M_PI * epsilon * minus_pow_6 * plus_pow_6 * sin_theta
                                   * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                   * cos_theta
                                   / pow((epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0),
                                         (3.0 / 2.0))
                           - 0.4096 * minus_pow_6 * plus_pow_6
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           + 0.4096 * minus_pow_6 * plus_pow_6
                                     * (2.0 * M_PI * e * (epsilon * epsilon) * r
                                                * ipow(sin_theta, 2) * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        - 4.0 * M_PI * e * (epsilon * epsilon) * r
                                                  * ipow(sin_theta, 2) * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        - 2.0 * M_PI * e * epsilon * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * epsilon * ipow(cos_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 1.6384 * (M_PI * M_PI) * minus_pow_6 * plus_pow_6 * sin_theta
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)
                           + 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon) * cos_theta
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           - 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                                     * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                     * sin_theta * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           + 2.4576 * minus_pow_6 * ipow(r + 1.0, 5)
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 4.9152 * M_PI * minus_pow_6 * ipow(r + 1.0, 5) * sin_theta
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           - 2.4576 * ipow(1.0 - r, 5) * plus_pow_6
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 4.9152 * M_PI * ipow(1.0 - r, 5) * plus_pow_6 * sin_theta
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1)
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              + r
                        * ((-2.0) * epsilon * ipow(sin_theta, 2) * cos_theta
                                   / pow((epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0), 2.0)
                           + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                      / (ipow(2.0 - tmp1, 2) * tmp1)
                              + e * cos_theta * xi / ((2.0 - tmp1)))
                                     * (2.0 * e * (epsilon * epsilon) * r * ipow(sin_theta, 2)
                                                * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        - 4.0 * e * (epsilon * epsilon) * r * ipow(sin_theta, 2)
                                                  * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        - 2.0 * e * epsilon * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * e * epsilon * ipow(cos_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)))
                        * (0.4096 * minus_pow_6 * plus_pow_6
                                   * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta / tmp1
                           + 2.4576 * minus_pow_6 * ipow(r + 1.0, 5)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           - 2.4576 * ipow(1.0 - r, 5) * plus_pow_6
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon))
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              - r
                        * (0.4096 * minus_pow_6 * plus_pow_6
                                   * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1)
                        * (2.0 * epsilon * sin_theta * ipow(cos_theta, 2)
                                   / pow((epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0), 2.0)
                           + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                      / (ipow(2.0 - tmp1, 2) * tmp1)
                              + e * cos_theta * xi / ((2.0 - tmp1)))
                                     * ((-e) * (epsilon * epsilon) * r * sin_theta
                                                * ipow(cos_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        + 2.0 * e * (epsilon * epsilon) * r * sin_theta
                                                  * ipow(cos_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        + 2.0 * e * epsilon * sin_theta * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1))
                           + (e * epsilon * r * sin_theta * cos_theta * xi
                                      / (ipow(2.0 - tmp1, 2) * tmp1)
                              + e * sin_theta * xi / ((2.0 - tmp1)))
                                     * (e * (epsilon * epsilon) * r * ipow(sin_theta, 2) * cos_theta
                                                * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        - 2.0 * e * (epsilon * epsilon) * r * ipow(sin_theta, 2)
                                                  * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        - e * epsilon * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * epsilon * ipow(cos_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)))
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              + r
                        * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                + e * cos_theta * xi / ((2.0 - tmp1))),
                               2.0)
                           + ipow(sin_theta, 2) / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (20.0 * pow(tanh_term, 2.0) - 20.0)
                        * (0.4096 * minus_pow_6 * plus_pow_6
                                   * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta / tmp1
                           + 2.4576 * minus_pow_6 * ipow(r + 1.0, 5)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           - 2.4576 * ipow(1.0 - r, 5) * plus_pow_6
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon))
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              + r
                        * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                + e * cos_theta * xi / ((2.0 - tmp1))),
                               2.0)
                           + ipow(sin_theta, 2) / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * ((-0.8192) * M_PI * epsilon * minus_pow_6 * plus_pow_6
                                   * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                   * ipow(cos_theta, 2)
                                   / pow((epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0),
                                         (3.0 / 2.0))
                           - 0.4096 * minus_pow_6 * plus_pow_6
                                     * pow((2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta
                                                    * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                            + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1))),
                                           2.0)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           + 0.4096 * minus_pow_6 * plus_pow_6
                                     * ((-2.0) * M_PI * e * (epsilon * epsilon) * r * sin_theta
                                                * ipow(cos_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        + 4.0 * M_PI * e * (epsilon * epsilon) * r * sin_theta
                                                  * ipow(cos_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        + 4.0 * M_PI * e * epsilon * sin_theta * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 1.6384 * (M_PI * M_PI) * minus_pow_6 * plus_pow_6
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * ipow(cos_theta, 2)
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)
                           + 1.6384 * M_PI * minus_pow_6 * plus_pow_6
                                     * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon) * cos_theta
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           + 4.9152 * minus_pow_6 * ipow(r + 1.0, 5)
                                     * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 9.8304 * M_PI * minus_pow_6 * ipow(r + 1.0, 5)
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta / tmp1
                           + 12.288 * minus_pow_6 * ipow(r + 1.0, 4)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           - 4.9152 * ipow(1.0 - r, 5) * plus_pow_6
                                     * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 9.8304 * M_PI * ipow(1.0 - r, 5) * plus_pow_6
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta / tmp1
                           - 29.4912 * ipow(1.0 - r, 5) * ipow(r + 1.0, 5)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           + 12.288 * ipow(1.0 - r, 4) * plus_pow_6
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon))
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              + r
                        * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                + e * cos_theta * xi / ((2.0 - tmp1))),
                               2.0)
                           + ipow(sin_theta, 2) / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (1.0 / 2.0
                                   * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                   * (4.0 * epsilon * sin_theta * ipow(cos_theta, 2)
                                              / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0),
                                                    2.0)
                                      + 2.0
                                                * ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                           / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * cos_theta * xi / ((2.0 - tmp1)))
                                                * ((-e) * (epsilon * epsilon) * r * sin_theta
                                                           * ipow(cos_theta, 2) * xi
                                                           / (ipow(2.0 - tmp1, 2)
                                                              * pow((epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0),
                                                                    (3.0 / 2.0)))
                                                   + 2.0 * e * (epsilon * epsilon) * r * sin_theta
                                                             * ipow(cos_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 3)
                                                                * (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))
                                                   + 2.0 * e * epsilon * sin_theta * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1))
                                      + 2.0
                                                * (e * epsilon * r * sin_theta * cos_theta * xi
                                                           / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * sin_theta * xi / ((2.0 - tmp1)))
                                                * (e * (epsilon * epsilon) * r * ipow(sin_theta, 2)
                                                           * cos_theta * xi
                                                           / (ipow(2.0 - tmp1, 2)
                                                              * pow((epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0),
                                                                    (3.0 / 2.0)))
                                                   - 2.0 * e * (epsilon * epsilon) * r
                                                             * ipow(sin_theta, 2) * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 3)
                                                                * (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))
                                                   - e * epsilon * ipow(sin_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * epsilon * ipow(cos_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)))
                           - 1.0 / 2.0
                                     * ((-2.0) * epsilon * ipow(cos_theta, 3)
                                                / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                       + 1.0),
                                                      2.0)
                                        + (e * epsilon * r * sin_theta * cos_theta * xi
                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                  * ((-2.0) * e * (epsilon * epsilon) * r
                                                             * sin_theta * ipow(cos_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2)
                                                                * pow((epsilon
                                                                               * (epsilon
                                                                                  + 2.0 * r * cos_theta)
                                                                       + 1.0),
                                                                      (3.0 / 2.0)))
                                                     + 4.0 * e * (epsilon * epsilon) * r * sin_theta
                                                               * ipow(cos_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 3)
                                                                  * (epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0))
                                                     + 4.0 * e * epsilon * sin_theta * cos_theta
                                                               * xi / (ipow(2.0 - tmp1, 2) * tmp1)))
                                     * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                             + e * cos_theta * xi / ((2.0 - tmp1))),
                                            2.0)
                                        + ipow(sin_theta, 2)
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0))
                           - 1.0 / 2.0
                                     * ((-2.0) * epsilon * ipow(sin_theta, 2) * cos_theta
                                                / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                       + 1.0),
                                                      2.0)
                                        + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                           + e * cos_theta * xi / ((2.0 - tmp1)))
                                                  * (2.0 * e * (epsilon * epsilon) * r
                                                             * ipow(sin_theta, 2) * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2)
                                                                * pow((epsilon
                                                                               * (epsilon
                                                                                  + 2.0 * r * cos_theta)
                                                                       + 1.0),
                                                                      (3.0 / 2.0)))
                                                     - 4.0 * e * (epsilon * epsilon) * r
                                                               * ipow(sin_theta, 2) * cos_theta * xi
                                                               / (ipow(2.0 - tmp1, 3)
                                                                  * (epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0))
                                                     - 2.0 * e * epsilon * ipow(sin_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + 2.0 * e * epsilon * ipow(cos_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)))
                                     * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                             + e * sin_theta * xi / ((2.0 - tmp1))),
                                            2.0)
                                        + ipow(cos_theta, 2)
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0)))
                        * (0.4096 * minus_pow_6 * plus_pow_6
                                   * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta / tmp1
                           + 2.4576 * minus_pow_6 * ipow(r + 1.0, 5)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           - 2.4576 * ipow(1.0 - r, 5) * plus_pow_6
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon))
                        * std::exp(-tanh_term)
                        / pow(((-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                   + 1.0)),
                                     2.0))
                               + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                      2.0)
                                  + ipow(sin_theta, 2)
                                            / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                         * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1))),
                                                2.0)
                                            + ipow(cos_theta, 2)
                                                      / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                         + 1.0))),
                              (3.0 / 2.0))
              - (((-e) * epsilon * r * ipow(sin_theta, 2) * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                  + e * cos_theta * xi / ((2.0 - tmp1)))
                         * (e * epsilon * r * sin_theta * cos_theta * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + e * sin_theta * xi / ((2.0 - tmp1)))
                 - sin_theta * cos_theta / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (0.4096 * minus_pow_6 * plus_pow_6
                                   * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1)
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              - (((-e) * epsilon * r * ipow(sin_theta, 2) * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                  + e * cos_theta * xi / ((2.0 - tmp1)))
                         * (e * epsilon * r * sin_theta * cos_theta * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + e * sin_theta * xi / ((2.0 - tmp1)))
                 - sin_theta * cos_theta / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (0.8192 * M_PI * epsilon * r * minus_pow_6 * plus_pow_6 * sin_theta
                                   * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                   * cos_theta
                                   / pow((epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0),
                                         (3.0 / 2.0))
                           - 0.4096 * r * minus_pow_6 * plus_pow_6
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           + 1.6384 * (M_PI * M_PI) * r * minus_pow_6 * plus_pow_6 * sin_theta
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)
                           + 0.8192 * M_PI * r * minus_pow_6 * plus_pow_6
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon) * cos_theta
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           - 0.8192 * M_PI * r * minus_pow_6 * plus_pow_6
                                     * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                     * sin_theta * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           + 2.4576 * r * minus_pow_6 * ipow(r + 1.0, 5)
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 4.9152 * M_PI * r * minus_pow_6 * ipow(r + 1.0, 5) * sin_theta
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           - 2.4576 * r * ipow(1.0 - r, 5) * plus_pow_6
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 4.9152 * M_PI * r * ipow(1.0 - r, 5) * plus_pow_6 * sin_theta
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           + 0.4096 * minus_pow_6 * plus_pow_6
                                     * (2.0 * M_PI * e * (epsilon * epsilon) * (r * r)
                                                * ipow(sin_theta, 2) * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        - 4.0 * M_PI * e * (epsilon * epsilon) * (r * r)
                                                  * ipow(sin_theta, 2) * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        - 4.0 * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * epsilon * r * ipow(cos_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1)
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              - (((-e) * epsilon * r * ipow(sin_theta, 2) * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                  + e * cos_theta * xi / ((2.0 - tmp1)))
                         * (e * epsilon * r * sin_theta * cos_theta * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + e * sin_theta * xi / ((2.0 - tmp1)))
                 - sin_theta * cos_theta / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (1.0 / 2.0
                                   * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                   * ((-4.0) * epsilon * r * ipow(sin_theta, 2) * cos_theta
                                              / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0),
                                                    2.0)
                                      + 2.0
                                                * ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                           / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * cos_theta * xi / ((2.0 - tmp1)))
                                                * (e * (epsilon * epsilon) * (r * r)
                                                           * ipow(sin_theta, 2) * cos_theta * xi
                                                           / (ipow(2.0 - tmp1, 2)
                                                              * pow((epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0),
                                                                    (3.0 / 2.0)))
                                                   - 2.0 * e * (epsilon * epsilon) * (r * r)
                                                             * ipow(sin_theta, 2) * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 3)
                                                                * (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))
                                                   - 2.0 * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * epsilon * r * ipow(cos_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * cos_theta * xi / ((2.0 - tmp1)))
                                      + 2.0
                                                * (e * epsilon * r * sin_theta * cos_theta * xi
                                                           / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * sin_theta * xi / ((2.0 - tmp1)))
                                                * ((-e) * (epsilon * epsilon) * (r * r)
                                                           * ipow(sin_theta, 3) * xi
                                                           / (ipow(2.0 - tmp1, 2)
                                                              * pow((epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0),
                                                                    (3.0 / 2.0)))
                                                   + 2.0 * e * (epsilon * epsilon) * (r * r)
                                                             * ipow(sin_theta, 3) * xi
                                                             / (ipow(2.0 - tmp1, 3)
                                                                * (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))
                                                   - 3.0 * e * epsilon * r * sin_theta * cos_theta
                                                             * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   - e * sin_theta * xi / ((2.0 - tmp1)))
                                      + 2.0 * ipow(sin_theta, 2)
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)
                                      - 2.0 * ipow(cos_theta, 2)
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                           - 1.0 / 2.0
                                     * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                             + e * cos_theta * xi / ((2.0 - tmp1))),
                                            2.0)
                                        + ipow(sin_theta, 2)
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0))
                                     * (2.0 * epsilon * r * sin_theta * ipow(cos_theta, 2)
                                                / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                       + 1.0),
                                                      2.0)
                                        + (e * epsilon * r * sin_theta * cos_theta * xi
                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                  * (2.0 * e * (epsilon * epsilon) * (r * r)
                                                             * ipow(sin_theta, 2) * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2)
                                                                * pow((epsilon
                                                                               * (epsilon
                                                                                  + 2.0 * r * cos_theta)
                                                                       + 1.0),
                                                                      (3.0 / 2.0)))
                                                     - 4.0 * e * (epsilon * epsilon) * (r * r)
                                                               * ipow(sin_theta, 2) * cos_theta * xi
                                                               / (ipow(2.0 - tmp1, 3)
                                                                  * (epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0))
                                                     - 4.0 * e * epsilon * r * ipow(sin_theta, 2)
                                                               * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + 2.0 * e * epsilon * r * ipow(cos_theta, 2)
                                                               * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + 2.0 * e * cos_theta * xi / ((2.0 - tmp1)))
                                        - 2.0 * sin_theta * cos_theta
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0))
                           - 1.0 / 2.0
                                     * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                             + e * sin_theta * xi / ((2.0 - tmp1))),
                                            2.0)
                                        + ipow(cos_theta, 2)
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0))
                                     * (2.0 * epsilon * r * ipow(sin_theta, 3)
                                                / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                       + 1.0),
                                                      2.0)
                                        + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                           + e * cos_theta * xi / ((2.0 - tmp1)))
                                                  * ((-2.0) * e * (epsilon * epsilon) * (r * r)
                                                             * ipow(sin_theta, 3) * xi
                                                             / (ipow(2.0 - tmp1, 2)
                                                                * pow((epsilon
                                                                               * (epsilon
                                                                                  + 2.0 * r * cos_theta)
                                                                       + 1.0),
                                                                      (3.0 / 2.0)))
                                                     + 4.0 * e * (epsilon * epsilon) * (r * r)
                                                               * ipow(sin_theta, 3) * xi
                                                               / (ipow(2.0 - tmp1, 3)
                                                                  * (epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0))
                                                     - 6.0 * e * epsilon * r * sin_theta * cos_theta
                                                               * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     - 2.0 * e * sin_theta * xi / ((2.0 - tmp1)))
                                        + 2.0 * sin_theta * cos_theta
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0)))
                        * (0.4096 * minus_pow_6 * plus_pow_6
                                   * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta / tmp1
                           + 2.4576 * minus_pow_6 * ipow(r + 1.0, 5)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           - 2.4576 * ipow(1.0 - r, 5) * plus_pow_6
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon))
                        * std::exp(-tanh_term)
                        / pow(((-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                   + 1.0)),
                                     2.0))
                               + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                      2.0)
                                  + ipow(sin_theta, 2)
                                            / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                         * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1))),
                                                2.0)
                                            + ipow(cos_theta, 2)
                                                      / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                         + 1.0))),
                              (3.0 / 2.0))
              + (0.4096 * minus_pow_6 * plus_pow_6
                         * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                         * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                         * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                 - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta
                           * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1))) / tmp1)
                        * (2.0 * epsilon * r * sin_theta * ipow(cos_theta, 2)
                                   / pow((epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0), 2.0)
                           + (e * epsilon * r * sin_theta * cos_theta * xi
                                      / (ipow(2.0 - tmp1, 2) * tmp1)
                              + e * sin_theta * xi / ((2.0 - tmp1)))
                                     * (2.0 * e * (epsilon * epsilon) * (r * r) * ipow(sin_theta, 2)
                                                * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        - 4.0 * e * (epsilon * epsilon) * (r * r)
                                                  * ipow(sin_theta, 2) * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        - 4.0 * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * e * epsilon * r * ipow(cos_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * e * cos_theta * xi / ((2.0 - tmp1)))
                           - 2.0 * sin_theta * cos_theta
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              + (0.4096 * minus_pow_6 * plus_pow_6
                         * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                         * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                         * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                 - 0.8192 * M_PI * minus_pow_6 * plus_pow_6 * sin_theta
                           * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1))) / tmp1)
                        * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                + e * sin_theta * xi / ((2.0 - tmp1))),
                               2.0)
                           + ipow(cos_theta, 2) / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (1.0 / 2.0
                                   * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                   * ((-4.0) * epsilon * r * ipow(sin_theta, 2) * cos_theta
                                              / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0),
                                                    2.0)
                                      + 2.0
                                                * ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                           / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * cos_theta * xi / ((2.0 - tmp1)))
                                                * (e * (epsilon * epsilon) * (r * r)
                                                           * ipow(sin_theta, 2) * cos_theta * xi
                                                           / (ipow(2.0 - tmp1, 2)
                                                              * pow((epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0),
                                                                    (3.0 / 2.0)))
                                                   - 2.0 * e * (epsilon * epsilon) * (r * r)
                                                             * ipow(sin_theta, 2) * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 3)
                                                                * (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))
                                                   - 2.0 * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * epsilon * r * ipow(cos_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * cos_theta * xi / ((2.0 - tmp1)))
                                      + 2.0
                                                * (e * epsilon * r * sin_theta * cos_theta * xi
                                                           / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   + e * sin_theta * xi / ((2.0 - tmp1)))
                                                * ((-e) * (epsilon * epsilon) * (r * r)
                                                           * ipow(sin_theta, 3) * xi
                                                           / (ipow(2.0 - tmp1, 2)
                                                              * pow((epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0),
                                                                    (3.0 / 2.0)))
                                                   + 2.0 * e * (epsilon * epsilon) * (r * r)
                                                             * ipow(sin_theta, 3) * xi
                                                             / (ipow(2.0 - tmp1, 3)
                                                                * (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))
                                                   - 3.0 * e * epsilon * r * sin_theta * cos_theta
                                                             * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                   - e * sin_theta * xi / ((2.0 - tmp1)))
                                      + 2.0 * ipow(sin_theta, 2)
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)
                                      - 2.0 * ipow(cos_theta, 2)
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                           - 1.0 / 2.0
                                     * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                             + e * cos_theta * xi / ((2.0 - tmp1))),
                                            2.0)
                                        + ipow(sin_theta, 2)
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0))
                                     * (2.0 * epsilon * r * sin_theta * ipow(cos_theta, 2)
                                                / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                       + 1.0),
                                                      2.0)
                                        + (e * epsilon * r * sin_theta * cos_theta * xi
                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                  * (2.0 * e * (epsilon * epsilon) * (r * r)
                                                             * ipow(sin_theta, 2) * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2)
                                                                * pow((epsilon
                                                                               * (epsilon
                                                                                  + 2.0 * r * cos_theta)
                                                                       + 1.0),
                                                                      (3.0 / 2.0)))
                                                     - 4.0 * e * (epsilon * epsilon) * (r * r)
                                                               * ipow(sin_theta, 2) * cos_theta * xi
                                                               / (ipow(2.0 - tmp1, 3)
                                                                  * (epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0))
                                                     - 4.0 * e * epsilon * r * ipow(sin_theta, 2)
                                                               * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + 2.0 * e * epsilon * r * ipow(cos_theta, 2)
                                                               * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + 2.0 * e * cos_theta * xi / ((2.0 - tmp1)))
                                        - 2.0 * sin_theta * cos_theta
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0))
                           - 1.0 / 2.0
                                     * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                             + e * sin_theta * xi / ((2.0 - tmp1))),
                                            2.0)
                                        + ipow(cos_theta, 2)
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0))
                                     * (2.0 * epsilon * r * ipow(sin_theta, 3)
                                                / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                       + 1.0),
                                                      2.0)
                                        + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                           + e * cos_theta * xi / ((2.0 - tmp1)))
                                                  * ((-2.0) * e * (epsilon * epsilon) * (r * r)
                                                             * ipow(sin_theta, 3) * xi
                                                             / (ipow(2.0 - tmp1, 2)
                                                                * pow((epsilon
                                                                               * (epsilon
                                                                                  + 2.0 * r * cos_theta)
                                                                       + 1.0),
                                                                      (3.0 / 2.0)))
                                                     + 4.0 * e * (epsilon * epsilon) * (r * r)
                                                               * ipow(sin_theta, 3) * xi
                                                               / (ipow(2.0 - tmp1, 3)
                                                                  * (epsilon
                                                                             * (epsilon
                                                                                + 2.0 * r * cos_theta)
                                                                     + 1.0))
                                                     - 6.0 * e * epsilon * r * sin_theta * cos_theta
                                                               * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     - 2.0 * e * sin_theta * xi / ((2.0 - tmp1)))
                                        + 2.0 * sin_theta * cos_theta
                                                  / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                     + 1.0)))
                        * std::exp(-tanh_term)
                        / pow(((-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1)))
                                              * (e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1)))
                                      - sin_theta * cos_theta
                                                / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                   + 1.0)),
                                     2.0))
                               + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                      2.0)
                                  + ipow(sin_theta, 2)
                                            / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                         * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * sin_theta * xi / ((2.0 - tmp1))),
                                                2.0)
                                            + ipow(cos_theta, 2)
                                                      / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                         + 1.0))),
                              (3.0 / 2.0))
              + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                      + e * cos_theta * xi / ((2.0 - tmp1))),
                     2.0)
                 + ipow(sin_theta, 2) / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * (0.4096 * minus_pow_6 * plus_pow_6
                                   * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                   * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           + 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta / tmp1
                           + 2.4576 * minus_pow_6 * ipow(r + 1.0, 5)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           - 2.4576 * ipow(1.0 - r, 5) * plus_pow_6
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon))
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              + (pow((e * epsilon * r * sin_theta * cos_theta * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                      + e * sin_theta * xi / ((2.0 - tmp1))),
                     2.0)
                 + ipow(cos_theta, 2) / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * ((-0.8192) * M_PI * epsilon * r * minus_pow_6 * plus_pow_6
                                   * ipow(sin_theta, 2)
                                   * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                   * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                   / pow((epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0),
                                         (3.0 / 2.0))
                           - 0.4096 * r * minus_pow_6 * plus_pow_6
                                     * pow(((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2)
                                                    * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                            + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1))),
                                           2.0)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           - 1.6384 * (M_PI * M_PI) * r * minus_pow_6 * plus_pow_6
                                     * ipow(sin_theta, 2)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)
                           - 1.6384 * M_PI * r * minus_pow_6 * plus_pow_6
                                     * ((-2.0) * M_PI * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + 2.0 * M_PI * e * cos_theta * xi / ((2.0 - tmp1)))
                                     * sin_theta * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     / tmp1
                           + 0.4096 * minus_pow_6 * plus_pow_6
                                     * ((-2.0) * M_PI * e * (epsilon * epsilon) * (r * r)
                                                * ipow(sin_theta, 3) * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        + 4.0 * M_PI * e * (epsilon * epsilon) * (r * r)
                                                  * ipow(sin_theta, 3) * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        - 6.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        - 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                                     * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::cos(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           - 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                                     * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                                     * std::sin(
                                             2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                                     * cos_theta / tmp1)
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0)))
              - (0.4096 * minus_pow_6 * plus_pow_6
                         * (2.0 * M_PI * e * epsilon * r * sin_theta * cos_theta * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + 2.0 * M_PI * e * sin_theta * xi / ((2.0 - tmp1)))
                         * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                         * std::cos(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                 + 0.8192 * M_PI * minus_pow_6 * plus_pow_6
                           * std::sin(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                           * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           * cos_theta / tmp1
                 + 2.4576 * minus_pow_6 * ipow(r + 1.0, 5)
                           * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon)
                 - 2.4576 * ipow(1.0 - r, 5) * plus_pow_6
                           * std::sin(2.0 * M_PI * e * r * sin_theta * xi / ((2.0 - tmp1)))
                           * std::cos(2.0 * M_PI * (1.0 - tmp1) / epsilon))
                        * ((-2.0) * epsilon * r * ipow(sin_theta, 2) * cos_theta
                                   / pow((epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0), 2.0)
                           + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                      / (ipow(2.0 - tmp1, 2) * tmp1)
                              + e * cos_theta * xi / ((2.0 - tmp1)))
                                     * (e * (epsilon * epsilon) * (r * r) * ipow(sin_theta, 2)
                                                * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        - 2.0 * e * (epsilon * epsilon) * (r * r)
                                                  * ipow(sin_theta, 2) * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        - 2.0 * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * epsilon * r * ipow(cos_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                           + (e * epsilon * r * sin_theta * cos_theta * xi
                                      / (ipow(2.0 - tmp1, 2) * tmp1)
                              + e * sin_theta * xi / ((2.0 - tmp1)))
                                     * ((-e) * (epsilon * epsilon) * (r * r) * ipow(sin_theta, 3)
                                                * xi
                                                / (ipow(2.0 - tmp1, 2)
                                                   * pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0),
                                                         (3.0 / 2.0)))
                                        + 2.0 * e * (epsilon * epsilon) * (r * r)
                                                  * ipow(sin_theta, 3) * xi
                                                  / (ipow(2.0 - tmp1, 3)
                                                     * (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                        - 3.0 * e * epsilon * r * sin_theta * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                        - e * sin_theta * xi / ((2.0 - tmp1)))
                           + ipow(sin_theta, 2) / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)
                           - ipow(cos_theta, 2) / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                        * std::exp(-tanh_term)
                        / std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))))
                     / (r
                        * std::sqrt(
                                (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1)))
                                       - sin_theta * cos_theta
                                                 / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0)),
                                      2.0))
                                + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1))),
                                       2.0)
                                   + ipow(sin_theta, 2)
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                          * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * sin_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(cos_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))));
}

//---------------------------------------------------------------------

template <>
double ManufacturedPoissonTest<CurvilinearSolution<CircularToCartesian<X, Y, R, Theta>>>::
        solution_at_pole(Coord<R, Theta> const& coord) const
{
    return 0.0;
}

template <>
double ManufacturedPoissonTest<CurvilinearSolution<CircularToCartesian<X, Y, R, Theta>>>::
        non_singular_solution(Coord<R, Theta> const& coord) const
{
    const double r = ddc::get<R>(coord);
    const double theta = ddc::get<Theta>(coord);

    const double cos_11_theta = std::cos(11.0 * theta);

    const double tanh_term = std::tanh(20.0 * r - 14.0);
    const double coeff_alpha = std::exp(tanh_term);

    const double minus_pow_6 = ipow(r - 1.0, 6);

    return 0.4096 * ipow(r, 6) * minus_pow_6 * coeff_alpha * cos_11_theta
           - ipow(r, 4)
                     * (r
                                * (12.288 * r * ipow(r - 1.0, 4) * cos_11_theta
                                   + 17.2032 * ipow(r - 1.0, 5) * cos_11_theta)
                                * std::exp(-tanh_term)
                        + r
                                  * (2.4576 * r * ipow(r - 1.0, 5) * cos_11_theta
                                     + 2.4576 * minus_pow_6 * cos_11_theta)
                                  * (20.0 * ipow(tanh_term, 2) - 20.0) * std::exp(-tanh_term)
                        - 49.5616 * minus_pow_6 * std::exp(-tanh_term) * cos_11_theta
                        + 6.0
                                  * (2.4576 * r * ipow(r - 1.0, 5) * cos_11_theta
                                     + 2.4576 * minus_pow_6 * cos_11_theta)
                                  * std::exp(-tanh_term));
}

//---------------------------------------------------------------------

template <>
double ManufacturedPoissonTest<CurvilinearSolution<CzarnyToCartesian<X, Y, R, Theta>>>::
        solution_at_pole(Coord<R, Theta> const& coord) const
{
    return 0.0;
}

template <>
double ManufacturedPoissonTest<CurvilinearSolution<CzarnyToCartesian<X, Y, R, Theta>>>::
        non_singular_solution(Coord<R, Theta> const& coord) const
{
    const double r = ddc::get<R>(coord);
    const double theta = ddc::get<Theta>(coord);
    const double epsilon = m_coordinate_converter.epsilon();
    const double e = m_coordinate_converter.e();

    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);

    const double tanh_term = std::tanh(20.0 * r - 14.0);
    const double coeff_alpha = std::exp(tanh_term);

    const double xi = 1. / std::sqrt(1. - epsilon * epsilon * 0.25);
    const double tmp1 = std::sqrt(epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0);

    return 0.4096 * ipow(r, 6) * ipow(r - 1.0, 6) * coeff_alpha * std::cos(11.0 * theta)
           - ipow(r, 4)
                     * (4.5056 * r * ipow(r - 1.0, 6)
                                * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                    + e * cos_theta * xi / ((2.0 - tmp1)))
                                           * (e * epsilon * r * sin_theta * cos_theta * xi
                                                      / (ipow(2.0 - tmp1, 2) * tmp1)
                                              + e * sin_theta * xi / ((2.0 - tmp1)))
                                   - sin_theta * cos_theta
                                             / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                * (20.0 * pow(tanh_term, 2.0) - 20.0) * std::exp(-tanh_term)
                                * std::sin(11.0 * theta)
                                / std::sqrt(
                                        (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                                + e * cos_theta * xi / ((2.0 - tmp1)))
                                                       * (e * epsilon * r * sin_theta * cos_theta
                                                                  * xi
                                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                                          + e * sin_theta * xi / ((2.0 - tmp1)))
                                               - sin_theta * cos_theta
                                                         / (epsilon
                                                                    * (epsilon
                                                                       + 2.0 * r * cos_theta)
                                                            + 1.0)),
                                              2.0))
                                        + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                                + e * cos_theta * xi / ((2.0 - tmp1))),
                                               2.0)
                                           + ipow(sin_theta, 2)
                                                     / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                        + 1.0))
                                                  * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                  * xi
                                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                                          + e * sin_theta * xi / ((2.0 - tmp1))),
                                                         2.0)
                                                     + ipow(cos_theta, 2)
                                                               / (epsilon
                                                                          * (epsilon
                                                                             + 2.0 * r * cos_theta)
                                                                  + 1.0)))
                        + 4.5056 * r * ipow(r - 1.0, 6)
                                  * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + e * cos_theta * xi / ((2.0 - tmp1)))
                                             * (e * epsilon * r * sin_theta * cos_theta * xi
                                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                                + e * sin_theta * xi / ((2.0 - tmp1)))
                                     - sin_theta * cos_theta
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * (1.0 / 2.0
                                             * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1)))
                                                        * (e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                - sin_theta * cos_theta
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0))
                                             * (4.0 * epsilon * sin_theta * ipow(cos_theta, 2)
                                                        / pow((epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0),
                                                              2.0)
                                                + 2.0
                                                          * ((-e) * epsilon * r * ipow(sin_theta, 2)
                                                                     * xi
                                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                                             + e * cos_theta * xi / ((2.0 - tmp1)))
                                                          * ((-e) * (epsilon * epsilon) * r
                                                                     * sin_theta
                                                                     * ipow(cos_theta, 2) * xi
                                                                     / (ipow(2.0 - tmp1, 2)
                                                                        * pow((epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0),
                                                                              (3.0 / 2.0)))
                                                             + 2.0 * e * (epsilon * epsilon) * r
                                                                       * sin_theta
                                                                       * ipow(cos_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 3)
                                                                          * (epsilon
                                                                                     * (epsilon
                                                                                        + 2.0 * r * cos_theta)
                                                                             + 1.0))
                                                             + 2.0 * e * epsilon * sin_theta
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1))
                                                + 2.0
                                                          * (e * epsilon * r * sin_theta * cos_theta
                                                                     * xi
                                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                                             + e * sin_theta * xi / ((2.0 - tmp1)))
                                                          * (e * (epsilon * epsilon) * r
                                                                     * ipow(sin_theta, 2)
                                                                     * cos_theta * xi
                                                                     / (ipow(2.0 - tmp1, 2)
                                                                        * pow((epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0),
                                                                              (3.0 / 2.0)))
                                                             - 2.0 * e * (epsilon * epsilon) * r
                                                                       * ipow(sin_theta, 2)
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 3)
                                                                          * (epsilon
                                                                                     * (epsilon
                                                                                        + 2.0 * r * cos_theta)
                                                                             + 1.0))
                                                             - e * epsilon * ipow(sin_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)
                                                             + e * epsilon * ipow(cos_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)))
                                     - 1.0 / 2.0
                                               * ((-2.0) * epsilon * ipow(cos_theta, 3)
                                                          / pow((epsilon
                                                                         * (epsilon
                                                                            + 2.0 * r * cos_theta)
                                                                 + 1.0),
                                                                2.0)
                                                  + (e * epsilon * r * sin_theta * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + e * sin_theta * xi / ((2.0 - tmp1)))
                                                            * ((-2.0) * e * (epsilon * epsilon) * r
                                                                       * sin_theta
                                                                       * ipow(cos_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * pow((epsilon
                                                                                         * (epsilon
                                                                                            + 2.0 * r * cos_theta)
                                                                                 + 1.0),
                                                                                (3.0 / 2.0)))
                                                               + 4.0 * e * (epsilon * epsilon) * r
                                                                         * sin_theta
                                                                         * ipow(cos_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 3)
                                                                            * (epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0))
                                                               + 4.0 * e * epsilon * sin_theta
                                                                         * cos_theta * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)))
                                               * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                                      2.0)
                                                  + ipow(sin_theta, 2)
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0))
                                     - 1.0 / 2.0
                                               * ((-2.0) * epsilon * ipow(sin_theta, 2) * cos_theta
                                                          / pow((epsilon
                                                                         * (epsilon
                                                                            + 2.0 * r * cos_theta)
                                                                 + 1.0),
                                                                2.0)
                                                  + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + e * cos_theta * xi / ((2.0 - tmp1)))
                                                            * (2.0 * e * (epsilon * epsilon) * r
                                                                       * ipow(sin_theta, 2)
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * pow((epsilon
                                                                                         * (epsilon
                                                                                            + 2.0 * r * cos_theta)
                                                                                 + 1.0),
                                                                                (3.0 / 2.0)))
                                                               - 4.0 * e * (epsilon * epsilon) * r
                                                                         * ipow(sin_theta, 2)
                                                                         * cos_theta * xi
                                                                         / (ipow(2.0 - tmp1, 3)
                                                                            * (epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0))
                                                               - 2.0 * e * epsilon
                                                                         * ipow(sin_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)
                                                               + 2.0 * e * epsilon
                                                                         * ipow(cos_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)))
                                               * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                       + e * sin_theta * xi / ((2.0 - tmp1))),
                                                      2.0)
                                                  + ipow(cos_theta, 2)
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0)))
                                  * std::exp(-tanh_term) * std::sin(11.0 * theta)
                                  / pow(((-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1)))
                                                        * (e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                - sin_theta * cos_theta
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0)),
                                               2.0))
                                         + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1))),
                                                2.0)
                                            + ipow(sin_theta, 2)
                                                      / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                         + 1.0))
                                                   * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1))),
                                                          2.0)
                                                      + ipow(cos_theta, 2)
                                                                / (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))),
                                        (3.0 / 2.0))
                        + 4.5056 * r * ipow(r - 1.0, 6)
                                  * (2.0 * epsilon * sin_theta * ipow(cos_theta, 2)
                                             / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0),
                                                   2.0)
                                     + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * ((-e) * (epsilon * epsilon) * r * sin_theta
                                                          * ipow(cos_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2)
                                                             * pow((epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0),
                                                                   (3.0 / 2.0)))
                                                  + 2.0 * e * (epsilon * epsilon) * r * sin_theta
                                                            * ipow(cos_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 3)
                                                               * (epsilon
                                                                          * (epsilon
                                                                             + 2.0 * r * cos_theta)
                                                                  + 1.0))
                                                  + 2.0 * e * epsilon * sin_theta * cos_theta * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1))
                                     + (e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * sin_theta * xi / ((2.0 - tmp1)))
                                               * (e * (epsilon * epsilon) * r * ipow(sin_theta, 2)
                                                          * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2)
                                                             * pow((epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0),
                                                                   (3.0 / 2.0)))
                                                  - 2.0 * e * (epsilon * epsilon) * r
                                                            * ipow(sin_theta, 2) * cos_theta * xi
                                                            / (ipow(2.0 - tmp1, 3)
                                                               * (epsilon
                                                                          * (epsilon
                                                                             + 2.0 * r * cos_theta)
                                                                  + 1.0))
                                                  - e * epsilon * ipow(sin_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * epsilon * ipow(cos_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1)))
                                  * std::exp(-tanh_term) * std::sin(11.0 * theta)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        + 27.0336 * r * ipow(r - 1.0, 5)
                                  * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + e * cos_theta * xi / ((2.0 - tmp1)))
                                             * (e * epsilon * r * sin_theta * cos_theta * xi
                                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                                + e * sin_theta * xi / ((2.0 - tmp1)))
                                     - sin_theta * cos_theta
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * std::exp(-tanh_term) * std::sin(11.0 * theta)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        + r
                                  * (12.288 * r * ipow(r - 1.0, 4) * std::cos(11.0 * theta)
                                     + 17.2032 * ipow(r - 1.0, 5) * std::cos(11.0 * theta))
                                  * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                          + e * cos_theta * xi / ((2.0 - tmp1))),
                                         2.0)
                                     + ipow(sin_theta, 2)
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * std::exp(-tanh_term)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta
                                                            / (std::sqrt(
                                                                       1.0
                                                                       - 1.0 / 4.0
                                                                                 * (epsilon
                                                                                    * epsilon))
                                                               * (2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2)
                                                          / (std::sqrt(
                                                                     1.0
                                                                     - 1.0 / 4.0
                                                                               * (epsilon
                                                                                  * epsilon))
                                                             * ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta
                                                            / (std::sqrt(
                                                                       1.0
                                                                       - 1.0 / 4.0
                                                                                 * (epsilon
                                                                                    * epsilon))
                                                               * (2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    / (std::sqrt(
                                                                               1.0
                                                                               - 1.0 / 4.0
                                                                                         * (epsilon
                                                                                            * epsilon))
                                                                       * ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        + r
                                  * (2.4576 * r * ipow(r - 1.0, 5) * std::cos(11.0 * theta)
                                     + 2.4576 * ipow(r - 1.0, 6) * std::cos(11.0 * theta))
                                  * ((-2.0) * epsilon * ipow(sin_theta, 2) * cos_theta
                                             / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0),
                                                   2.0)
                                     + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (2.0 * e * (epsilon * epsilon) * r
                                                          * ipow(sin_theta, 2) * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2)
                                                             * pow((epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0),
                                                                   (3.0 / 2.0)))
                                                  - 4.0 * e * (epsilon * epsilon) * r
                                                            * ipow(sin_theta, 2) * cos_theta * xi
                                                            / (ipow(2.0 - tmp1, 3)
                                                               * (epsilon
                                                                          * (epsilon
                                                                             + 2.0 * r * cos_theta)
                                                                  + 1.0))
                                                  - 2.0 * e * epsilon * ipow(sin_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + 2.0 * e * epsilon * ipow(cos_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1)))
                                  * std::exp(-tanh_term)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        + r
                                  * (2.4576 * r * ipow(r - 1.0, 5) * std::cos(11.0 * theta)
                                     + 2.4576 * ipow(r - 1.0, 6) * std::cos(11.0 * theta))
                                  * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                          + e * cos_theta * xi / ((2.0 - tmp1))),
                                         2.0)
                                     + ipow(sin_theta, 2)
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * (20.0 * pow(tanh_term, 2.0) - 20.0) * std::exp(-tanh_term)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2)
                                                          / (std::sqrt(
                                                                     1.0
                                                                     - 1.0 / 4.0
                                                                               * (epsilon
                                                                                  * epsilon))
                                                             * ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        + r
                                  * (2.4576 * r * ipow(r - 1.0, 5) * std::cos(11.0 * theta)
                                     + 2.4576 * ipow(r - 1.0, 6) * std::cos(11.0 * theta))
                                  * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                          + e * cos_theta * xi / ((2.0 - tmp1))),
                                         2.0)
                                     + ipow(sin_theta, 2)
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * (1.0 / 2.0
                                             * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1)))
                                                        * (e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                - sin_theta * cos_theta
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0))
                                             * (4.0 * epsilon * sin_theta * ipow(cos_theta, 2)
                                                        / pow((epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0),
                                                              2.0)
                                                + 2.0
                                                          * ((-e) * epsilon * r * ipow(sin_theta, 2)
                                                                     * xi
                                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                                             + e * cos_theta * xi / ((2.0 - tmp1)))
                                                          * ((-e) * (epsilon * epsilon) * r
                                                                     * sin_theta
                                                                     * ipow(cos_theta, 2) * xi
                                                                     / (ipow(2.0 - tmp1, 2)
                                                                        * pow((epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0),
                                                                              (3.0 / 2.0)))
                                                             + 2.0 * e * (epsilon * epsilon) * r
                                                                       * sin_theta
                                                                       * ipow(cos_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 3)
                                                                          * (epsilon
                                                                                     * (epsilon
                                                                                        + 2.0 * r * cos_theta)
                                                                             + 1.0))
                                                             + 2.0 * e * epsilon * sin_theta
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1))
                                                + 2.0
                                                          * (e * epsilon * r * sin_theta * cos_theta
                                                                     * xi
                                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                                             + e * sin_theta * xi / ((2.0 - tmp1)))
                                                          * (e * (epsilon * epsilon) * r
                                                                     * ipow(sin_theta, 2)
                                                                     * cos_theta * xi
                                                                     / (ipow(2.0 - tmp1, 2)
                                                                        * pow((epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0),
                                                                              (3.0 / 2.0)))
                                                             - 2.0 * e * (epsilon * epsilon) * r
                                                                       * ipow(sin_theta, 2)
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 3)
                                                                          * (epsilon
                                                                                     * (epsilon
                                                                                        + 2.0 * r * cos_theta)
                                                                             + 1.0))
                                                             - e * epsilon * ipow(sin_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)
                                                             + e * epsilon * ipow(cos_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)))
                                     - 1.0 / 2.0
                                               * ((-2.0) * epsilon * ipow(cos_theta, 3)
                                                          / pow((epsilon
                                                                         * (epsilon
                                                                            + 2.0 * r * cos_theta)
                                                                 + 1.0),
                                                                2.0)
                                                  + (e * epsilon * r * sin_theta * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + e * sin_theta * xi / ((2.0 - tmp1)))
                                                            * ((-2.0) * e * (epsilon * epsilon) * r
                                                                       * sin_theta
                                                                       * ipow(cos_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * pow((epsilon
                                                                                         * (epsilon
                                                                                            + 2.0 * r * cos_theta)
                                                                                 + 1.0),
                                                                                (3.0 / 2.0)))
                                                               + 4.0 * e * (epsilon * epsilon) * r
                                                                         * sin_theta
                                                                         * ipow(cos_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 3)
                                                                            * (epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0))
                                                               + 4.0 * e * epsilon * sin_theta
                                                                         * cos_theta * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)))
                                               * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                                      2.0)
                                                  + ipow(sin_theta, 2)
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0))
                                     - 1.0 / 2.0
                                               * ((-2.0) * epsilon * ipow(sin_theta, 2) * cos_theta
                                                          / pow((epsilon
                                                                         * (epsilon
                                                                            + 2.0 * r * cos_theta)
                                                                 + 1.0),
                                                                2.0)
                                                  + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + e * cos_theta * xi / ((2.0 - tmp1)))
                                                            * (2.0 * e * (epsilon * epsilon) * r
                                                                       * ipow(sin_theta, 2)
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * pow((epsilon
                                                                                         * (epsilon
                                                                                            + 2.0 * r * cos_theta)
                                                                                 + 1.0),
                                                                                (3.0 / 2.0)))
                                                               - 4.0 * e * (epsilon * epsilon) * r
                                                                         * ipow(sin_theta, 2)
                                                                         * cos_theta * xi
                                                                         / (ipow(2.0 - tmp1, 3)
                                                                            * (epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0))
                                                               - 2.0 * e * epsilon
                                                                         * ipow(sin_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)
                                                               + 2.0 * e * epsilon
                                                                         * ipow(cos_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)))
                                               * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                       + e * sin_theta * xi / ((2.0 - tmp1))),
                                                      2.0)
                                                  + ipow(cos_theta, 2)
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0)))
                                  * std::exp(-tanh_term)
                                  / pow(((-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1)))
                                                        * (e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                - sin_theta * cos_theta
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0)),
                                               2.0))
                                         + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1))),
                                                2.0)
                                            + ipow(sin_theta, 2)
                                                      / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                         + 1.0))
                                                   * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1))),
                                                          2.0)
                                                      + ipow(cos_theta, 2)
                                                                / (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))),
                                        (3.0 / 2.0))
                        + 27.0336 * ipow(r - 1.0, 6)
                                  * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + e * cos_theta * xi / ((2.0 - tmp1)))
                                             * (e * epsilon * r * sin_theta * cos_theta * xi
                                                        / (ipow(2.0 - tmp1, 2) * tmp1)
                                                + e * sin_theta * xi / ((2.0 - tmp1)))
                                     - sin_theta * cos_theta
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * std::exp(-tanh_term) * std::sin(11.0 * theta)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        - 49.5616 * ipow(r - 1.0, 6)
                                  * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                          + e * sin_theta * xi / ((2.0 - tmp1))),
                                         2.0)
                                     + ipow(cos_theta, 2)
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * std::exp(-tanh_term) * std::cos(11.0 * theta)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        - 4.5056 * ipow(r - 1.0, 6)
                                  * (2.0 * epsilon * r * sin_theta * ipow(cos_theta, 2)
                                             / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0),
                                                   2.0)
                                     + (e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * sin_theta * xi / ((2.0 - tmp1)))
                                               * (2.0 * e * (epsilon * epsilon) * (r * r)
                                                          * ipow(sin_theta, 2) * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2)
                                                             * pow((epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0),
                                                                   (3.0 / 2.0)))
                                                  - 4.0 * e * (epsilon * epsilon) * (r * r)
                                                            * ipow(sin_theta, 2) * cos_theta * xi
                                                            / (ipow(2.0 - tmp1, 3)
                                                               * (epsilon
                                                                          * (epsilon
                                                                             + 2.0 * r * cos_theta)
                                                                  + 1.0))
                                                  - 4.0 * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + 2.0 * e * epsilon * r * ipow(cos_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + 2.0 * e * cos_theta * xi / ((2.0 - tmp1)))
                                     - 2.0 * sin_theta * cos_theta
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * std::exp(-tanh_term) * std::sin(11.0 * theta)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta
                                                            / (std::sqrt(
                                                                       1.0
                                                                       - 1.0 / 4.0
                                                                                 * (epsilon
                                                                                    * epsilon))
                                                               * (2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        - 4.5056 * ipow(r - 1.0, 6)
                                  * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                          + e * sin_theta * xi / ((2.0 - tmp1))),
                                         2.0)
                                     + ipow(cos_theta, 2)
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * (1.0 / 2.0
                                             * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1)))
                                                        * (e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                - sin_theta * cos_theta
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0))
                                             * ((-4.0) * epsilon * r * ipow(sin_theta, 2)
                                                        * cos_theta
                                                        / pow((epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0),
                                                              2.0)
                                                + 2.0
                                                          * ((-e) * epsilon * r * ipow(sin_theta, 2)
                                                                     * xi
                                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                                             + e * cos_theta * xi / ((2.0 - tmp1)))
                                                          * (e * (epsilon * epsilon) * (r * r)
                                                                     * ipow(sin_theta, 2)
                                                                     * cos_theta * xi
                                                                     / (ipow(2.0 - tmp1, 2)
                                                                        * pow((epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0),
                                                                              (3.0 / 2.0)))
                                                             - 2.0 * e * (epsilon * epsilon)
                                                                       * (r * r)
                                                                       * ipow(sin_theta, 2)
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 3)
                                                                          * (epsilon
                                                                                     * (epsilon
                                                                                        + 2.0 * r * cos_theta)
                                                                             + 1.0))
                                                             - 2.0 * e * epsilon * r
                                                                       * ipow(sin_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)
                                                             + e * epsilon * r * ipow(cos_theta, 2)
                                                                       * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)
                                                             + e * cos_theta * xi / ((2.0 - tmp1)))
                                                + 2.0
                                                          * (e * epsilon * r * sin_theta * cos_theta
                                                                     * xi
                                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                                             + e * sin_theta * xi / ((2.0 - tmp1)))
                                                          * ((-e) * (epsilon * epsilon) * (r * r)
                                                                     * ipow(sin_theta, 3) * xi
                                                                     / (ipow(2.0 - tmp1, 2)
                                                                        * pow((epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0),
                                                                              (3.0 / 2.0)))
                                                             + 2.0 * e * (epsilon * epsilon)
                                                                       * (r * r)
                                                                       * ipow(sin_theta, 3) * xi
                                                                       / (ipow(2.0 - tmp1, 3)
                                                                          * (epsilon
                                                                                     * (epsilon
                                                                                        + 2.0 * r * cos_theta)
                                                                             + 1.0))
                                                             - 3.0 * e * epsilon * r * sin_theta
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)
                                                             - e * sin_theta * xi / ((2.0 - tmp1)))
                                                + 2.0 * ipow(sin_theta, 2)
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0)
                                                - 2.0 * ipow(cos_theta, 2)
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0))
                                     - 1.0 / 2.0
                                               * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                                      2.0)
                                                  + ipow(sin_theta, 2)
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0))
                                               * (2.0 * epsilon * r * sin_theta * ipow(cos_theta, 2)
                                                          / pow((epsilon
                                                                         * (epsilon
                                                                            + 2.0 * r * cos_theta)
                                                                 + 1.0),
                                                                2.0)
                                                  + (e * epsilon * r * sin_theta * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + e * sin_theta * xi / ((2.0 - tmp1)))
                                                            * (2.0 * e * (epsilon * epsilon)
                                                                       * (r * r)
                                                                       * ipow(sin_theta, 2)
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * pow((epsilon
                                                                                         * (epsilon
                                                                                            + 2.0 * r * cos_theta)
                                                                                 + 1.0),
                                                                                (3.0 / 2.0)))
                                                               - 4.0 * e * (epsilon * epsilon)
                                                                         * (r * r)
                                                                         * ipow(sin_theta, 2)
                                                                         * cos_theta * xi
                                                                         / (ipow(2.0 - tmp1, 3)
                                                                            * (epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0))
                                                               - 4.0 * e * epsilon * r
                                                                         * ipow(sin_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)
                                                               + 2.0 * e * epsilon * r
                                                                         * ipow(cos_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)
                                                               + 2.0 * e * cos_theta * xi
                                                                         / ((2.0 - tmp1)))
                                                  - 2.0 * sin_theta * cos_theta
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0))
                                     - 1.0 / 2.0
                                               * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                       + e * sin_theta * xi / ((2.0 - tmp1))),
                                                      2.0)
                                                  + ipow(cos_theta, 2)
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0))
                                               * (2.0 * epsilon * r * ipow(sin_theta, 3)
                                                          / pow((epsilon
                                                                         * (epsilon
                                                                            + 2.0 * r * cos_theta)
                                                                 + 1.0),
                                                                2.0)
                                                  + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + e * cos_theta * xi / ((2.0 - tmp1)))
                                                            * ((-2.0) * e * (epsilon * epsilon)
                                                                       * (r * r)
                                                                       * ipow(sin_theta, 3) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * pow((epsilon
                                                                                         * (epsilon
                                                                                            + 2.0 * r * cos_theta)
                                                                                 + 1.0),
                                                                                (3.0 / 2.0)))
                                                               + 4.0 * e * (epsilon * epsilon)
                                                                         * (r * r)
                                                                         * ipow(sin_theta, 3) * xi
                                                                         / (ipow(2.0 - tmp1, 3)
                                                                            * (epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0))
                                                               - 6.0 * e * epsilon * r * sin_theta
                                                                         * cos_theta * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)
                                                               - 2.0 * e * sin_theta * xi
                                                                         / ((2.0 - tmp1)))
                                                  + 2.0 * sin_theta * cos_theta
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0)))
                                  * std::exp(-tanh_term) * std::sin(11.0 * theta)
                                  / pow(((-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1)))
                                                        * (e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                - sin_theta * cos_theta
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0)),
                                               2.0))
                                         + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1))),
                                                2.0)
                                            + ipow(sin_theta, 2)
                                                      / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                         + 1.0))
                                                   * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1))),
                                                          2.0)
                                                      + ipow(cos_theta, 2)
                                                                / (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))),
                                        (3.0 / 2.0))
                        - (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + e * cos_theta * xi / ((2.0 - tmp1)))
                                   * (e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + e * sin_theta * xi / ((2.0 - tmp1)))
                           - sin_theta * cos_theta
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * ((-27.0336) * r * ipow(r - 1.0, 5) * std::sin(11.0 * theta)
                                     - 27.0336 * ipow(r - 1.0, 6) * std::sin(11.0 * theta))
                                  * std::exp(-tanh_term)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        - (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                            + e * cos_theta * xi / ((2.0 - tmp1)))
                                   * (e * epsilon * r * sin_theta * cos_theta * xi
                                              / (ipow(2.0 - tmp1, 2) * tmp1)
                                      + e * sin_theta * xi / ((2.0 - tmp1)))
                           - sin_theta * cos_theta
                                     / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * (2.4576 * r * ipow(r - 1.0, 5) * std::cos(11.0 * theta)
                                     + 2.4576 * ipow(r - 1.0, 6) * std::cos(11.0 * theta))
                                  * (1.0 / 2.0
                                             * (((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1)))
                                                        * (e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                - sin_theta * cos_theta
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0))
                                             * ((-4.0) * epsilon * r * ipow(sin_theta, 2)
                                                        * cos_theta
                                                        / pow((epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0),
                                                              2.0)
                                                + 2.0
                                                          * ((-e) * epsilon * r * ipow(sin_theta, 2)
                                                                     * xi
                                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                                             + e * cos_theta * xi / ((2.0 - tmp1)))
                                                          * (e * (epsilon * epsilon) * (r * r)
                                                                     * ipow(sin_theta, 2)
                                                                     * cos_theta * xi
                                                                     / (ipow(2.0 - tmp1, 2)
                                                                        * pow((epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0),
                                                                              (3.0 / 2.0)))
                                                             - 2.0 * e * (epsilon * epsilon)
                                                                       * (r * r)
                                                                       * ipow(sin_theta, 2)
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 3)
                                                                          * (epsilon
                                                                                     * (epsilon
                                                                                        + 2.0 * r * cos_theta)
                                                                             + 1.0))
                                                             - 2.0 * e * epsilon * r
                                                                       * ipow(sin_theta, 2) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)
                                                             + e * epsilon * r * ipow(cos_theta, 2)
                                                                       * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)
                                                             + e * cos_theta * xi / ((2.0 - tmp1)))
                                                + 2.0
                                                          * (e * epsilon * r * sin_theta * cos_theta
                                                                     * xi
                                                                     / (ipow(2.0 - tmp1, 2) * tmp1)
                                                             + e * sin_theta * xi / ((2.0 - tmp1)))
                                                          * ((-e) * (epsilon * epsilon) * (r * r)
                                                                     * ipow(sin_theta, 3) * xi
                                                                     / (ipow(2.0 - tmp1, 2)
                                                                        * pow((epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0),
                                                                              (3.0 / 2.0)))
                                                             + 2.0 * e * (epsilon * epsilon)
                                                                       * (r * r)
                                                                       * ipow(sin_theta, 3) * xi
                                                                       / (ipow(2.0 - tmp1, 3)
                                                                          * (epsilon
                                                                                     * (epsilon
                                                                                        + 2.0 * r * cos_theta)
                                                                             + 1.0))
                                                             - 3.0 * e * epsilon * r * sin_theta
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * tmp1)
                                                             - e * sin_theta * xi / ((2.0 - tmp1)))
                                                + 2.0 * ipow(sin_theta, 2)
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0)
                                                - 2.0 * ipow(cos_theta, 2)
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0))
                                     - 1.0 / 2.0
                                               * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                       + e * cos_theta * xi / ((2.0 - tmp1))),
                                                      2.0)
                                                  + ipow(sin_theta, 2)
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0))
                                               * (2.0 * epsilon * r * sin_theta * ipow(cos_theta, 2)
                                                          / pow((epsilon
                                                                         * (epsilon
                                                                            + 2.0 * r * cos_theta)
                                                                 + 1.0),
                                                                2.0)
                                                  + (e * epsilon * r * sin_theta * cos_theta * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + e * sin_theta * xi / ((2.0 - tmp1)))
                                                            * (2.0 * e * (epsilon * epsilon)
                                                                       * (r * r)
                                                                       * ipow(sin_theta, 2)
                                                                       * cos_theta * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * pow((epsilon
                                                                                         * (epsilon
                                                                                            + 2.0 * r * cos_theta)
                                                                                 + 1.0),
                                                                                (3.0 / 2.0)))
                                                               - 4.0 * e * (epsilon * epsilon)
                                                                         * (r * r)
                                                                         * ipow(sin_theta, 2)
                                                                         * cos_theta * xi
                                                                         / (ipow(2.0 - tmp1, 3)
                                                                            * (epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0))
                                                               - 4.0 * e * epsilon * r
                                                                         * ipow(sin_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)
                                                               + 2.0 * e * epsilon * r
                                                                         * ipow(cos_theta, 2) * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)
                                                               + 2.0 * e * cos_theta * xi
                                                                         / ((2.0 - tmp1)))
                                                  - 2.0 * sin_theta * cos_theta
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0))
                                     - 1.0 / 2.0
                                               * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                               / (ipow(2.0 - tmp1, 2) * tmp1)
                                                       + e * sin_theta * xi / ((2.0 - tmp1))),
                                                      2.0)
                                                  + ipow(cos_theta, 2)
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0))
                                               * (2.0 * epsilon * r * ipow(sin_theta, 3)
                                                          / pow((epsilon
                                                                         * (epsilon
                                                                            + 2.0 * r * cos_theta)
                                                                 + 1.0),
                                                                2.0)
                                                  + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                                     + e * cos_theta * xi / ((2.0 - tmp1)))
                                                            * ((-2.0) * e * (epsilon * epsilon)
                                                                       * (r * r)
                                                                       * ipow(sin_theta, 3) * xi
                                                                       / (ipow(2.0 - tmp1, 2)
                                                                          * pow((epsilon
                                                                                         * (epsilon
                                                                                            + 2.0 * r * cos_theta)
                                                                                 + 1.0),
                                                                                (3.0 / 2.0)))
                                                               + 4.0 * e * (epsilon * epsilon)
                                                                         * (r * r)
                                                                         * ipow(sin_theta, 3) * xi
                                                                         / (ipow(2.0 - tmp1, 3)
                                                                            * (epsilon
                                                                                       * (epsilon
                                                                                          + 2.0 * r * cos_theta)
                                                                               + 1.0))
                                                               - 6.0 * e * epsilon * r * sin_theta
                                                                         * cos_theta * xi
                                                                         / (ipow(2.0 - tmp1, 2)
                                                                            * tmp1)
                                                               - 2.0 * e * sin_theta * xi
                                                                         / ((2.0 - tmp1)))
                                                  + 2.0 * sin_theta * cos_theta
                                                            / (epsilon
                                                                       * (epsilon
                                                                          + 2.0 * r * cos_theta)
                                                               + 1.0)))
                                  * std::exp(-tanh_term)
                                  / pow(((-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1)))
                                                        * (e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1)))
                                                - sin_theta * cos_theta
                                                          / (epsilon
                                                                     * (epsilon
                                                                        + 2.0 * r * cos_theta)
                                                             + 1.0)),
                                               2.0))
                                         + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                         / (ipow(2.0 - tmp1, 2) * tmp1)
                                                 + e * cos_theta * xi / ((2.0 - tmp1))),
                                                2.0)
                                            + ipow(sin_theta, 2)
                                                      / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                         + 1.0))
                                                   * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                   * xi
                                                                   / (ipow(2.0 - tmp1, 2) * tmp1)
                                                           + e * sin_theta * xi / ((2.0 - tmp1))),
                                                          2.0)
                                                      + ipow(cos_theta, 2)
                                                                / (epsilon
                                                                           * (epsilon
                                                                              + 2.0 * r * cos_theta)
                                                                   + 1.0))),
                                        (3.0 / 2.0))
                        + 6.0
                                  * (2.4576 * r * ipow(r - 1.0, 5) * std::cos(11.0 * theta)
                                     + 2.4576 * ipow(r - 1.0, 6) * std::cos(11.0 * theta))
                                  * (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                  / (ipow(2.0 - tmp1, 2) * tmp1)
                                          + e * cos_theta * xi / ((2.0 - tmp1))),
                                         2.0)
                                     + ipow(sin_theta, 2)
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * std::exp(-tanh_term)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0)))
                        - (2.4576 * r * ipow(r - 1.0, 5) * std::cos(11.0 * theta)
                           + 2.4576 * ipow(r - 1.0, 6) * std::cos(11.0 * theta))
                                  * ((-2.0) * epsilon * r * ipow(sin_theta, 2) * cos_theta
                                             / pow((epsilon * (epsilon + 2.0 * r * cos_theta)
                                                    + 1.0),
                                                   2.0)
                                     + ((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * cos_theta * xi / ((2.0 - tmp1)))
                                               * (e * (epsilon * epsilon) * (r * r)
                                                          * ipow(sin_theta, 2) * cos_theta * xi
                                                          / (ipow(2.0 - tmp1, 2)
                                                             * pow((epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0),
                                                                   (3.0 / 2.0)))
                                                  - 2.0 * e * (epsilon * epsilon) * (r * r)
                                                            * ipow(sin_theta, 2) * cos_theta * xi
                                                            / (ipow(2.0 - tmp1, 3)
                                                               * (epsilon
                                                                          * (epsilon
                                                                             + 2.0 * r * cos_theta)
                                                                  + 1.0))
                                                  - 2.0 * e * epsilon * r * ipow(sin_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * epsilon * r * ipow(cos_theta, 2) * xi
                                                            / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                     + (e * epsilon * r * sin_theta * cos_theta * xi
                                                / (ipow(2.0 - tmp1, 2) * tmp1)
                                        + e * sin_theta * xi / ((2.0 - tmp1)))
                                               * ((-e) * (epsilon * epsilon) * (r * r)
                                                          * ipow(sin_theta, 3) * xi
                                                          / (ipow(2.0 - tmp1, 2)
                                                             * pow((epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0),
                                                                   (3.0 / 2.0)))
                                                  + 2.0 * e * (epsilon * epsilon) * (r * r)
                                                            * ipow(sin_theta, 3) * xi
                                                            / (ipow(2.0 - tmp1, 3)
                                                               * (epsilon
                                                                          * (epsilon
                                                                             + 2.0 * r * cos_theta)
                                                                  + 1.0))
                                                  - 3.0 * e * epsilon * r * sin_theta * cos_theta
                                                            * xi / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  - e * sin_theta * xi / ((2.0 - tmp1)))
                                     + ipow(sin_theta, 2)
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)
                                     - ipow(cos_theta, 2)
                                               / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                  * std::exp(-tanh_term)
                                  / std::sqrt(
                                          (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1)))
                                                         * (e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1)))
                                                 - sin_theta * cos_theta
                                                           / (epsilon
                                                                      * (epsilon
                                                                         + 2.0 * r * cos_theta)
                                                              + 1.0)),
                                                2.0))
                                          + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                                          / (ipow(2.0 - tmp1, 2) * tmp1)
                                                  + e * cos_theta * xi / ((2.0 - tmp1))),
                                                 2.0)
                                             + ipow(sin_theta, 2)
                                                       / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                          + 1.0))
                                                    * (pow((e * epsilon * r * sin_theta * cos_theta
                                                                    * xi
                                                                    / (ipow(2.0 - tmp1, 2) * tmp1)
                                                            + e * sin_theta * xi / ((2.0 - tmp1))),
                                                           2.0)
                                                       + ipow(cos_theta, 2)
                                                                 / (epsilon
                                                                            * (epsilon
                                                                               + 2.0 * r * cos_theta)
                                                                    + 1.0))))
                     / std::sqrt(
                             (-pow((((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                     + e * cos_theta * xi / ((2.0 - tmp1)))
                                            * (e * epsilon * r * sin_theta * cos_theta * xi
                                                       / (ipow(2.0 - tmp1, 2) * tmp1)
                                               + e * sin_theta * xi / ((2.0 - tmp1)))
                                    - sin_theta * cos_theta
                                              / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0)),
                                   2.0))
                             + (pow(((-e) * epsilon * r * ipow(sin_theta, 2) * xi
                                             / (ipow(2.0 - tmp1, 2) * tmp1)
                                     + e * cos_theta * xi / ((2.0 - tmp1))),
                                    2.0)
                                + ipow(sin_theta, 2)
                                          / (epsilon * (epsilon + 2.0 * r * cos_theta) + 1.0))
                                       * (pow((e * epsilon * r * sin_theta * cos_theta * xi
                                                       / (ipow(2.0 - tmp1, 2) * tmp1)
                                               + e * sin_theta * xi / ((2.0 - tmp1))),
                                              2.0)
                                          + ipow(cos_theta, 2)
                                                    / (epsilon * (epsilon + 2.0 * r * cos_theta)
                                                       + 1.0)));
}

template class ManufacturedPoissonTest<CurvilinearSolution<CircularToCartesian<X, Y, R, Theta>>>;
template class ManufacturedPoissonTest<CurvilinearSolution<CzarnyToCartesian<X, Y, R, Theta>>>;
template class ManufacturedPoissonTest<CartesianSolution<CircularToCartesian<X, Y, R, Theta>>>;
template class ManufacturedPoissonTest<CartesianSolution<CzarnyToCartesian<X, Y, R, Theta>>>;
