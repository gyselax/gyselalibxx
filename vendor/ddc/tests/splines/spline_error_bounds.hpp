#pragma once

#include <algorithm>
#include <array>

#include <ddc/kernels/splines.hpp>

template <class Evaluator>
class SplineErrorBounds
{
private:
    static constexpr std::array<double, 10> tihomirov_error_bound_array = std::array<double, 10>(
            {1.0 / 2.0,
             1.0 / 8.0,
             1.0 / 24.0,
             5.0 / 384.0,
             1.0 / 240.0,
             61.0 / 46080.0,
             17.0 / 40320.0,
             277.0 / 2064384.0,
             31.0 / 725760.0,
             50521.0 / 3715891200.0});

private:
    Evaluator const& m_evaluator;

private:
    /*******************************************************************************
     * Error bound in max norm for spline interpolation of periodic functions from:
     *
     * V M Tihomirov 1969 Math. USSR Sb. 9 275
     * https://doi.org/10.1070/SM1969v009n02ABEH002052 (page 286, bottom)
     *
     * Yu. S. Volkov and Yu. N. Subbotin
     * https://doi.org/10.1134/S0081543815020236 (equation 14)
     *
     * Also applicable to first derivative by passing deg-1 instead of deg
     * Volkov & Subbotin 2015, eq. 15
     *******************************************************************************/
    static double tihomirov_error_bound(double cell_width, int degree, double max_norm)
    {
        degree = std::min(degree, 9);
        return tihomirov_error_bound_array[degree] * ddc::detail::ipow(cell_width, degree + 1)
               * max_norm;
    }

public:
    SplineErrorBounds(Evaluator const& evaluator) : m_evaluator(evaluator) {}

    double error_bound(double cell_width, int degree)
    {
        return tihomirov_error_bound(cell_width, degree, m_evaluator.max_norm(degree + 1));
    }
    double error_bound(double cell_width1, double cell_width2, int degree1, int degree2)
    {
        double norm1 = m_evaluator.max_norm(degree1 + 1, 0);
        double norm2 = m_evaluator.max_norm(0, degree2 + 1);
        return tihomirov_error_bound(cell_width1, degree1, norm1)
               + tihomirov_error_bound(cell_width2, degree2, norm2);
    }
    double error_bound_on_deriv(double cell_width, int degree)
    {
        return tihomirov_error_bound(cell_width, degree - 1, m_evaluator.max_norm(degree + 1));
    }
    double error_bound_on_deriv_1(double cell_width1, double cell_width2, int degree1, int degree2)
    {
        double norm1 = m_evaluator.max_norm(degree1 + 1, 0);
        double norm2 = m_evaluator.max_norm(0, degree2 + 1);
        return tihomirov_error_bound(cell_width1, degree1 - 1, norm1)
               + tihomirov_error_bound(cell_width2, degree2, norm2);
    }
    double error_bound_on_deriv_2(double cell_width1, double cell_width2, int degree1, int degree2)
    {
        double norm1 = m_evaluator.max_norm(degree1 + 1, 0);
        double norm2 = m_evaluator.max_norm(0, degree2 + 1);
        return tihomirov_error_bound(cell_width1, degree1, norm1)
               + tihomirov_error_bound(cell_width2, degree2 - 1, norm2);
    }

    /*******************************************************************************
     * NOTE: The following estimates have no theoretical justification but capture
     *       the correct asympthotic rate of convergence.
     *       The error constant may be overestimated.
     *******************************************************************************/
    double error_bound_on_deriv_12(double cell_width1, double cell_width2, int degree1, int degree2)
    {
        double norm1 = m_evaluator.max_norm(degree1 + 1, 1);
        double norm2 = m_evaluator.max_norm(1, degree2 + 1);
        return tihomirov_error_bound(cell_width1, degree1 - 1, norm1)
               + tihomirov_error_bound(cell_width2, degree2 - 1, norm2);
    }

    double error_bound_on_int(double cell_width, int degree)
    {
        return tihomirov_error_bound(cell_width, degree + 1, m_evaluator.max_norm(degree + 1));
    }
};
