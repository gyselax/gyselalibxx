#pragma once

#include <cmath>
#include <random>
#include <vector>

#include <ddc/ddc.hpp>

struct PolynomialEvaluator
{
    template <class DDim>
    class Evaluator
    {
        static inline constexpr double s_2_pi = 2. * M_PI;

    private:
        std::vector<double> m_coeffs;
        int const m_degree;

    public:
        Evaluator(int max_degree) : m_coeffs(max_degree + 1, 0.0), m_degree(max_degree)
        {
            for (int i(0); i < max_degree + 1; ++i) {
                m_coeffs[i] = double(rand() % 100) / 100.0;
            }
        }

        double operator()(double const x) const noexcept
        {
            return eval(x, 0);
        }

        void operator()(ChunkSpan<double, DiscreteDomain<DDim>> chunk) const
        {
            auto const& domain = chunk.domain();

            for (DiscreteCoordinate<DDim> const i : domain) {
                chunk(i) = eval(to_real(i), 0);
            }
        }

        double deriv(double const x, int const derivative) const noexcept
        {
            return eval(x, derivative);
        }

        void deriv(ChunkSpan<double, DiscreteDomain<DDim>> chunk, int const derivative) const
        {
            auto const& domain = chunk.domain();

            for (DiscreteCoordinate<DDim> const i : domain) {
                chunk(i) = eval(to_real(i), derivative);
            }
        }

    private:
        double eval(double const x, int const derivative) const
        {
            double result(0.0);
            int start = derivative < 0 ? 0 : derivative;
            for (int i(start); i < m_degree; ++i) {
                double v = double(falling_factorial(i, derivative)) * std::pow(x, i - derivative);
                result += m_coeffs[i] * v;
            }
            return result;
        }

        double falling_factorial(int i, int d) const
        {
            double c = 1.0;
            if (d >= 0) {
                for (int k(0); k < d; ++k) {
                    c *= (i - k);
                }
            } else {
                for (int k(-1); k > d - 1; --k) {
                    c /= (i - k);
                }
            }
            return c;
        }
    };
};
