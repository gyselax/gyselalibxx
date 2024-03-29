#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

#include <ddc/ddc.hpp>

#include <sll/math_tools.hpp>

struct PolynomialEvaluator
{
    template <class DDim, std::size_t Degree>
    class Evaluator
    {
    public:
        using Dim = DDim;

    private:
        std::array<double, Degree + 1> m_coeffs;
        int const m_degree;
        double const m_xN;

    public:
        template <class Domain>
        Evaluator(Domain domain)
            : m_degree(Degree)
            , m_xN(std::max(std::abs(rmin(domain)), std::abs(rmax(domain))))
        {
            for (int i(0); i < m_degree + 1; ++i) {
                m_coeffs[i] = double(rand() % 100) / 100.0;
            }
        }

        double operator()(double const x) const noexcept
        {
            return eval(x, 0);
        }

        void operator()(ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> chunk) const
        {
            auto const& domain = chunk.domain();

            for (ddc::DiscreteElement<DDim> const i : domain) {
                chunk(i) = eval(ddc::coordinate(i), 0);
            }
        }

        double deriv(double const x, int const derivative) const noexcept
        {
            return eval(x, derivative);
        }

        void deriv(ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> chunk, int const derivative)
                const
        {
            auto const& domain = chunk.domain();

            for (ddc::DiscreteElement<DDim> const i : domain) {
                chunk(i) = eval(ddc::coordinate(i), derivative);
            }
        }

        double max_norm(int diff = 0) const
        {
            return std::abs(deriv(m_xN, diff));
        }

    private:
        double eval(double const x, int const derivative) const
        {
            double result(0.0);
            int start = derivative < 0 ? 0 : derivative;
            for (int i(start); i < m_degree + 1; ++i) {
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
