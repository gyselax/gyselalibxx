#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

struct CosineEvaluator
{
    template <class DDim>
    class Evaluator
    {
        static inline constexpr double s_2_pi = 2. * M_PI;

    private:
        double m_c0 = 1.;

        double m_c1 = 0.;

    public:
        Evaluator() = default;

        Evaluator(double c0, double c1) : m_c0(c0), m_c1(c1) {}

        double operator()(double const x) const noexcept
        {
            return eval(x, 0);
        }

        void operator()(ChunkSpan<double, DiscreteDomain<DDim>> chunk) const
        {
            auto const& domain = chunk.domain();

            for (DiscreteElement<DDim> const i : domain) {
                chunk(i) = eval(coordinate(i), 0);
            }
        }

        double deriv(double const x, int const derivative) const noexcept
        {
            return eval(x, derivative);
        }

        void deriv(ChunkSpan<double, DiscreteDomain<DDim>> chunk, int const derivative) const
        {
            auto const& domain = chunk.domain();

            for (DiscreteElement<DDim> const i : domain) {
                chunk(i) = eval(coordinate(i), derivative);
            }
        }

    private:
        double eval(double const x, int const derivative) const noexcept
        {
            return std::pow(s_2_pi * m_c0, derivative)
                   * std::cos(M_PI_2 * derivative + s_2_pi * (m_c0 * x + m_c1));
        }
    };
};
