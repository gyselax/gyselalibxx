#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/math_tools.hpp>

struct CosineEvaluator
{
    template <class DDim>
    class Evaluator
    {
    public:
        using Dim = DDim;

    private:
        static inline constexpr double s_2_pi = 2. * M_PI;

    private:
        double m_c0;

        double m_c1;

    public:
        template <class Domain>
        Evaluator(Domain domain) : m_c0(1.0)
                                 , m_c1(0.0)
        {
        }

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

        double max_norm(int diff = 0) const
        {
            return ipow(s_2_pi * m_c0, diff);
        }

    private:
        double eval(double const x, int const derivative) const noexcept
        {
            return ipow(s_2_pi * m_c0, derivative)
                   * std::cos(M_PI_2 * derivative + s_2_pi * (m_c0 * x + m_c1));
        }
    };
};
