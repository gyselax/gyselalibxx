#pragma once

#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

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

        KOKKOS_FUNCTION double operator()(double const x) const noexcept
        {
            return eval(x, 0);
        }

        KOKKOS_FUNCTION void operator()(
                ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> chunk) const
        {
            auto const& domain = chunk.domain();

            for (ddc::DiscreteElement<DDim> const i : domain) {
                chunk(i) = eval(ddc::coordinate(i), 0);
            }
        }

        KOKKOS_FUNCTION double deriv(double const x, int const derivative) const noexcept
        {
            return eval(x, derivative);
        }

        KOKKOS_FUNCTION void deriv(
                ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> chunk,
                int const derivative) const
        {
            auto const& domain = chunk.domain();

            for (ddc::DiscreteElement<DDim> const i : domain) {
                chunk(i) = eval(ddc::coordinate(i), derivative);
            }
        }

        KOKKOS_FUNCTION double max_norm(int diff = 0) const
        {
            return ddc::detail::ipow(s_2_pi * m_c0, diff);
        }

    private:
        KOKKOS_FUNCTION double eval(double const x, int const derivative) const noexcept
        {
            return ddc::detail::ipow(s_2_pi * m_c0, derivative)
                   * Kokkos::cos(M_PI_2 * derivative + s_2_pi * (m_c0 * x + m_c1));
        }
    };
};
