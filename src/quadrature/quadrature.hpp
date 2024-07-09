// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

/**
 * @brief A class providing an operator for integrating functions defined on a discrete domain.
 */
template <class ExecSpace, class... IDim>
class Quadrature
{
private:
    ddc::ChunkSpan<
            const double,
            ddc::DiscreteDomain<IDim...>,
            std::experimental::layout_right,
            typename ExecSpace::memory_space>
            m_coefficients;

public:
    /**
     * @brief Create a Quadrature object.
     * @param coeffs
     * 	      The coefficients of the quadrature.
     */
    explicit Quadrature(ddc::ChunkSpan<
                        double,
                        ddc::DiscreteDomain<IDim...>,
                        std::experimental::layout_right,
                        typename ExecSpace::memory_space> const& coeffs)
        : m_coefficients(coeffs)
    {
    }

    /**
     * @brief An operator for calculating the integral of a function defined on a discrete domain.
     * @param[in] values
     *        The values of the function on the points of the discrete domain.
     *
     * @returns The integral of the function over the domain.
     */
    double operator()(ddc::ChunkView<
                      double,
                      ddc::DiscreteDomain<IDim...>,
                      std::experimental::layout_right,
                      typename ExecSpace::memory_space> const values) const
    {
        assert(ddc::get_domain<IDim...>(values) == ddc::get_domain<IDim...>(m_coefficients));
        auto coeff_proxy = m_coefficients;
        return ddc::parallel_transform_reduce(
                Kokkos::DefaultExecutionSpace(),
                values.domain(),
                0.0,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(ddc::DiscreteElement<IDim...> const ix) {
                    return coeff_proxy(ix) * values(ix);
                });
    }
};
