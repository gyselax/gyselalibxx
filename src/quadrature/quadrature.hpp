// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

/**
 * @brief A class providing an operator for integrating functions defined on a discrete domain.
 */
template <class... IDim>
class Quadrature
{
private:
    ddc::Chunk<double, ddc::DiscreteDomain<IDim...>> m_coefficients;

public:
    /**
     * @brief Create a Quadrature object.
     * @param coeffs
     * 	      The coefficients of the quadrature.
     */
    explicit Quadrature(ddc::Chunk<double, ddc::DiscreteDomain<IDim...>>&& coeffs)
        : m_coefficients(std::move(coeffs))
    {
    }

    /**
     * @brief Create a Quadrature object by copy.
     * @param rhs The object being copied.
     */
    Quadrature(Quadrature&& rhs) = default;

    ~Quadrature() = default;

    /**
     * @brief An operator for calculating the integral of a function defined on a discrete domain.
     * @param[in] values
     *        The values of the function on the points of the discrete domain.
     *
     * @returns The integral of the function over the domain.
     */
    double operator()(ddc::ChunkSpan<const double, ddc::DiscreteDomain<IDim...>> const values) const
    {
        assert(ddc::get_domain<IDim...>(values) == ddc::get_domain<IDim...>(m_coefficients));
        return ddc::transform_reduce(
                values.domain(),
                0.0,
                ddc::reducer::sum<double>(),
                [&](ddc::DiscreteElement<IDim...> const ix) {
                    return m_coefficients(ix) * values(ix);
                });
    }
};
