// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

template <class IDim>
class Quadrature
{
private:
    Chunk<double, DiscreteDomain<IDim>> m_coefficients;

public:
    explicit Quadrature(Chunk<double, DiscreteDomain<IDim>>&& coeffs)
        : m_coefficients(std::move(coeffs))
    {
    }

    ~Quadrature() = default;

    double operator()(ChunkSpan<const double, DiscreteDomain<IDim>> const values) const
    {
        assert(get_domain<IDim>(values) == get_domain<IDim>(m_coefficients));
        return transform_reduce(
                policies::parallel_host,
                values.domain(),
                0.0,
                reducer::sum<double>(),
                [&](DiscreteElement<IDim> const ix) { return m_coefficients(ix) * values(ix); });
    }
};
