// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

template <class IDim>
class Quadrature
{
private:
    ddc::Chunk<double, ddc::DiscreteDomain<IDim>> m_coefficients;

public:
    explicit Quadrature(ddc::Chunk<double, ddc::DiscreteDomain<IDim>>&& coeffs)
        : m_coefficients(std::move(coeffs))
    {
    }

    Quadrature(Quadrature&& rhs) = default;

    ~Quadrature() = default;

    double operator()(ddc::ChunkSpan<const double, ddc::DiscreteDomain<IDim>> const values) const
    {
        assert(ddc::get_domain<IDim>(values) == ddc::get_domain<IDim>(m_coefficients));
        return ddc::transform_reduce(
                ddc::policies::parallel_host,
                values.domain(),
                0.0,
                ddc::reducer::sum<double>(),
                [&](ddc::DiscreteElement<IDim> const ix) {
                    return m_coefficients(ix) * values(ix);
                });
    }
};
