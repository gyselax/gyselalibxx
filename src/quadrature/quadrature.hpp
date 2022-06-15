// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>

#include <ddc/ddc.hpp>

template <class IDim>
class Quadrature
{
private:
    using IDomain = DiscreteDomain<IDim>;

    Chunk<double, IDomain> coefficients;

public:
    Quadrature(Chunk<double, IDomain>&& coeffs) : coefficients(std::move(coeffs)) {}
    ~Quadrature() = default;

    double operator()(ChunkSpan<double, IDomain> const& values);
};


template <class IDim>
double Quadrature<IDim>::operator()(ChunkSpan<double, DiscreteDomain<IDim>> const& values)
{
    auto domain = get_domain<IDim>(coefficients);
    assert(get_domain<IDim>(values) == domain);
    return transform_reduce(
            policies::parallel_host,
            domain,
            0.0,
            reducer::sum<double>(),
            [&](DiscreteCoordinate<IDim> const ix) { return coefficients(ix) * values(ix); });
}
