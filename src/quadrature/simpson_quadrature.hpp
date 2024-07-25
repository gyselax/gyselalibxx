// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>


/**
 * @brief Get the Simpson coefficients in 1D.
 *
 * Calculate the quadrature coefficients for the Simpson method defined on the provided domain.
 *
 * @tparam ExecSpace Execution space, depends on Kokkos.
 *
 * @param[in] domain The domain on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the Simpson method defined on the provided domain.
 *         The allocation place (host or device ) will depend on the ExecSpace.
 */
template <class ExecSpace, class IDim>
ddc::Chunk<
        double,
        ddc::DiscreteDomain<IDim>,
        ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
simpson_quadrature_coefficients_1d(ddc::DiscreteDomain<IDim> const& domain)
{
    ddc::Chunk<
            double,
            ddc::DiscreteDomain<IDim>,
            ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
            coefficients_alloc(domain);
    ddc::ChunkSpan const coefficients = coefficients_alloc.span_view();
    double const dx_l = distance_at_left(domain.back());
    double const dx_r = distance_at_right(domain.front());
    Kokkos::parallel_for(
            "bounds",
            Kokkos::RangePolicy<ExecSpace>(0, 1),
            KOKKOS_LAMBDA(const int i) {
                coefficients(domain.front()) = 1. / 3. * dx_r;
                coefficients(domain.back()) = 1. / 3. * dx_l;
                for (auto it = domain.begin() + 1; it < domain.end() - 1; it += 2) {
                    ddc::DiscreteElement<IDim> idx = *it;
                    coefficients(idx) = 2. / 3. * (distance_at_left(idx) + distance_at_right(idx));
                    idx += ddc::DiscreteVector<IDim>(1);
                    coefficients(idx) = 1. / 3. * (distance_at_left(idx) + distance_at_right(idx));
                }
                if (IDim::continuous_dimension_type::PERIODIC) {
                    coefficients(domain.front()) += 2. / 3. * dx_l;
                    coefficients(domain.back()) += 2. / 3. * dx_r;
                }
            });
    return coefficients_alloc;
}