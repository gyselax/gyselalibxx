#include <cmath>

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>
#include <geometry.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

#include "mask_tanh.hpp"

host_t<DFieldMemX> mask_tanh(
        IdxRangeX const& gridx,
        double const extent,
        double const stiffness,
        MaskType const type,
        bool const normalized)
{
    host_t<DFieldMemX> mask(gridx);

    IdxStepX const Nx(gridx.size());
    CoordX const x_min(ddc::coordinate(gridx.front()));
    double const Lx = ddcHelper::total_interval_length(gridx);
    CoordX const x_left(x_min + Lx * extent);
    CoordX const x_right(x_min + Lx - Lx * extent);

    switch (type) {
    case MaskType::Normal:
        ddc::for_each(gridx, [&](IdxX const ix) {
            CoordX const coordx = ddc::coordinate(ix);
            mask(ix) = 0.5
                       * (std::tanh((coordx - x_left) / stiffness)
                          - std::tanh((coordx - x_right) / stiffness));
        });
        break;

    case MaskType::Inverted:
        ddc::for_each(gridx, [&](IdxX const ix) {
            CoordX const coordx = ddc::coordinate(ix);
            mask(ix) = 1
                       - 0.5
                                 * (std::tanh((coordx - x_left) / stiffness)
                                    - std::tanh((coordx - x_right) / stiffness));
        });
        break;
    }

    if (normalized) {
        host_t<DFieldMemX> const quadrature_coeffs = trapezoid_quadrature_coefficients(gridx);
        host_t<Quadrature<IdxRangeX>> const integrate_x(quadrature_coeffs);
        double const coeff_norm
                = integrate_x(Kokkos::DefaultHostExecutionSpace(), get_const_field(mask));
        ddc::for_each(gridx, [&](IdxX const ix) { mask(ix) = mask(ix) / coeff_norm; });
    }

    return mask;
}
