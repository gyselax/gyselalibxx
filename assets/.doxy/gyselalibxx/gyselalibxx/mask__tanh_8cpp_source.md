

# File mask\_tanh.cpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**mask\_tanh.cpp**](mask__tanh_8cpp.md)

[Go to the documentation of this file](mask__tanh_8cpp.md)


```C++
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "geometry.hpp"
#include "mask_tanh.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

host_t<DFieldMemX> mask_tanh(
        IdxRangeX const& gridx,
        double const extent,
        double const stiffness,
        MaskType const type,
        bool const normalised)
{
    if (stiffness <= 0) {
        throw std::runtime_error("Invalid stiffness, cannot be negative");
    }
    if (extent < 0 || extent > 0.5) {
        throw std::runtime_error("Invalid extent, cannot be more than 0.5");
    }
    host_t<DFieldMemX> mask(gridx);
    if (extent == 0) {
        ddc::parallel_fill(mask, 0.);
    } else {
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

        if (normalised) {
            host_t<DFieldMemX> const quadrature_coeffs
                    = trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(gridx);
            host_t<Quadrature<IdxRangeX>> const integrate_x(get_const_field(quadrature_coeffs));
            double const coeff_norm
                    = integrate_x(Kokkos::DefaultHostExecutionSpace(), get_const_field(mask));
            ddc::for_each(gridx, [&](IdxX const ix) { mask(ix) = mask(ix) / coeff_norm; });
        }
    }

    return mask;
}
```


