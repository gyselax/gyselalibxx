

# File spline\_quadrature.hpp

[**File List**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**spline\_quadrature.hpp**](spline__quadrature_8hpp.md)

[Go to the documentation of this file](spline__quadrature_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"



namespace {
template <class Grid1D>
using CoefficientFieldMem1D_host = host_t<DFieldMem<IdxRange<Grid1D>>>;
}


template <class Grid1D, class SplineBuilder>
host_t<DFieldMem<IdxRange<Grid1D>>> spline_quadrature_coefficients_1d(
        IdxRange<Grid1D> const& idx_range,
        SplineBuilder const& builder)
{
    static_assert(
            SplineBuilder::s_nbc_xmin == 0,
            "The spline quadrature requires a builder which can construct the coefficients using "
            "only the values at the interpolation points.");
    static_assert(
            SplineBuilder::s_nbc_xmax == 0,
            "The spline quadrature requires a builder which can construct the coefficients using "
            "only the values at the interpolation points.");
    // Since spline builder quadrature coeffs are not available on device, need host allocated builder.
    // See https://github.com/CExA-project/ddc/issues/598
    static_assert(
            (std::is_same_v<typename SplineBuilder::memory_space, Kokkos::HostSpace>),
            "SplineBuilder must be host allocated.");

    DFieldMem<IdxRange<Grid1D>, typename SplineBuilder::memory_space> quadrature_coefficients(
            builder.interpolation_domain());
    std::tie(std::ignore, quadrature_coefficients, std::ignore) = builder.quadrature_coefficients();
    return ddc::create_mirror_and_copy(quadrature_coefficients[idx_range]);
}



template <class ExecSpace, class... DDims, class... SplineBuilders>
DFieldMem<IdxRange<DDims...>, typename ExecSpace::memory_space> spline_quadrature_coefficients(
        IdxRange<DDims...> const& idx_range,
        SplineBuilders const&... builders)
{
    assert((std::is_same_v<
                    typename DDims::continuous_dimension_type,
                    typename SplineBuilders::continuous_dimension_type> and ...));

    // Get coefficients for each dimension
    std::tuple<CoefficientFieldMem1D_host<DDims>...> current_dim_coeffs(
            spline_quadrature_coefficients_1d(ddc::select<DDims>(idx_range), builders)...);

    // Allocate ND coefficients
    DFieldMem<IdxRange<DDims...>, typename ExecSpace::memory_space> coefficients(idx_range);
    auto coefficients_host = ddc::create_mirror(get_field(coefficients));
    // Serial loop is used due to nvcc bug concerning functions with variadic template arguments
    // (see https://github.com/kokkos/kokkos/pull/7059)
    ddc::for_each(idx_range, [&](Idx<DDims...> const idim) {
        // multiply the 1D coefficients by one another
        coefficients_host(idim)
                = (std::get<CoefficientFieldMem1D_host<DDims>>(current_dim_coeffs)(
                           ddc::select<DDims>(idim))
                   * ... * 1);
    });
    ddc::parallel_deepcopy(coefficients, coefficients_host);
    return coefficients;
}
```


