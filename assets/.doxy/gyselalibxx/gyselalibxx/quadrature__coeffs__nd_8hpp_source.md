

# File quadrature\_coeffs\_nd.hpp

[**File List**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**quadrature\_coeffs\_nd.hpp**](quadrature__coeffs__nd_8hpp.md)

[Go to the documentation of this file](quadrature__coeffs__nd_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


namespace {
template <class ExecSpace, class Grid1D>
using CoefficientFieldMem1D = DFieldMem<IdxRange<Grid1D>, typename ExecSpace::memory_space>;
template <class ExecSpace, class Grid1D>
using CoefficientField1D = DField<IdxRange<Grid1D>, typename ExecSpace::memory_space>;
} // namespace

template <class ExecSpace, class... DDims>
DFieldMem<IdxRange<DDims...>, typename ExecSpace::memory_space> quadrature_coeffs_nd(
        IdxRange<DDims...> const& idx_range,
        std::function<DFieldMem<IdxRange<DDims>, typename ExecSpace::memory_space>(
                IdxRange<DDims>)>... funcs)
{
    DFieldMem<IdxRange<DDims...>, typename ExecSpace::memory_space> coefficients_alloc(idx_range);
    DField<IdxRange<DDims...>, typename ExecSpace::memory_space> coefficients(
            get_field(coefficients_alloc));
    // Get coefficients for each dimension
    std::tuple<CoefficientFieldMem1D<ExecSpace, DDims>...> current_dim_coeffs_alloc(
            funcs(ddc::select<DDims>(idx_range))...);
    std::tuple<CoefficientField1D<ExecSpace, DDims>...> current_dim_coeffs(get_field(
            std::get<CoefficientFieldMem1D<ExecSpace, DDims>>(current_dim_coeffs_alloc))...);

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range,
            KOKKOS_LAMBDA(Idx<DDims...> const idim) {
                // multiply the 1D coefficients by one another

                coefficients(idim)
                        = (std::get<CoefficientField1D<ExecSpace, DDims>>(current_dim_coeffs)(
                                   ddc::select<DDims>(idim))
                           * ... * 1);
            });
    return coefficients_alloc;
}
```


