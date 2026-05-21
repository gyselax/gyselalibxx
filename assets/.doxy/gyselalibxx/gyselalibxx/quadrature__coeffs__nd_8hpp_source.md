

# File quadrature\_coeffs\_nd.hpp

[**File List**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**quadrature\_coeffs\_nd.hpp**](quadrature__coeffs__nd_8hpp.md)

[Go to the documentation of this file](quadrature__coeffs__nd_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


namespace {
template <class ExecSpace, class Grid1D, std::floating_point DataType>
using CoefficientFieldMem1D
        = FieldMem<DataType, IdxRange<Grid1D>, typename ExecSpace::memory_space>;
template <class ExecSpace, class Grid1D, std::floating_point DataType>
using CoefficientField1D = Field<DataType, IdxRange<Grid1D>, typename ExecSpace::memory_space>;
} // namespace

template <class ExecSpace, class DataType, class... DDims>
FieldMem<DataType, IdxRange<DDims...>, typename ExecSpace::memory_space> quadrature_coeffs_nd(
        IdxRange<DDims...> const& idx_range,
        std::function<FieldMem<DataType, IdxRange<DDims>, typename ExecSpace::memory_space>(
                IdxRange<DDims>)>... funcs)
{
    FieldMem<DataType, IdxRange<DDims...>, typename ExecSpace::memory_space>
            coefficients_alloc("coefficients (quadrature_coeffs_nd)", idx_range);
    Field<DataType, IdxRange<DDims...>, typename ExecSpace::memory_space> coefficients(
            get_field(coefficients_alloc));
    // Get coefficients for each dimension
    std::tuple<CoefficientFieldMem1D<ExecSpace, DDims, DataType>...> current_dim_coeffs_alloc(
            funcs(ddc::select<DDims>(idx_range))...);
    std::tuple<CoefficientField1D<ExecSpace, DDims, DataType>...> current_dim_coeffs(
            get_field(std::get<CoefficientFieldMem1D<ExecSpace, DDims, DataType>>(
                    current_dim_coeffs_alloc))...);

    const std::source_location location = std::source_location::current();
    ddc::parallel_for_each(
            location.function_name(),
            ExecSpace(),
            idx_range,
            KOKKOS_LAMBDA(Idx<DDims...> const idim) {
                // multiply the 1D coefficients by one another

                coefficients(idim)
                        = (std::get<CoefficientField1D<ExecSpace, DDims, DataType>>(
                                   current_dim_coeffs)(ddc::select<DDims>(idim))
                           * ... * 1);
            });
    return coefficients_alloc;
}
```


