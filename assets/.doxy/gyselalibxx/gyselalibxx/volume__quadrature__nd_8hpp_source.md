

# File volume\_quadrature\_nd.hpp

[**File List**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**volume\_quadrature\_nd.hpp**](volume__quadrature__nd_8hpp.md)

[Go to the documentation of this file](volume__quadrature__nd_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "quadrature.hpp"



template <class Mapping, class IdxRangeCoeffs, class ExecSpace>
DFieldMem<IdxRangeCoeffs, typename ExecSpace::memory_space> compute_coeffs_on_mapping(
        ExecSpace exec_space,
        Mapping& mapping,
        DFieldMem<IdxRangeCoeffs, typename ExecSpace::memory_space>&& coefficients_alloc)
{
    static_assert(is_curvilinear_2d_mapping_v<Mapping>);

    using R = typename Mapping::curvilinear_tag_r;
    using Theta = typename Mapping::curvilinear_tag_theta;

    static_assert(has_jacobian_v<Mapping, Coord<R, Theta>>);

    using IdxCoeffs = typename IdxRangeCoeffs::discrete_element_type;
    IdxRangeCoeffs grid = get_idx_range(coefficients_alloc);
    DField<IdxRangeCoeffs, typename ExecSpace::memory_space> coefficients(coefficients_alloc);
    ddc::parallel_for_each(
            exec_space,
            grid,
            KOKKOS_LAMBDA(IdxCoeffs const idx) {
                Coord<R, Theta> coord(ddc::coordinate(idx));
                coefficients(idx) *= fabs(mapping.jacobian(coord));
            });
    return std::move(coefficients_alloc);
}
```


