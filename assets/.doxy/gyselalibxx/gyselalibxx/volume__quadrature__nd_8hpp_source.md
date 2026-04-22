

# File volume\_quadrature\_nd.hpp

[**File List**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**volume\_quadrature\_nd.hpp**](volume__quadrature__nd_8hpp.md)

[Go to the documentation of this file](volume__quadrature__nd_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "quadrature.hpp"



template <concepts::MappingWithJacobian Mapping, class IdxRangeCoeffs, class ExecSpace>
DFieldMem<IdxRangeCoeffs, typename ExecSpace::memory_space> compute_coeffs_on_mapping(
        ExecSpace exec_space,
        Mapping& mapping,
        DFieldMem<IdxRangeCoeffs, typename ExecSpace::memory_space>&& coefficients_alloc)
{
    // Get type information about Jacobian definition
    using CoordJacobian = typename Mapping::CoordJacobian;
    using IdxJacobian = find_idx_t<CoordJacobian, IdxRangeCoeffs>;
    using TypeSeqJacobian = ddc::to_type_seq_t<IdxJacobian>;
    // Strip any Grids not found in the index range
    using TypeSeqJ = ddc::type_seq_remove_t<TypeSeqJacobian, ddc::detail::TypeSeq<void>>;
    using IdxRangeJ = ddc::detail::convert_type_seq_to_discrete_domain_t<TypeSeqJ>;
    using IdxJ = typename IdxRangeJ::discrete_element_type;
    using CoordJ = ddc::coordinate_of_t<IdxJ>;
    // Check that the resulting coordinate is valid for retrieving the determinant of the Jacobian
    static_assert(
            concepts::MappingWithIntegrationCoord<Mapping, CoordJ>,
            "The provided index range does not give enough positional information to define the "
            "determinant of the Jacobian");

    using IdxCoeffs = typename IdxRangeCoeffs::discrete_element_type;
    IdxRangeCoeffs grid = get_idx_range(coefficients_alloc);
    DField<IdxRangeCoeffs, typename ExecSpace::memory_space> coefficients(coefficients_alloc);
    const std::source_location location = std::source_location::current();
    ddc::parallel_for_each(
            location.function_name(),
            exec_space,
            grid,
            KOKKOS_LAMBDA(IdxCoeffs const idx) {
                CoordJ coord(ddc::coordinate(IdxJ(idx)));
                coefficients(idx) *= fabs(mapping.jacobian(coord));
            });
    return std::move(coefficients_alloc);
}
```


