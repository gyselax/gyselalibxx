

# File ipoisson\_solver.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**ipoisson\_solver.hpp**](ipoisson__solver_8hpp.md)

[Go to the documentation of this file](ipoisson__solver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "vector_field.hpp"

template <class IdxRangeLaplacian, class IdxRangeFull, class MemorySpace, class LayoutSpace>
class IPoissonSolver;

template <class... ODims, class IdxRangeFull, class MemorySpace, class LayoutSpace>
class IPoissonSolver<IdxRange<ODims...>, IdxRangeFull, MemorySpace, LayoutSpace>
{
protected:
    using real_laplacian_tags = ddc::detail::TypeSeq<typename ODims::continuous_dimension_type...>;
    using laplacian_tags = ddc::detail::TypeSeq<ODims...>;
    using space_tags = ddc::to_type_seq_t<IdxRangeFull>;
    using batch_tags = ddc::type_seq_remove_t<space_tags, laplacian_tags>;

protected:
    static constexpr bool using_vector_field = ddc::type_seq_size_v<laplacian_tags> == 1;

public:
    using field_type = DField<IdxRangeFull, MemorySpace, LayoutSpace>;
    using const_field_type = DConstField<IdxRangeFull, MemorySpace, LayoutSpace>;

    using vector_field_type = std::conditional_t<
            ddc::type_seq_size_v<laplacian_tags> == 1,
            field_type,
            VectorField<double, IdxRangeFull, real_laplacian_tags, MemorySpace, LayoutSpace>>;

    using batch_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<batch_tags>;
    using batch_index_type = typename batch_idx_range_type::discrete_element_type;

    using laplacian_idx_range_type = IdxRange<ODims...>;

    using layout_space = LayoutSpace;
    using memory_space = MemorySpace;

public:
    virtual field_type operator()(field_type phi, field_type rho) const = 0;

    virtual field_type operator()(field_type phi, vector_field_type E, field_type rho) const = 0;
};
```


