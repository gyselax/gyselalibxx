

# File ipolar\_poisson\_like\_solver.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**ipolar\_poisson\_like\_solver.hpp**](ipolar__poisson__like__solver_8hpp.md)

[Go to the documentation of this file](ipolar__poisson__like__solver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

template <
        class IdxRangeLaplacian,
        class IdxRangeFull,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutSpace = Kokkos::layout_right>
class IPolarPoissonLikeSolver;

template <class... ODims, class IdxRangeFull, class MemorySpace, class LayoutSpace>
class IPolarPoissonLikeSolver<IdxRange<ODims...>, IdxRangeFull, MemorySpace, LayoutSpace>
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
    virtual ~IPolarPoissonLikeSolver() = default;

    virtual void update_coefficients(const_field_type alpha, const_field_type beta) = 0;

    virtual void operator()(field_type phi, const_field_type rho) const = 0;
};
```


