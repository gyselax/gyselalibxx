

# File identity\_coordinate\_change.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**identity\_coordinate\_change.hpp**](identity__coordinate__change_8hpp.md)

[Go to the documentation of this file](identity__coordinate__change_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

template <class ArgBasis, class ResultBasis>
class IdentityCoordinateChange
{
public:
    using CoordArg = ddcHelper::to_coord_t<ArgBasis>;
    using CoordResult = ddcHelper::to_coord_t<ResultBasis>;
    using CoordJacobian = ddcHelper::to_coord_t<ArgBasis>;

private:
    using ArgBasisCov = vector_index_set_dual_t<ArgBasis>;

public:
    template <class... ArgDims>
    KOKKOS_FUNCTION CoordResult operator()(Coord<ArgDims...> const& coord) const
    {
        static_assert(std::is_same_v<CoordArg, Coord<ArgDims...>>);
        return CoordResult((ddc::get<ArgDims>(coord))...);
    }

    KOKKOS_INLINE_FUNCTION double jacobian(CoordArg const& coord) const
    {
        return 1;
    }

    KOKKOS_FUNCTION DTensor<ResultBasis, ArgBasisCov> jacobian_matrix(CoordArg const& coord) const
    {
        DTensor<ResultBasis, ArgBasisCov> jacobian_matrix(0.0);
        fill_diagonal_elements(jacobian_matrix);
        return jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ResultBasis>);
        static_assert(ddc::in_tags_v<IndexTag2, ArgBasisCov>);
        if constexpr (
                (ddc::type_seq_rank_v<IndexTag1, ResultBasis>)
                == (ddc::type_seq_rank_v<IndexTag2, ArgBasisCov>)) {
            return 1;
        } else {
            return 0;
        }
    }

    KOKKOS_INLINE_FUNCTION double inv_jacobian(CoordArg const& coord) const
    {
        return 1;
    }

    KOKKOS_FUNCTION DTensor<ResultBasis, ArgBasisCov> inv_jacobian_matrix(
            CoordArg const& coord) const
    {
        DTensor<ResultBasis, ArgBasisCov> inv_jacobian_matrix(0.0);
        fill_diagonal_elements(inv_jacobian_matrix);
        return inv_jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double inv_jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ResultBasis>);
        static_assert(ddc::in_tags_v<IndexTag2, ArgBasisCov>);
        if constexpr (
                (ddc::type_seq_rank_v<IndexTag1, ResultBasis>)
                == (ddc::type_seq_rank_v<IndexTag2, ArgBasisCov>)) {
            return 1;
        } else {
            return 0;
        }
    }

    KOKKOS_INLINE_FUNCTION IdentityCoordinateChange<ResultBasis, ArgBasis> get_inverse_mapping()
            const
    {
        return IdentityCoordinateChange<ResultBasis, ArgBasis>();
    }

private:
    template <class... RowDims, class... ColDims>
    void fill_diagonal_elements(
            DTensor<VectorIndexSet<RowDims...>, VectorIndexSet<ColDims...>>& matrix) const
    {
        ((ddcHelper::get<RowDims, ColDims>(matrix) = 1.0), ...);
    }
};

namespace mapping_detail {
template <class ArgBasis, class ResultBasis, class ExecSpace>
struct MappingAccessibility<ExecSpace, IdentityCoordinateChange<ArgBasis, ResultBasis>>
    : std::true_type
{
};
} // namespace mapping_detail
```


