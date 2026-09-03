

# File discrete\_mapping.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**discrete\_mapping.hpp**](discrete__mapping_8hpp.md)

[Go to the documentation of this file](discrete__mapping_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "i_interpolation_evaluator.hpp"
#include "math_tools.hpp"
#include "tensor.hpp"
#include "tensor_common.hpp"
#include "vector_field.hpp"
#include "vector_field_evaluation.hpp"

namespace details {

template <
        class RowDim,
        class DataType,
        class... RDim,
        class... ADim,
        class Mapping,
        class CoordJacobian>
KOKKOS_INLINE_FUNCTION void fill_jacobian_matrix_row(
        Tensor<DataType, VectorIndexSet<RDim...>, VectorIndexSet<ADim...>>& jacobian_matrix,
        Mapping const& mapping,
        CoordJacobian const& coord)
{
    ((ddcHelper::get<RowDim, ADim>(jacobian_matrix)
      = mapping.template jacobian_component<RowDim, ADim>(coord)),
     ...);
}

template <class DataType, class... RDim, class... ADim, class Mapping, class CoordJacobian>
KOKKOS_INLINE_FUNCTION void fill_jacobian_matrix(
        Tensor<DataType, VectorIndexSet<RDim...>, VectorIndexSet<ADim...>>& jacobian_matrix,
        Mapping const& mapping,
        CoordJacobian const& coord)
{
    ((fill_jacobian_matrix_row<RDim>(jacobian_matrix, mapping, coord)), ...);
}

} // namespace details

template <class StartCoord, class EndCoord, concepts::InterpolationEvaluator NDEvaluator>
class DiscreteMapping
{
    static_assert(InterpolationEvaluatorTraits<NDEvaluator>::rank() == StartCoord::size());

public:
    using CoordArg = StartCoord;
    using CoordResult = EndCoord;
    using CoordJacobian = StartCoord;

    using DataType = typename InterpolationEvaluatorTraits<NDEvaluator>::data_type;
    using ArgBasis = ddc::to_type_seq_t<StartCoord>;
    using ResultBasis = ddc::to_type_seq_t<EndCoord>;

private:
    using CoeffField = VectorConstField<
            DataType,
            typename InterpolationEvaluatorTraits<NDEvaluator>::coeff_idx_range_type,
            ResultBasis,
            typename NDEvaluator::memory_space>;

private:
    CoeffField m_coeff_representation;
    NDEvaluator m_evaluator;

public:
    DiscreteMapping(CoeffField coeff_representation, NDEvaluator const& evaluator)
        : m_coeff_representation(coeff_representation)
        , m_evaluator(evaluator)
    {
    }

    KOKKOS_DEFAULTED_FUNCTION DiscreteMapping(DiscreteMapping const& other) = default;

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        return CoordResult(
                ddcHelper::to_coord(ndEval::evaluate(m_evaluator, coord, m_coeff_representation)));
    }

    template <class ExecSpace, class... GridType>
    void operator()(ExecSpace exec_space, Field<CoordResult, IdxRange<GridType...>> coords)
    {
        static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename NDEvaluator::memory_space>::
                              accessible);
        static_assert(((ddc::in_tags_v<
                        typename GridType::continuous_dimension_type,
                        ddc::to_type_seq_t<StartCoord>>)&&...));
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space,
                get_idx_range(coords),
                KOKKOS_CLASS_LAMBDA(Idx<GridType...> idx) {
                    coords(idx) = (*this)(ddc::coordinate(idx));
                });
    }

    KOKKOS_FUNCTION Tensor<DataType, ResultBasis, get_covariant_dims_t<ArgBasis>> jacobian_matrix(
            CoordJacobian const& coord) const
    {
        Tensor<DataType, ResultBasis, get_covariant_dims_t<ArgBasis>> jacobian_matrix;
        details::fill_jacobian_matrix(jacobian_matrix, *this, coord);
        return jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(CoordJacobian coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ddc::to_type_seq_t<EndCoord>>);
        static_assert(
                ddc::in_tags_v<IndexTag2, get_covariant_dims_t<ddc::to_type_seq_t<StartCoord>>>);

        return m_evaluator
                .deriv(Idx<ddc::Deriv<typename IndexTag2::Dual>>(1),
                       coord,
                       get_const_field(ddcHelper::get<IndexTag1>(m_coeff_representation)));
    }

    KOKKOS_FUNCTION double jacobian(CoordJacobian const& coord) const
    {
        Tensor J = jacobian_matrix(coord);
        return determinant(J);
    }
};


namespace mapping_detail {
template <class StartCoord, class EndCoord, class NDEvaluator, class ExecSpace>
struct MappingAccessibility<ExecSpace, DiscreteMapping<StartCoord, EndCoord, NDEvaluator>>
{
    static constexpr bool value
            = Kokkos::SpaceAccessibility<ExecSpace, typename NDEvaluator::memory_space>::accessible;
};

} // namespace mapping_detail
```


