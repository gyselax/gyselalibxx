

# File vector\_mapper.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**vector\_mapper.hpp**](vector__mapper_8hpp.md)

[Go to the documentation of this file](vector__mapper_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "inverse_jacobian_matrix.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

namespace vector_mapper_detail {
template <class OutVectorSpace, class InVectorSpace>
static constexpr bool is_variance_change_only_v
        = !std::is_same_v<
                  OutVectorSpace,
                  InVectorSpace> && std::is_same_v<OutVectorSpace, vector_index_set_dual_t<InVectorSpace>>;

template <class OutVectorSpace, class InVectorSpace>
static constexpr bool is_coordinate_change_only_v
        = !std::is_same_v<
                  OutVectorSpace,
                  InVectorSpace> && !is_variance_change_only_v<OutVectorSpace, InVectorSpace> && has_same_variance_v<InVectorSpace, OutVectorSpace>;

template <class OutVectorSpace, class InVectorSpace, class Mapping>
static constexpr bool is_arg_to_result_v
        = is_coordinate_change_only_v<OutVectorSpace, InVectorSpace>&& std::is_same_v<
                get_contravariant_dims_t<InVectorSpace>,
                ddc::to_type_seq_t<typename Mapping::CoordArg>>&&
                std::is_same_v<
                        get_contravariant_dims_t<OutVectorSpace>,
                        ddc::to_type_seq_t<typename Mapping::CoordResult>>;

template <class OutVectorSpace, class InVectorSpace, class Mapping>
static constexpr bool is_result_to_arg_v
        = is_coordinate_change_only_v<OutVectorSpace, InVectorSpace>&& std::is_same_v<
                get_contravariant_dims_t<InVectorSpace>,
                ddc::to_type_seq_t<typename Mapping::CoordResult>>&&
                std::is_same_v<
                        get_contravariant_dims_t<OutVectorSpace>,
                        ddc::to_type_seq_t<typename Mapping::CoordArg>>;

template <class OutVectorSpace, class InVectorSpace>
static constexpr bool are_contravariant_v = is_contravariant_vector_index_set_v<InVectorSpace>&&
        is_contravariant_vector_index_set_v<OutVectorSpace>;

template <class MappingOrMetric>
struct EvaluationCoord;

template <concepts::MetricTensor MetricTensorType>
struct EvaluationCoord<MetricTensorType>
{
    using type = typename MetricTensorType::CoordArg;
};

template <concepts::Mapping Mapping>
requires requires
{
    typename Mapping::CoordJacobian;
}
struct EvaluationCoord<Mapping>
{
    using type = typename Mapping::CoordJacobian;
};

template <class MappingOrMetric>
using evaluation_coord_t = typename EvaluationCoord<MappingOrMetric>::type;
} // namespace vector_mapper_detail

template <
        class OutVectorSpace,
        class MappingOrMetric,
        class CoordType,
        class ElementType,
        class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, OutVectorSpace> to_vector_space(
        [[maybe_unused]] MappingOrMetric const& mapping_or_metric,
        [[maybe_unused]] CoordType const& coord,
        Tensor<ElementType, InVectorSpace> const&
                in_vector) requires(std::is_same_v<OutVectorSpace, InVectorSpace>)
{
    return in_vector;
}

template <
        class OutVectorSpace,
        concepts::MetricTensor MetricTensorType,
        class CoordType,
        class ElementType,
        class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, OutVectorSpace> to_vector_space(
        MetricTensorType const& metric_tensor,
        CoordType const& coord,
        Tensor<ElementType, InVectorSpace> const&
                in_vector) requires(vector_mapper_detail::
                                            is_variance_change_only_v<
                                                    OutVectorSpace,
                                                    InVectorSpace>&&
                                                    is_contravariant_vector_index_set_v<
                                                            InVectorSpace>)
{
    return tensor_mul(index<'i', 'j'>(metric_tensor(coord)), index<'j'>(in_vector));
}

template <
        class OutVectorSpace,
        concepts::MetricTensor MetricTensorType,
        class CoordType,
        class ElementType,
        class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, OutVectorSpace>
to_vector_space(MetricTensorType const& metric_tensor, CoordType const& coord, Tensor<ElementType, InVectorSpace> const& in_vector) requires(
        vector_mapper_detail::is_variance_change_only_v<
                OutVectorSpace,
                InVectorSpace> && !is_contravariant_vector_index_set_v<InVectorSpace>)
{
    return tensor_mul(index<'i', 'j'>(metric_tensor.inverse(coord)), index<'j'>(in_vector));
}

template <
        class OutVectorSpace,
        concepts::MappingWithJacobian Mapping,
        class CoordType,
        class ElementType,
        class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, OutVectorSpace> to_vector_space(
        Mapping const& mapping,
        CoordType const& coord,
        Tensor<ElementType, InVectorSpace> const&
                in_vector) requires(vector_mapper_detail::
                                            is_arg_to_result_v<
                                                    OutVectorSpace,
                                                    InVectorSpace,
                                                    Mapping>&&
                                                    vector_mapper_detail::are_contravariant_v<
                                                            OutVectorSpace,
                                                            InVectorSpace>)
{
    // A_{\{p\}}^i = J_{\{q\rightarrow p\}}^i_j A_{\{q\}}^j
    return tensor_mul(index<'i', 'j'>(mapping.jacobian_matrix(coord)), index<'j'>(in_vector));
}

template <
        class OutVectorSpace,
        concepts::MappingWithJacobian Mapping,
        class CoordType,
        class ElementType,
        class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, OutVectorSpace>
to_vector_space(Mapping const& mapping, CoordType const& coord, Tensor<ElementType, InVectorSpace> const& in_vector) requires(
        vector_mapper_detail::is_arg_to_result_v<
                OutVectorSpace,
                InVectorSpace,
                Mapping> && !vector_mapper_detail::are_contravariant_v<OutVectorSpace, InVectorSpace>)
{
    // A_{\{p\} i} = (J_{\{q\rightarrow p\}}^{-T})_i^j A_{\{q\} j}
    //             = (J_{\{q\rightarrow p\}}^{-1})_j^i A_{\{q\} j}
    InverseJacobianMatrix inv_jacobian(mapping);
    return tensor_mul(index<'j', 'i'>(inv_jacobian(coord)), index<'j'>(in_vector));
}

template <
        class OutVectorSpace,
        concepts::MappingWithJacobian Mapping,
        class CoordType,
        class ElementType,
        class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, OutVectorSpace> to_vector_space(
        Mapping const& mapping,
        CoordType const& coord,
        Tensor<ElementType, InVectorSpace> const&
                in_vector) requires(vector_mapper_detail::
                                            is_result_to_arg_v<
                                                    OutVectorSpace,
                                                    InVectorSpace,
                                                    Mapping>&&
                                                    vector_mapper_detail::are_contravariant_v<
                                                            OutVectorSpace,
                                                            InVectorSpace>)
{
    // A_{\{q\}}^i = (J_{\{q\rightarrow p\}}^{-1})^i_j A_{\{q\}}^j
    InverseJacobianMatrix inv_jacobian(mapping);
    DTensor<OutVectorSpace, get_covariant_dims_t<InVectorSpace>> I_J = inv_jacobian(coord);
    return tensor_mul(index<'i', 'j'>(I_J), index<'j'>(in_vector));
}

template <
        class OutVectorSpace,
        concepts::MappingWithJacobian Mapping,
        class CoordType,
        class ElementType,
        class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, OutVectorSpace>
to_vector_space(Mapping const& mapping, CoordType const& coord, Tensor<ElementType, InVectorSpace> const& in_vector) requires(
        vector_mapper_detail::is_result_to_arg_v<
                OutVectorSpace,
                InVectorSpace,
                Mapping> && !vector_mapper_detail::are_contravariant_v<OutVectorSpace, InVectorSpace>)
{
    // A_{\{p\} i} = (J_{\{q\rightarrow p\}}^{-T})_i^j A_{\{q\} j}
    //             = (J_{\{q\rightarrow p\}}^{-1})_j^i A_{\{q\} j}
    //             = (J_{\{p\rightarrow q\}})_j^i A_{\{q\} j}
    return tensor_mul(index<'j', 'i'>(mapping.jacobian_matrix(coord)), index<'j'>(in_vector));
}

template <
        class OutVectorSpace,
        class ExecSpace,
        class MappingOrMetric,
        class ElementType,
        class IdxRangeType,
        class InVectorSpace,
        class LayoutStridedPolicy>
void copy_to_vector_space(
        ExecSpace exec_space,
        VectorField<
                ElementType,
                IdxRangeType,
                OutVectorSpace,
                typename ExecSpace::memory_space,
                LayoutStridedPolicy> vector_field_out,
        MappingOrMetric mapping_or_metric,
        VectorConstField<
                ElementType,
                IdxRangeType,
                InVectorSpace,
                typename ExecSpace::memory_space,
                LayoutStridedPolicy> vector_field)
{
    using IdxType = typename IdxRangeType::discrete_element_type;
    using CoordType = vector_mapper_detail::evaluation_coord_t<MappingOrMetric>;
    using IdxCoordType = find_idx_t<CoordType, IdxRangeType>;
    const std::source_location location = std::source_location::current();
    ddc::parallel_for_each(
            location.function_name(),
            exec_space,
            get_idx_range(vector_field),
            KOKKOS_LAMBDA(IdxType idx) {
                IdxCoordType coord_idx(idx);
                CoordType coord = ddc::coordinate(coord_idx);
                ddcHelper::assign_vector_field_element(
                        vector_field_out,
                        idx,
                        to_vector_space<
                                OutVectorSpace>(mapping_or_metric, coord, vector_field(idx)));
            });
}

template <
        class OutVectorSpace,
        class ExecSpace,
        class MappingOrMetric,
        class ElementType,
        class IdxRangeType,
        class InVectorSpace,
        class LayoutStridedPolicy>
auto create_mirror_view_and_copy_on_vector_space(
        std::string const& label,
        ExecSpace exec_space,
        VectorField<
                ElementType,
                IdxRangeType,
                InVectorSpace,
                typename ExecSpace::memory_space,
                LayoutStridedPolicy> vector_field,
        MappingOrMetric mapping_or_metric)
{
    if constexpr (std::is_same_v<InVectorSpace, OutVectorSpace>) {
        return vector_field;
    } else {
        VectorFieldMem<
                std::remove_const_t<ElementType>,
                IdxRangeType,
                OutVectorSpace,
                typename ExecSpace::memory_space>
                vector_field_out(label, get_idx_range(vector_field));
        copy_to_vector_space(
                exec_space,
                get_field(vector_field_out),
                mapping_or_metric,
                get_const_field(vector_field));
        return vector_field_out;
    }
}

template <
        class OutVectorSpace,
        class ExecSpace,
        class MappingOrMetric,
        class ElementType,
        class IdxRangeType,
        class InVectorSpace,
        class LayoutStridedPolicy>
auto create_mirror_view_and_copy_on_vector_space(
        ExecSpace exec_space,
        VectorField<
                ElementType,
                IdxRangeType,
                InVectorSpace,
                typename ExecSpace::memory_space,
                LayoutStridedPolicy> vector_field,
        MappingOrMetric mapping_or_metric)
{
    return create_mirror_view_and_copy_on_vector_space<
            OutVectorSpace,
            ExecSpace,
            MappingOrMetric,
            ElementType,
            IdxRangeType,
            InVectorSpace,
            LayoutStridedPolicy>(
            "_ (create_mirror_view_and_copy_on_vector_space)",
            exec_space,
            vector_field,
            mapping_or_metric);
}
```


