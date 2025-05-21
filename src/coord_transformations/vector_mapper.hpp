// SPDX-License-Identifier: MIT
#pragma once

#include "inverse_jacobian_matrix.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

/**
 * @brief A helper method to get a vector on a different vector space.
 *
 * @param[in] mapping A mapping describing the relation between the 2 vector
 *              spaces (or describing the relation between the vector space and
 *              the Cartesian space if a change of variance is required.
 * @param[in] coord The coordinate where the vector is defined.
 * @param[in] in_vector The vector on the original vector space.
 *
 * @tparam OutVectorSpace The vector space of the final vector.
 *
 * @returns A VectorField or VectorFieldMem containing the vectors in the
 *      requested vector space.
 */
template <
        class OutVectorSpace,
        class Mapping,
        class CoordType,
        class ElementType,
        class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, OutVectorSpace> to_vector_space(
        Mapping const& mapping,
        CoordType const& coord,
        Tensor<ElementType, InVectorSpace> const& in_vector)
{
    if constexpr (std::is_same_v<OutVectorSpace, InVectorSpace>) {
        // If in and out tensor are already defined on the same basis
        return in_vector;

    } else if constexpr (std::is_same_v<OutVectorSpace, vector_index_set_dual_t<InVectorSpace>>) {
        // If in and out tensor are defined on co/contra variant bases associated with the same coordinate system
        if constexpr (is_contravariant_vector_index_set_v<InVectorSpace>) {
            MetricTensorEvaluator<Mapping, CoordType> metric(mapping);
            return tensor_mul(index<'i', 'j'>(metric(coord)), index<'j'>(in_vector));
        } else {
            MetricTensorEvaluator<Mapping, CoordType> metric(mapping);
            return tensor_mul(index<'i', 'j'>(metric.inverse(coord)), index<'j'>(in_vector));
        }

    } else if constexpr (!has_same_variance_v<InVectorSpace, OutVectorSpace>) {
        // If different variance (co/contra)
        using OutVectorSpaceDual = vector_index_set_dual_t<OutVectorSpace>;
        Tensor<ElementType, OutVectorSpaceDual> dual_out_vector
                = to_vector_space<OutVectorSpaceDual>(mapping, coord, in_vector);
        return to_vector_space<OutVectorSpace>(mapping, coord, dual_out_vector);

    } else {
        using ArgBasis = ddc::to_type_seq_t<typename Mapping::CoordArg>;
        using ResultBasis = ddc::to_type_seq_t<typename Mapping::CoordResult>;
        using InContraVectorSpace = get_contravariant_dims_t<InVectorSpace>;
        using OutContraVectorSpace = get_contravariant_dims_t<OutVectorSpace>;
        if constexpr ((std::is_same_v<InContraVectorSpace, ArgBasis>)&&(
                              std::is_same_v<OutContraVectorSpace, ResultBasis>)) {
            if constexpr ((is_contravariant_vector_index_set_v<InVectorSpace>)&&(
                                  is_contravariant_vector_index_set_v<OutVectorSpace>)) {
                // A_{\{p\}}^i = J_{\{q\rightarrow p\}}^i_j A_{\{q\}}^j
                return tensor_mul(
                        index<'i', 'j'>(mapping.jacobian_matrix(coord)),
                        index<'j'>(in_vector));
            } else {
                // A_{\{p\} i} = (J_{\{q\rightarrow p\}}^{-T})_i^j A_{\{q\} j}
                //             = (J_{\{q\rightarrow p\}}^{-1})_j^i A_{\{q\} j}
                InverseJacobianMatrix inv_jacobian(mapping);
                return tensor_mul(index<'j', 'i'>(inv_jacobian(coord)), index<'j'>(in_vector));
            }
        } else {
            static_assert((std::is_same_v<InContraVectorSpace, ResultBasis>)&&(
                    std::is_same_v<OutContraVectorSpace, ArgBasis>));
            if constexpr ((is_contravariant_vector_index_set_v<InVectorSpace>)&&(
                                  is_contravariant_vector_index_set_v<OutVectorSpace>)) {
                static_assert(is_contravariant_vector_index_set_v<OutVectorSpace>);
                // A_{\{q\}}^i = (J_{\{q\rightarrow p\}}^{-1})^i_j A_{\{q\}}^j
                InverseJacobianMatrix inv_jacobian(mapping);
                DTensor<OutVectorSpace, get_covariant_dims_t<InVectorSpace>> I_J
                        = inv_jacobian(coord);
                return tensor_mul(index<'i', 'j'>(I_J), index<'j'>(in_vector));
            } else {
                // A_{\{p\} i} = (J_{\{q\rightarrow p\}}^{-T})_i^j A_{\{q\} j}
                //             = (J_{\{q\rightarrow p\}}^{-1})_j^i A_{\{q\} j}
                //             = (J_{\{p\rightarrow q\}})_j^i A_{\{q\} j}
                return tensor_mul(
                        index<'j', 'i'>(mapping.jacobian_matrix(coord)),
                        index<'j'>(in_vector));
            }
        }
    }
}

/**
 * @brief A helper method to get a vector field on a different vector space.
 *
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] vector_field The vector field to be mapped to the new vector space.
 * @param[in] mapping A mapping describing the relation between the 2 vector
 *              spaces (or describing the relation between the vector space and
 *              the Cartesian space if a change of variance is required.
 *
 * @returns A VectorField or VectorFieldMem containing the vectors in the
 *      requested vector space.
 */
template <
        class OutVectorSpace,
        class ExecSpace,
        class Mapping,
        class ElementType,
        class IdxRangeType,
        class InVectorSpace,
        class LayoutStridedPolicy>
void copy_to_vector_space(
        ExecSpace exec_space,
        VectorField<ElementType, IdxRangeType, OutVectorSpace, typename ExecSpace::memory_space>
                vector_field_out,
        Mapping mapping,
        VectorConstField<
                ElementType,
                IdxRangeType,
                InVectorSpace,
                typename ExecSpace::memory_space,
                LayoutStridedPolicy> vector_field)
{
    ddc::parallel_for_each(
            exec_space,
            get_idx_range(vector_field),
            KOKKOS_LAMBDA(IdxType idx) {
                IdxJacobianType coord_idx(idx);
                CoordType coord = ddc::coordinate(coord_idx);
                ddcHelper::assign_vector_field_element(
                        vector_field_out,
                        idx,
                        to_vector_space<OutVectorSpace>(mapping, coord, vector_field(idx)));
            });
}

/**
 * @brief A helper method to get a vector field on a different vector space.
 * If the requested vector space is the same as the current vector space
 * then the same vector field is returned.
 * If the vector space is different then the vectors in the vector
 * field are mapped to the new vector space and a VectorFieldMem is returned.
 *
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] vector_field The vector field to be mapped to the new vector space.
 * @param[in] mapping A mapping describing the relation between the 2 vector
 *              spaces (or describing the relation between the vector space and
 *              the Cartesian space if a change of variance is required.
 *
 * @tparam OutVectorSpace The vector space of the final vector.
 *
 * @returns A VectorField or VectorFieldMem containing the vectors in the
 *      requested vector space.
 */
template <
        class OutVectorSpace,
        class ExecSpace,
        class Mapping,
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
        Mapping mapping)
{
    if constexpr (std::is_same_v<InVectorSpace, OutVectorSpace>) {
        return vector_field;
    } else {
        using IdxType = typename IdxRangeType::discrete_element_type;
        using CoordType = typename Mapping::CoordJacobian;
        using IdxJacobianType = find_idx_t<CoordType, IdxRangeType>;
        VectorFieldMem<
                std::remove_const_t<ElementType>,
                IdxRangeType,
                OutVectorSpace,
                typename ExecSpace::memory_space>
                vector_field_out(get_idx_range(vector_field));
        copy_to_vector_space(
                exec_space,
                get_field(vector_field_out),
                get_const_field(vector_field),
                mapping);
        return vector_field_out;
    }
}
