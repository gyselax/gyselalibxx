// SPDX-License-Identifier: MIT
#pragma once

#include "inverse_jacobian_matrix.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

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

/** The general predeclaration of VectorMapper.
 * @see @ref VectorMapperImplementation
 */
template <class InVectorSpace, class OutVectorSpace, class Mapping, class ExecSpace>
class VectorMapper;

/**
 * @brief A class to map vector fields from one coordinate system to another
 *
 * @anchor VectorMapperImplementation
 *
 * @tparam InVectorSpace A VectorIndexSet<XIn, YIn> describing the dimensions of the coordinate system taken as input.
 * @tparam OutVectorSpace A VectorIndexSet<XOut, YOut> describing the dimensions of the coordinate system returned as output.
 * @tparam Mapping A class describing a mapping system.
 * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
 */
template <class XIn, class YIn, class XOut, class YOut, class Mapping, class ExecSpace>
class VectorMapper<VectorIndexSet<XIn, YIn>, VectorIndexSet<XOut, YOut>, Mapping, ExecSpace>
{
    static_assert(is_accessible_v<ExecSpace, Mapping>);
    static_assert(
            (std::is_same_v<Coord<XIn, YIn>, typename Mapping::CoordArg>)
            || (std::is_same_v<Coord<XIn, YIn>, typename Mapping::CoordResult>));
    static_assert(
            (std::is_same_v<Coord<XOut, YOut>, typename Mapping::CoordArg>)
            || (std::is_same_v<Coord<XOut, YOut>, typename Mapping::CoordResult>));

public:
    /// @brief The type of the memory space where the field is saved (CPU vs GPU).
    using memory_space = typename ExecSpace::memory_space;

private:
    Mapping m_mapping;

public:
    /**
     * @brief A constructor for the VectorMapper.
     * @param[in] mapping The mapping description.
     */
    explicit VectorMapper(Mapping mapping) : m_mapping(mapping) {}

    /**
     * @brief Convert vectors defined in the input coordinate system to equivalent vectors in the
     * output coordinate system.
     *
     * @param[in] exec_space The space on which the function is executed (CPU/GPU).
     * @param[out] vector_field_output The vector field containing the vectors in the output coordinate system.
     * @param[out] vector_field_input The vector field containing the vectors in the input coordinate system.
     */
    template <class IdxRangeType, class LayoutStridedPolicy1, class LayoutStridedPolicy2>
    void operator()(
            ExecSpace exec_space,
            VectorField<
                    double,
                    IdxRangeType,
                    VectorIndexSet<XOut, YOut>,
                    memory_space,
                    LayoutStridedPolicy1> vector_field_output,
            VectorConstField<
                    double,
                    IdxRangeType,
                    VectorIndexSet<XIn, YIn>,
                    memory_space,
                    LayoutStridedPolicy2> vector_field_input)
    {
        using IdxType = typename IdxRangeType::discrete_element_type;

        if constexpr (std::is_same_v<Coord<XIn, YIn>, typename Mapping::CoordArg>) {
            Mapping mapping_proxy = m_mapping;
            ddc::parallel_for_each(
                    exec_space,
                    get_idx_range(vector_field_input),
                    KOKKOS_LAMBDA(IdxType idx) {
                        Tensor jacobian = mapping_proxy.jacobian_matrix(ddc::coordinate(idx));
                        DVector<XOut, YOut> vector_field_out = tensor_mul(
                                index<'i', 'j'>(jacobian),
                                index<'j'>(vector_field_input));
                        ddcHelper::get<XOut>(vector_field_output)(idx)
                                = ddcHelper::get<XOut>(vector_field_out);
                        ddcHelper::get<YOut>(vector_field_output)(idx)
                                = ddcHelper::get<YOut>(vector_field_out);
                    });
        } else {
            InverseJacobianMatrix<Mapping> inv_mapping(m_mapping);
            ddc::parallel_for_each(
                    exec_space,
                    get_idx_range(vector_field_input),
                    KOKKOS_LAMBDA(IdxType idx) {
                        Tensor map_J = inv_mapping(ddc::coordinate(idx));

                        DVector<XOut, YOut> vector_out = tensor_mul(
                                index<'i', 'j'>(map_J),
                                index<'j'>(vector_field_input(idx)));
                        ddcHelper::get<XOut>(vector_field_output)(idx)
                                = ddcHelper::get<XOut>(vector_out);
                        ddcHelper::get<YOut>(vector_field_output)(idx)
                                = ddcHelper::get<YOut>(vector_out);
                    });
        }
    }
};


/**
 * @brief A helper class to get a vector field on a pseudo-Cartesian geometry.
 * If the pseudo-Cartesian geometry is the same as the Cartesian geometry then
 * the same vector field is returned.
 * If the pseudo-Cartesian geometry is different then the vectors in the vector
 * field are mapped to the new geometry and a VectorFieldMem is returned.
 *
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] vector_field The vector field to be mapped to the pseudo-Cartesian
 *      geometry.
 * @param[in] mapping A mapping describing the relation between the Cartesian
 *      and pseudo-Cartesian geometries.
 *
 * @returns A VectorField or VectorFieldMem containing the vectors in the
 *      pseudo-Cartesian geometry.
*/
template <
        class ExecSpace,
        class Mapping,
        class ElementType,
        class IdxRangeType,
        class X,
        class Y,
        class LayoutStridedPolicy>
auto create_geometry_mirror_view(
        ExecSpace exec_space,
        VectorField<
                ElementType,
                IdxRangeType,
                VectorIndexSet<X, Y>,
                typename ExecSpace::memory_space,
                LayoutStridedPolicy> vector_field,
        Mapping mapping)
{
    using CoordOut = std::conditional_t<
            std::is_same_v<typename Mapping::CoordArg, Coord<X, Y>>,
            typename Mapping::CoordResult,
            typename Mapping::CoordArg>;
    if constexpr (std::is_same_v<CoordOut, Coord<X, Y>>) {
        return vector_field;
    } else {
        using X_out = ddc::type_seq_element_t<0, ddc::to_type_seq_t<CoordOut>>;
        using Y_out = ddc::type_seq_element_t<1, ddc::to_type_seq_t<CoordOut>>;
        VectorFieldMem<
                std::remove_const_t<ElementType>,
                IdxRangeType,
                VectorIndexSet<X_out, Y_out>,
                typename ExecSpace::memory_space>
                vector_field_out(get_idx_range(vector_field));
        VectorMapper<VectorIndexSet<X, Y>, VectorIndexSet<X_out, Y_out>, Mapping, ExecSpace>
                vector_mapping(mapping);
        vector_mapping(exec_space, get_field(vector_field_out), get_const_field(vector_field));
        return vector_field_out;
    }
}
