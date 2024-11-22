// SPDX-License-Identifier: MIT
#pragma once

#include <sll/view.hpp>

#include "directional_tag.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

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
 * @tparam InVectorSpace A NDTag<XIn, YIn> describing the dimensions of the coordinate system taken as input.
 * @tparam OutVectorSpace A NDTag<XOut, YOut> describing the dimensions of the coordinate system returned as output.
 * @tparam Mapping A class describing a mapping system.
 * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
 */
template <class XIn, class YIn, class XOut, class YOut, class Mapping, class ExecSpace>
class VectorMapper<NDTag<XIn, YIn>, NDTag<XOut, YOut>, Mapping, ExecSpace>
{
    static_assert(is_accessible_v<ExecSpace, Mapping>);

public:
    /// @brief The type of the memory space where the field is saved (CPU vs GPU).
    using memory_space = typename ExecSpace::memory_space;

    /// The vector type in the coordinate system taken as input.
    using vector_element_type_in = typename ddc::detail::TaggedVector<double, XIn, YIn>;
    /// The vector type in the coordinate system returned as output.
    using vector_element_type_out = typename ddc::detail::TaggedVector<double, XOut, YOut>;

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
            VectorField<double, IdxRangeType, NDTag<XOut, YOut>, memory_space, LayoutStridedPolicy1>
                    vector_field_output,
            VectorConstField<
                    double,
                    IdxRangeType,
                    NDTag<XIn, YIn>,
                    memory_space,
                    LayoutStridedPolicy2> vector_field_input)
    {
        using IdxType = typename IdxRangeType::discrete_element_type;

        Mapping mapping_proxy = m_mapping;

        ddc::parallel_for_each(
                exec_space,
                get_idx_range(vector_field_input),
                KOKKOS_LAMBDA(IdxType idx) {
                    Matrix_2x2 map_J;
                    mapping_proxy.jacobian_matrix(ddc::coordinate(idx), map_J);

                    vector_element_type_out vector_out;
                    vector_out.array() = mat_vec_mul(map_J, vector_field_input(idx).array());
                    ddcHelper::get<XOut>(vector_field_output)(idx) = ddc::get<XOut>(vector_out);
                    ddcHelper::get<YOut>(vector_field_output)(idx) = ddc::get<YOut>(vector_out);
                });
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
        class PseudoCartToCartMapping,
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
                NDTag<X, Y>,
                typename ExecSpace::memory_space,
                LayoutStridedPolicy> vector_field,
        PseudoCartToCartMapping mapping)
{
    static_assert(
            std::is_same_v<X, typename PseudoCartToCartMapping::cartesian_tag_x>,
            "The mapping must describe the input vector space.");
    static_assert(
            std::is_same_v<Y, typename PseudoCartToCartMapping::cartesian_tag_y>,
            "The mapping must describe the input vector space.");
    using X_pC = typename PseudoCartToCartMapping::pseudo_cartesian_tag_x;
    using Y_pC = typename PseudoCartToCartMapping::pseudo_cartesian_tag_y;
    if constexpr (std::is_same_v<X, X_pC> && std::is_same_v<Y, Y_pC>) {
        return vector_field;
    } else {
        VectorFieldMem<
                std::remove_const_t<ElementType>,
                IdxRangeType,
                NDTag<X_pC, Y_pC>,
                typename ExecSpace::memory_space>
                vector_field_pC(get_idx_range(vector_field));
        VectorMapper<NDTag<X, Y>, NDTag<X_pC, Y_pC>, PseudoCartToCartMapping, ExecSpace>
                vector_mapping(mapping);
        vector_mapping(exec_space, get_field(vector_field_pC), get_const_field(vector_field));
        return vector_field_pC;
    }
}
