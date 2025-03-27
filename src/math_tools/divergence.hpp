// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "mapping.hpp"
#include "vector_field.hpp"


/**
 * @brief A class which implements a gradient operator.
 * @tparam Mapping A mapping.
 */
template <class Mapping>
class Divergence
{
    static_assert(is_mapping_v<Mapping>);
  
    using CoordArg = typename Mapping::CoordArg;

    using Dims = ddc::to_type_seq_t<CoordArg>;
    using Dims_cov = vector_index_set_dual_t<Dims>;

    using DVectorType = typename MetricTensorType::ContravariantVectorType;
    using DVectorCov = typename MetricTensorType::CovariantVectorType;

    template <class IdxRange>
    using DVectorFieldType = DVectorField<IdxRange, Dims>;

    template <class IdxRange>
    using DVectorFieldCovType = DVectorField<IdxRange, Dims_cov>;

    template <class IdxRange>
    using DVectorConstFieldCovType = DVectorConstField<IdxRange, Dims_cov>;

    Mapping m_mapping;

public:
    /**
     * @brief Construct an instance of the class Divergence.
     *
     * @param mapping A Mapping.
     */
    explicit KOKKOS_FUNCTION Divergence(Mapping mapping)
        : m_mapping(mapping)
    {
    }

    /**
     * @brief Compute the gradient of a scalar field at a given coordinate, 
     * using partial derivatives of the field. The gradient is expressed 
     * on the contravariant basis.  
     *
     * @param[in] partial_derivatives A vector that contains the partial 
     * derivatives of the scalar field expressed at a given coordinate.
     * @param[in] coord The coordinate at which the gradient should be evaluated.
     *
     * @return The components of the gradient at a given coordinate, 
     * expressed on the contravariant basis.
     */
    KOKKOS_INLINE_FUNCTION DVectorType
    operator()(DVectorCov const& partial_derivatives, CoordArg const& coord) const
    {
        return tensor_mul(
                index<'i', 'j'>(m_metric_tensor.inverse(coord)),
                index<'j'>(partial_derivatives));
    }

    /**
     * @brief Compute the gradient of a scalar field using partial derivatives
     * of the field. The gradient is expressed on the contravariant basis.  
     *
     * @param[out] gradient A vector field that contains on output the 
     * value of the gradient components expressed on the contravariant
     * @param[in] partial_derivatives A vector field that contains the partial 
     * derivatives of the scalar field.
     * basis. 
     */
    template <class IdxRange>
    KOKKOS_FUNCTION void operator()(
            DVectorFieldType<IdxRange> const gradient,
            DVectorConstFieldCovType<IdxRange> const partial_derivatives) const
    {
        using IdxType = typename IdxRange::discrete_element_type;
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(gradient),
                KOKKOS_CLASS_LAMBDA(IdxType const idx) {
                    DVectorType gradient_proxy
                            = (*this)(partial_derivatives(idx), ddc::coordinate(idx));
                    ddcHelper::assign_vector_field_element(gradient, idx, gradient_proxy);
                });
    }
};
