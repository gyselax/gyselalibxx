// SPDX-License-Identifier: MIT
#pragma once
#include <tuple>

#include "coord_transformation_tools.hpp"
#include "ipartial_derivative.hpp"
#include "vector_mapper.hpp"

template <
        concepts::MappingWithJacobian MappingType,
        typename IdxRangeFull,
        typename... DerivativeDims>
class GradientCreator
{
    static_assert((DerivativeDims::IS_CONTRAVARIANT && ...));
    static_assert(is_accessible_v<Kokkos::DefaultExecutionSpace, MappingType>);

private:
    MappingType const& m_mapping;
    std::tuple<IPartialDerivativeCreator<IdxRangeFull, DerivativeDims> const&...>
            m_derivative_creators;

public:
    GradientCreator(
            MappingType const& mapping,
            IPartialDerivativeCreator<
                    IdxRangeFull,
                    DerivativeDims> const&... partial_derivative_operator)
        : m_mapping(mapping)
        , m_derivative_creators(std::tie(partial_derivative_operator...))
    {
    }

    /**
     * @brief Fill a VectorField with the values of the gradient of a function.
     *
     * @param[in] field A field to be passed to the constructor of IPartialDerivative.
     *
     * @return A pointer to an IPartialDerivative object.
     *
     * @see IPartialDerivative
     */
    void operator()(
            DVectorField<IdxRangeFull, get_covariant_dims_t<VectorIndexSet<DerivativeDims...>>>
                    grad_func_cov,
            DConstField<IdxRangeFull> func) const
    {
        (((*std::get<IPartialDerivativeCreator<IdxRangeFull, DerivativeDims> const&>(
                    m_derivative_creators)
                    .create_instance(get_const_field(func)))(
                 ddcHelper::get<typename DerivativeDims::Dual>(grad_func_cov))),
         ...);
    }
};
