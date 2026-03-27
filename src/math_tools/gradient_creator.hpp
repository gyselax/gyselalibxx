// SPDX-License-Identifier: MIT
#pragma once
#include <tuple>

#include "coord_transformation_tools.hpp"
#include "ipartial_derivative.hpp"
#include "vector_mapper.hpp"

/**
 * @brief Operator to calculate the gradient of a function.
 *
 * @tparam IdxRangeFull The index range on which the function being differentiated
 * should be defined.
 * @tparam DerivativeDims The dimensions along which the function should be
 * differentiated to construct the gradient.
 */
template <typename IdxRangeFull, typename... DerivativeDims>
class GradientCreator
{
    static_assert((DerivativeDims::IS_CONTRAVARIANT && ...));

private:
    std::tuple<IPartialDerivativeCreator<IdxRangeFull, DerivativeDims> const&...>
            m_derivative_creators;

public:
    /**
     * @brief Create a GradientCreator.
     * @param[in] partial_derivative_operator The operators used to calculate the
     * derivative along each of the dimensions.
     */
    explicit GradientCreator(IPartialDerivativeCreator<
                             IdxRangeFull,
                             DerivativeDims> const&... partial_derivative_operator)
        : m_derivative_creators(std::tie(partial_derivative_operator...))
    {
    }

    /**
     * @brief Fill a VectorField with the values of the gradient of a function.
     *
     * @param[out] grad_func_cov The gradient of the function expressed as a covariant vector.
     * @param[in] func The function whose gradient is calculated.
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
