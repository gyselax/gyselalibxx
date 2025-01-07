
// SPDX-License-Identifier: MIT
/**
 * @file partial_derivatives.hpp
 * File containing functions to compute the partial derivatives
 */

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"


/**
 * @brief A super class for a partial derivative operator
 */
template <class DFieldValue, class DerivativeDirection>
class IPartialDerivative
{
    static_assert(is_borrowed_chunk_v<DFieldValue>);

public:
    /// The dimension Xi on which the partial derivative is calculated.
    using DerivativeDirection = typename FieldXiBuilderBatched::continuous_dimension_type;

    /// The index range on which this operator acts.
    using IdxRangeFieldVal = typename DFieldValue::discrete_domain_type;

    /// The type of the object that will be differentiated.
    using DFieldVal = DFieldValue;

    /// The type of the calculated derivative.
    using DConstFieldVal = typename DFieldValue::view_type;

public:
    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[out] dfieldval_dxi Partial derivatives in Xi direction.
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    void operator()(DFieldVal dfieldval_dxi, DConstFieldVal fieldval) = 0;
};
