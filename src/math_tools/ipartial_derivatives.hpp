// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "type_seq_tools.hpp"


/**
 * @brief A super class for a partial derivative operator.
 * @tparam DFieldValue The type of the field on which the operator acts.
 * @tparam DerivDirection The dimension Xi on which the partial derivative is calculated.
 */
template <class DFieldValue, class DerivDirection>
class IPartialDerivative
{
    static_assert(ddc::is_borrowed_chunk_v<DFieldValue>);

public:
    /// The dimension Xi on which the partial derivative is calculated.
    using DerivativeDirection = DerivDirection;

    /// The index range on which this operator acts.
    using IdxRangeFieldVal = typename DFieldValue::discrete_domain_type;

    /// The type of the object that will be differentiated.
    using DFieldVal = DFieldValue;

    /// The type of the calculated derivative.
    using DConstFieldVal = typename DFieldValue::view_type;

    /// The type of the grid on the dimension Xi on which the partial derivative is calculated.
    using GridDerivativeDirection
            = find_grid_t<DerivativeDirection, ddc::to_type_seq_t<IdxRangeFieldVal>>;

    /// The index range of the dimension Xi on which the partial derivative is calculated.
    using IdxRangeDeriv = IdxRange<GridDerivativeDirection>;

    /// The index range of all dimensions except Xi.
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFieldVal, GridDerivativeDirection>;

public:
    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[out] dfieldval_dxi Partial derivatives in Xi direction.
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    virtual void operator()(DFieldVal dfieldval_dxi, DConstFieldVal fieldval) const = 0;
};
