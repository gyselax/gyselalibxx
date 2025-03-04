// SPDX-License-Identifier: MIT
#pragma once

#include <iostream>

#include "ddc/uniform_point_sampling.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivatives.hpp"


/**
 * @brief A class which implementes a partial derivative operator
 * using a central finite differences calculation. Forward and
 * backward finite differences are used at the boundaries.
 */
template <class DFieldValue, class DerivDirection>
class CentralFDMPartialDerivative : public IPartialDerivative<DFieldValue, DerivDirection>
{
private:
    using base_type = IPartialDerivative<DFieldValue, DerivDirection>;

public:
    /// The dimension Xi on which the partial derivative is calculated.
    using typename base_type::DerivativeDirection;

    /// The index range on which this operator acts.
    using typename base_type::IdxRangeFieldVal;

    /// The type of the object that will be differentiated.
    using typename base_type::DFieldVal;

    /// The type of the calculated derivative.
    using typename base_type::DConstFieldVal;

    /// The type of the grid on the dimension Xi on which the partial derivative is calculated.
    using typename base_type::GridDerivativeDirection;

    /// The index range of the dimension Xi on which the partial derivative is calculated.
    using typename base_type::IdxRangeDeriv;

    /// The index range of all dimensions except Xi.
    using typename base_type::IdxRangeBatch;

public:
    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[out] dfieldval_dxi Partial derivatives in Xi direction.
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    void operator()(DFieldVal dfieldval_dxi, DConstFieldVal fieldval) const final
    {
        using IdxFieldVal = typename IdxRangeFieldVal::discrete_element_type;
        using IdxDeriv = typename IdxRangeDeriv::discrete_element_type;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;
        using IdxStepDeriv = typename IdxRangeDeriv::discrete_vector_type;

        IdxRangeFieldVal idxrange_full = get_idx_range(fieldval);
        IdxRangeDeriv idxrange_deriv(idxrange_full);
        IdxRangeBatch idxrange_batch(idxrange_full);

        IdxStepDeriv step(1);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange_full,
                KOKKOS_LAMBDA(IdxFieldVal ibx) {
                    IdxBatch ib(ibx);
                    IdxDeriv ix(ibx);
                    if (ix == idxrange_deriv.front()) {
                        // Calculate forward differences at left boundary
                        dfieldval_dxi(ibx) = (fieldval(ix + step, ib) - fieldval(ibx))
                                             / (ddc::coordinate(ix + step) - ddc::coordinate(ix));
                    } else if (ix == idxrange_deriv.back()) {
                        // Calculate forward differences at left boundary
                        dfieldval_dxi(ix, ib)
                                = (fieldval(ix, ib) - fieldval(ix - step, ib))
                                  / (ddc::coordinate(ix) - ddc::coordinate(ix-step));
                    } else {
                        dfieldval_dxi(ibx)
                                = (fieldval(ib, ix + step) - fieldval(ib, ix - step))
                                  / (ddc::coordinate(ix + step) - ddc::coordinate(ix - step));
                    }
                });
    }
};
