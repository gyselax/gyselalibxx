// SPDX-License-Identifier: MIT
#pragma once

#include <iostream>

#include "ddc/discrete_domain.hpp"
#include "ddc/uniform_point_sampling.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivatives.hpp"


/**
 * @brief A class which implements a partial derivative operator
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
                        double const h1 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
                        double const h2
                                = ddc::coordinate(ix + 2 * step) - ddc::coordinate(ix + step);
                        double const c3 = -h1 / (h2 * (h1 + h2));
                        double const c2 = 1. / h1 + 1. / h2;
                        double const c1 = -c3 - c2;
                        dfieldval_dxi(ibx) = c1 * fieldval(ibx) + c2 * fieldval(ib, ix + step)
                                             + c3 * fieldval(ib, ix + 2 * step);
                    } else if (ix == idxrange_deriv.back()) {
                        double const h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
                        double const h2
                                = ddc::coordinate(ix - step) - ddc::coordinate(ix - 2 * step);
                        double const c3 = h1 / (h2 * (h1 + h2));
                        double const c2 = -(h1 + h2) / (h1 * h2);
                        double const c1 = -c3 - c2;
                        dfieldval_dxi(ibx) = c1 * fieldval(ibx) + c2 * fieldval(ib, ix - step)
                                             + c3 * fieldval(ib, ix - 2 * step);
                    } else {
                        // forward FDM
                        double const forward_fdm
                                = (fieldval(ib, ix + step) - fieldval(ibx))
                                  / (ddc::coordinate(ix + step) - ddc::coordinate(ix));
                        // backward FDM
                        double const backward_fdm
                                = (fieldval(ibx) - fieldval(ib, ix - step))
                                  / (ddc::coordinate(ix) - ddc::coordinate(ix - step));
                        // mean of the two
                        dfieldval_dxi(ibx) = (backward_fdm + forward_fdm) / 2.;
                    }
                });
    }
};
