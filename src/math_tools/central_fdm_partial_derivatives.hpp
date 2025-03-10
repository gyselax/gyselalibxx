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
 * using a finite differences calculation of order two. A decentered
 * scheme is used at the boundary, whereas centered finite difference
 * are used inside the domain.
 */
template <class DFieldValue, class DerivDirection>
class CentralFDMPartialDerivative : public IPartialDerivative<DFieldValue, DerivDirection>
{
private:
    using base_type = IPartialDerivative<DFieldValue, DerivDirection>;

public:
    /// The index range on which this operator acts.
    using typename base_type::IdxRangeFieldVal;

    /// The type of the object that will be differentiated.
    using typename base_type::DFieldVal;

    /// The type of the calculated derivative.
    using typename base_type::DConstFieldVal;

    /// The index range of the dimension Xi on which the partial derivative is calculated.
    using typename base_type::IdxRangeDeriv;

    /// The index range of all dimensions except Xi.
    using typename base_type::IdxRangeBatch;

public:
    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    * For more information about the coefficients, see `./README.md`
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

        IdxStepDeriv const step(1);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange_full,
                KOKKOS_LAMBDA(IdxFieldVal ibx) {
                    IdxBatch ib(ibx);
                    IdxDeriv ix(ibx);
                    double h1, h2, c1, c2, c3;
                    if (ix == idxrange_deriv.front()) {
                        h1 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
                        h2 = ddc::coordinate(ix + 2 * step) - ddc::coordinate(ix + step);
                        c3 = -h1 / (h2 * (h1 + h2));
                        c2 = 1. / h1 + 1. / h2;
                        c1 = -c3 - c2;
                        dfieldval_dxi(ibx) = c1 * fieldval(ibx) + c2 * fieldval(ib, ix + step)
                                             + c3 * fieldval(ib, ix + 2 * step);
                    } else if (ix == idxrange_deriv.back()) {
                        h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
                        h2 = ddc::coordinate(ix - step) - ddc::coordinate(ix - 2 * step);
                        c3 = h1 / (h2 * (h1 + h2));
                        c2 = -(h1 + h2) / (h1 * h2);
                        c1 = -c3 - c2;
                        dfieldval_dxi(ibx) = c1 * fieldval(ibx) + c2 * fieldval(ib, ix - step)
                                             + c3 * fieldval(ib, ix - 2 * step);
                    } else {
                        h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
                        h2 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
                        c3 = h1 / (h2 * (h1 + h2));
                        c2 = 1. / h1 - 1. / h2;
                        c1 = -c3 - c2;
                        dfieldval_dxi(ibx) = c1 * fieldval(ib, ix - step) + c2 * fieldval(ibx)
                                             + c3 * fieldval(ib, ix + step);
                    }
                });
    }
};
