// SPDX-License-Identifier: MIT
#pragma once

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

        // Calculate forward differences at left boundary
        double inv_dx0 = 1.0
                         / (ddc::coordinate(idxrange_deriv.front() + step)
                            - ddc::coordinate(idxrange_deriv.front()));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange_batch,
                KOKKOS_LAMBDA(IdxBatch ib) {
                    dfieldval_dxi(idxrange_deriv.front(), ib)
                            = (fieldval(idxrange_deriv.front() + step, ib)
                               - fieldval(idxrange_deriv.front(), ib))
                              * inv_dx0;
                });

        // Calculate backward differences at right boundary
        double inv_dxN = 1.0
                         / (ddc::coordinate(idxrange_deriv.back())
                            - ddc::coordinate(idxrange_deriv.back() - step));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange_batch,
                KOKKOS_LAMBDA(IdxBatch ib) {
                    dfieldval_dxi(idxrange_deriv.front(), ib)
                            = (fieldval(idxrange_deriv.back(), ib)
                               - fieldval(idxrange_deriv.back() - step, ib))
                              * inv_dxN;
                });

        IdxRangeDeriv idxrange_deriv_central = idxrange_deriv.remove(step, step);
        IdxRangeFieldVal idxrange(idxrange_deriv_central, idxrange_batch);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange,
                KOKKOS_LAMBDA(IdxFieldVal ibx) {
                    IdxBatch ib(ibx);
                    IdxDeriv ix(ibx);
                    dfieldval_dxi(ibx)
                            = (fieldval(ib, ix + step) - fieldval(ib, ix - step))
                              / (ddc::coordinate(ix + step) - ddc::coordinate(ix - step));
                });
    }
};
