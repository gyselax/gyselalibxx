// SPDX-License-Identifier: MIT
#pragma once

#include <iostream>

#include "ddc/discrete_domain.hpp"
#include "ddc/uniform_point_sampling.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"


/**
 * @brief A class which implements a partial derivative operator
 * using a finite differences calculation of order two. A decentered
 * scheme is used at the boundary, whereas centred finite difference
 * are used inside the domain.
 */
template <class IdxRangeBatched, class DerivativeDimension>
class CentralFDMPartialDerivative : public IPartialDerivative<IdxRangeBatched, DerivativeDimension>
{
private:
    using base_type = IPartialDerivative<IdxRangeBatched, DerivativeDimension>;

public:
    /// The type of a reference to the field to be differentiated.
    using DFieldType = DField<IdxRangeBatched>;

    /// The type of the calculated derivative.
    using typename DFieldType::DConstFieldVal;

    /// The index range of the dimension Xi on which the partial derivative is calculated.
    using IdxRangeDeriv = base_type:: ddc::remove_dims_of_t<IdxRangeFieldVal, GridDerivativeDirection>;

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
    void operator()(DFieldType dfieldval_dxi, DConstFieldVal fieldval) const final
    {
        using IdxFieldVal = typename IdxRangeBatched::discrete_element_type;
        using IdxDeriv = typename IdxRangeDeriv::discrete_element_type;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;
        using IdxStepDeriv = typename IdxRangeDeriv::discrete_vector_type;

        IdxRangeBatched idxrange_full = get_idx_range(fieldval);
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

/**
 * @brief A class which stores information necessary to create a pointer to 
 * an instance of the CentralFDMPartialDerivative class.
 *
 * This class allows an instance of the CentralFDMPartialDerivative class to be instantiated where necessary.
 * Typically, the CentralFDMPartialDerivativeCreator is instantiated in the initialisation of the simulation, 
 * and the corresponding CentralFDMPartialDerivative object is instantiated where computing partial derivatives
 * is required. 
 * @tparam Spline1DBuilder A 1D spline builder.
 * @tparam Spline1DEvaluator A 1D spline evaluator.
 */
template <class DFieldTypeue, class DerivDirection>
class CentralFDMPartialDerivativeCreator
    : public IPartialDerivativeCreator<DFieldValue, DerivDirection>
{
private:
    using DConstFieldType = typename DFieldValue::view_type;

    Spline1DBuilder const& m_builder;
    Spline1DEvaluator const& m_evaluator;

public:
    /**
     * @brief Construct an instance of the CentralFDMPartialDerivativeCreator class.
     * @param[in] builder A 1d spline builder.
     * @param[in] evaluator A 1d spline evaluator.
     */
    CentralFDMPartialDerivativeCreator(
            Spline1DBuilder const& builder,
            Spline1DEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    /**
     * Create a pointer to an instance of the abstract class IPartialDerivative.
     * The type of the returned object will be determined when the pointer is 
     * dereferenced.
     *
     * @param[in] field A field to be differentiated.
     *
     * @return A pointer to an instance of the IPartialDerivative class.
     */
    std::unique_ptr<IPartialDerivative<
            typename Spline1DBuilder::batched_interpolation_domain_type,
            typename Spline1DBuilder::continuous_dimension_type>>
    create_instance(DConstFieldType field) const
    {
        return std::make_unique<CentralFDMPartialDerivative<
                Spline1DBuilder,
                Spline1DEvaluator>>(m_builder, m_evaluator, field);
    }
};
