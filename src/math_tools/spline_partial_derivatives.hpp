// SPDX-License-Identifier: MIT
/**
 * @file partial_derivatives.hpp
 * File containing functions to compute the partial derivatives
 */

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"


/**
 * @brief A class which implementes a partial derivative operator
 * using a spline interpolation.
 */
template <class FieldXiBuilderBatched, class FieldXiEvaluatorBatched>
class SplinePartialDerivative
    : IPartialDerivative<
              DField<typename FieldXiBuilderBatched::batched_interpolation_domain_type>,
              typename FieldXiBuilderBatched::continuous_dimension_type>
{
    static_assert(std::is_same_v<
                  typename FieldXiBuilderBatched::batched_spline_domain_type,
                  typename FieldXiEvaluatorBatched::batched_spline_domain_type>);
    static_assert(std::is_same_v<
                  typename FieldXiBuilderBatched::batched_interpolation_domain_type,
                  typename FieldXiEvaluatorBatched::batched_evaluation_domain_type>);

private:
    using base_type = IPartialDerivative<
            DField<typename FieldXiBuilderBatched::batched_interpolation_domain_type>,
            typename FieldXiBuilderBatched::continuous_dimension_type>;

public:
    /// The dimension Xi on which the partial derivative is calculated.
    using base_type::DerivativeDirection;

    /// The index range on which this operator acts.
    using base_type::IdxRangeFieldVal;

    /// The type of the object that will be differentiated.
    using base_type::DFieldVal;

    /// The type of the calculated derivative.
    using base_type::DConstFieldVal;

private:
    // Type for spline representation of the field
    using IdxRangeBSFieldXi = typename FieldXiBuilderBatched::batched_spline_domain_type;
    using FieldXiSplineMem = DFieldMem<IdxRangeBSFieldXi>;
    using FieldXiSplineCoeffs = DField<IdxRangeBSFieldXi>;

    FieldXiBuilderBatched const& m_fieldxi_builder;
    FieldXiEvaluatorBatched const& m_fieldxi_evaluator;

public:
    /**
    * @brief Construct an instance of the class SplinePartialDerivative.
    *
    * @param fieldxi_builder Builder for intermediate interpolation representation.
    * @param fieldxi_evaluator Evaluator for intermediate interpolation representation.
    */
    explicit SplinePartialDerivative(
            FieldXiBuilderBatched const& fieldxi_builder,
            FieldXiEvaluatorBatched const& fieldxi_evaluator)
        : m_fieldxi_builder(fieldxi_builder)
        , m_fieldxi_evaluator(fieldxi_evaluator)
    {
    }

    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[out] dfieldval_dxi Partial derivatives in Xi direction.
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    void operator()(DFieldVal dfieldval_dxi, DConstFieldVal fieldval) final
    {
        // Build spline representation of the field ....................................
        FieldXiSplineMem fieldxi_coefs_alloc(m_fieldxi_builder.batched_spline_domain());
        FieldXiSplineCoeffs fieldxi_coefs = get_field(fieldxi_coefs_alloc);

        m_fieldxi_builder(fieldxi_coefs, get_const_field(fieldval));
        m_fieldxi_evaluator.deriv(dfieldval_dxi, get_const_field(fieldxi_coefs));
    }
};
