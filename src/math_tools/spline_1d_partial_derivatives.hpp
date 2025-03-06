// SPDX-License-Identifier: MIT
/**
 * @file spline_1d_partial_derivatives_creator.hpp
 * File containing a classes to compute partial derivatives of a field using its 1d spline representation. 
 */

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivatives.hpp"


/**
 * @brief A class which implements a partial derivative operator
 * using a 1d spline interpolation.
 */
template <class Spline1DBuilder, class Spline1DEvaluator>
class Spline1DPartialDerivative
    : public IPartialDerivative<
              typename Spline1DBuilder::batched_interpolation_domain_type>,
              typename Spline1DBuilder::continuous_dimension_type>
{
    static_assert(std::is_same_v<
                  typename Spline1DBuilder::batched_spline_domain_type,
                  typename Spline1DEvaluator::batched_spline_domain_type>);
    static_assert(std::is_same_v<
                  typename Spline1DBuilder::batched_interpolation_domain_type,
                  typename Spline1DEvaluator::batched_evaluation_domain_type>);

private:
    using base_type = IPartialDerivative<
            DField<typename Spline1DBuilder::batched_interpolation_domain_type>,
            typename Spline1DBuilder::continuous_dimension_type>;

public:
    /// The dimension Xi on which the partial derivative is calculated.
    using typename base_type::DerivativeDirection;

    /// The index range on which this operator acts.
    using typename base_type::IdxRange;

    /// The type of the object that will be differentiated.
    using typename base_type::DFieldType;

    /// The type of the calculated derivative.
    using typename base_type::DConstFieldType;

private:
    // Type for spline representation of the field
    using IdxRangeBS = typename Spline1DBuilder::batched_spline_domain_type;
    using DFieldBSMem = DFieldMem<IdxRangeBS>;
    using DFieldBS = DField<IdxRangeBS>;

    Spline1DBuilder const& m_field_builder;
    Spline1DEvaluator const& m_field_evaluator;

public:
    /**
    * @brief Construct an instance of the class Spline1DPartialDerivative.
    *
    * @param field_builder Builder for intermediate interpolation representation.
    * @param field_evaluator Evaluator for intermediate interpolation representation.
    */
    explicit Spline1DPartialDerivative(
            Spline1DBuilder const& field_builder,
            Spline1DEvaluator const& field_evaluator)
        : m_field_builder(field_builder)
        , m_field_evaluator(field_evaluator)
        // TODO coeffs ?
    {
    }

    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[out] dfieldval_dxi Partial derivatives in Xi direction.
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    void operator()(DFieldVal dfieldval_dxi, DConstFieldVal fieldval) const final
    {
        // Build spline representation of the field ....................................
        FieldXiSplineMem fieldxi_coefs_alloc(m_fieldxi_builder.batched_spline_domain());
        FieldXiSplineCoeffs fieldxi_coefs = get_field(fieldxi_coefs_alloc);

        m_fieldxi_builder(fieldxi_coefs, get_const_field(fieldval));
        m_fieldxi_evaluator.deriv(dfieldval_dxi, get_const_field(fieldxi_coefs));
    }
};


/**
 * @brief A class which stores information necessary to create a pointer to an instance of the class.Spline1DPartialDerivative
 *
 * This class allows an instance of the Spline1DPartialDerivative class to be instantiated where necessary. This allows the
 * memory allocated in the private members of the Spline1DPartialDerivative to be freed when the object is not in use.
 * These objects are: m_coefs.
 */
template <class Spline1Dbuilder, class Spline1DEvaluator>
class PreallocatableSpline1DPartialDerivative
: public IPreallocatablePartialDerivative<
              typename Spline1DBuilder::interpolation_domain_type,
              typename Spline1DBuilder::batched_interpolation_domain_type>
{
private:
    Spline1Dbuilder const& m_builder;

    Spline1DEvaluator const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating Spline1DPartialDerivative objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSpline1DPartialDerivative(
            Spline1Dbuilder const& builder,
            Spline1DEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    /**
     * Create a pointer to an instance of the Spline1DPartialDerivative class.
     *
     * @return A pointer to an instance of the Spline1DPartialDerivative class.
     */
    std::unique_ptr<Spline1DPartialDerivative<
            typename Spline1Dbuilder::interpolation_domain_type,
            typename Spline1Dbuilder::batched_interpolation_domain_type>>
    preallocate() const
    {
        return std::make_unique<
                Spline1DPartialDerivative<Spline1Dbuilder, Spline1DEvaluator>>(m_builder, m_evaluator);
    }
};
