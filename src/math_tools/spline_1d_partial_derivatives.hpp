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
              typename Spline1DBuilder::batched_interpolation_domain_type,
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
            typename Spline1DBuilder::batched_interpolation_domain_type,
            typename Spline1DBuilder::continuous_dimension_type>;

public:
    /// The type of the object that will be differentiated.
    using typename base_type::DFieldMemType;
    using typename base_type::DFieldType;

    /// The type of the calculated derivative.
    using typename base_type::DConstFieldType;

private:
    // Type for spline representation of the field
    using IdxRangeBS = typename Spline1DBuilder::batched_spline_domain_type;
    using DFieldBSMem = DFieldMem<IdxRangeBS>;
    using DFieldBS = DField<IdxRangeBS>;

    Spline1DBuilder const& m_builder;
    Spline1DEvaluator const& m_evaluator;

public:
    /**
    * @brief Construct an instance of the class Spline1DPartialDerivative.
    *
    * @param field_builder Builder for intermediate interpolation representation.
    * @param field_evaluator Evaluator for intermediate interpolation representation.
    */
    explicit Spline1DPartialDerivative(
            Spline1DBuilder const& builder,
            Spline1DEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[out] dfieldval_dxi Partial derivatives in Xi direction.
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    DFieldType operator()(DConstFieldType fieldval) const final
    {
        // Build spline representation of the field ....................................
        DFieldBSMem spline_coefs_alloc(m_builder.batched_spline_domain());
        DFieldBS spline_coefs = get_field(spline_coefs_alloc);

        m_builder(spline_coefs, fieldval);
        DFieldMemType differentiated_fieldval(get_idx_range(fieldval));
        m_evaluator.deriv(get_field(differentiated_fieldval), get_const_field(spline_coefs));
        return get_field(differentiated_fieldval);
    }
};


/**
 * @brief A class which stores information necessary to create a pointer to an instance of the class.Spline1DPartialDerivative
 *
 * This class allows an instance of the Spline1DPartialDerivative class to be instantiated where necessary. This allows the
 * memory allocated in the private members of the Spline1DPartialDerivative to be freed when the object is not in use.
 * These objects are: m_coefs.
 */
template <class Spline1DBuilder, class Spline1DEvaluator>
class Spline1DPartialDerivativeCreator
    : public IPartialDerivativeCreator<
              typename Spline1DBuilder::batched_interpolation_domain_type,
              typename Spline1DBuilder::continuous_dimension_type>
{
private:
    Spline1DBuilder const& m_builder;

    Spline1DEvaluator const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating Spline1DPartialDerivative objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    Spline1DPartialDerivativeCreator(
            Spline1DBuilder const& builder,
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
    std::unique_ptr<IPartialDerivative<
            typename Spline1DBuilder::batched_interpolation_domain_type,
            typename Spline1DBuilder::continuous_dimension_type>>
    preallocate() const
    {
        return std::make_unique<Spline1DPartialDerivative<
                Spline1DBuilder,
                Spline1DEvaluator>>(m_builder, m_evaluator);
    }
};
