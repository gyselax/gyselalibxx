// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"


/**
 * @brief A class which implements a gradient operator 
 *
 * @tparam Spline1DBuilder A 1D spline builder.
 * @tparam Spline1DEvaluator A 1D spline evaluator.
 */
class Gradient
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

    using typename base_type::DConstFieldType;
    using typename base_type::DFieldType;

    using IdxRangeBS = typename Spline1DBuilder::batched_spline_domain_type;
    using DFieldBSMem = DFieldMem<IdxRangeBS>;
    using DFieldBS = DField<IdxRangeBS>;

    Spline1DBuilder const& m_builder;
    Spline1DEvaluator const& m_evaluator;
    DFieldBSMem m_spline_coefs;

public:
    /**
     * @brief Construct an instance of the class Spline1DPartialDerivative.
     *
     * @param builder A 1D spline builder.
     * @param evaluator A 1D spline evaluator.
     * @param field The field to be differentiated.
     */
    explicit Spline1DPartialDerivative(
            Spline1DBuilder const& builder,
            Spline1DEvaluator const& evaluator,
            DConstFieldType const field)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_spline_coefs(builder.batched_spline_domain())
    {
        m_builder(get_field(m_spline_coefs), field);
    }

    /**
     * @brief Compute the partial derivative of a field in the direction 
     * where the field is represented using 1d splines.
     *
     * @param[out] differentiated_field Contains on output the value of the differentiated field.
     */
    void operator()(DFieldType differentiated_field) const final
    {
        m_evaluator.deriv(differentiated_field, get_const_field(m_spline_coefs));
    }
};
