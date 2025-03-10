// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"


/**
 * @brief A class which implements a partial derivative operator
 * using a 1d spline interpolation.
 * @tparam SplineBuilder1D A 1D spline builder.
 * @tparam Spline1DEvaluator A 1D spline evaluator.
 */
template <class SplineBuilder1D, class Spline1DEvaluator>
class Spline1DPartialDerivative
    : public IPartialDerivative<
              typename SplineBuilder1D::batched_interpolation_domain_type,
              typename SplineBuilder1D::continuous_dimension_type>
{
    static_assert(std::is_same_v<
                  typename SplineBuilder1D::batched_spline_domain_type,
                  typename Spline1DEvaluator::batched_spline_domain_type>);
    static_assert(std::is_same_v<
                  typename SplineBuilder1D::batched_interpolation_domain_type,
                  typename Spline1DEvaluator::batched_evaluation_domain_type>);

private:
    using base_type = IPartialDerivative<
            typename SplineBuilder1D::batched_interpolation_domain_type,
            typename SplineBuilder1D::continuous_dimension_type>;

    using typename base_type::DConstFieldType;
    using typename base_type::DFieldType;

    using IdxRangeBS = typename SplineBuilder1D::batched_spline_domain_type;
    using DFieldBSMem = DFieldMem<IdxRangeBS>;
    using DFieldBS = DField<IdxRangeBS>;

    SplineBuilder1D const& m_builder;
    Spline1DEvaluator const& m_evaluator;
    DConstFieldType const m_field;
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
            SplineBuilder1D const& builder,
            Spline1DEvaluator const& evaluator,
            DConstFieldType const field)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_field(field)
        , m_spline_coefs(builder.batched_spline_domain())
    {
        m_builder(get_field(m_spline_coefs), m_field);
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


/**
 * @brief A class which stores information necessary to create a pointer to 
 * an instance of the Spline1DPartialDerivative class.
 *
 * This class allows an instance of the Spline1DPartialDerivative class to be instantiated where necessary.
 * Typically, the Spline1DPartialDerivativeCreator is instantiated in the initialisation of the simulation, 
 * and the corresponding Spline1DPartialDerivative object is instantiated where computing partial derivatives
 * is required. 
 * @tparam SplineBuilder1D A 1D spline builder.
 * @tparam Spline1DEvaluator A 1D spline evaluator.
 */
template <class SplineBuilder1D, class Spline1DEvaluator>
class Spline1DPartialDerivativeCreator
    : public IPartialDerivativeCreator<
              typename SplineBuilder1D::batched_interpolation_domain_type,
              typename SplineBuilder1D::continuous_dimension_type>
{
private:
    using DConstFieldType
            = DConstField<typename SplineBuilder1D::batched_interpolation_domain_type>;

    SplineBuilder1D const& m_builder;
    Spline1DEvaluator const& m_evaluator;

public:
    /**
     * @brief Construct an instance of the Spline1DPartialDerivativeCreator class.
     * @param[in] builder A 1d spline builder.
     * @param[in] evaluator A 1d spline evaluator.
     */
    Spline1DPartialDerivativeCreator(
            SplineBuilder1D const& builder,
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
            typename SplineBuilder1D::batched_interpolation_domain_type,
            typename SplineBuilder1D::continuous_dimension_type>>
    create_instance(DConstFieldType field) const
    {
        return std::make_unique<Spline1DPartialDerivative<
                SplineBuilder1D,
                Spline1DEvaluator>>(m_builder, m_evaluator, field);
    }
};
