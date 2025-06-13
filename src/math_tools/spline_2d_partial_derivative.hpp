// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"
#include "spline_builder_2d_cache.hpp"


/**
 * @brief A class which implements a partial derivative operator
 * using a 2d spline interpolation.
 * @tparam SplineBuilder2D A 2D spline builder.
 * @tparam SplineEvaluator2D A 2D spline evaluator.
 * @tparam DerivativeDimension The dimension along which we are differentiating.
 * @tparam IdxRangeBatched The type af the index range over which this operator will operate.
 */
template <
        class SplineBuilder2D,
        class SplineEvaluator2D,
        class DerivativeDimension,
        class IdxRangeBatched>
class Spline2DPartialDerivative : public IPartialDerivative<IdxRangeBatched, DerivativeDimension>
{
    static_assert(
            (std::is_same_v<
                    typename SplineEvaluator2D::continuous_dimension_type1,
                    DerivativeDimension>)
            || (std::is_same_v<
                    typename SplineEvaluator2D::continuous_dimension_type2,
                    DerivativeDimension>));

private:
    using base_type = IPartialDerivative<IdxRangeBatched, DerivativeDimension>;

    using typename base_type::DConstFieldType;
    using typename base_type::DFieldType;

    using IdxRangeBS = IdxRangeBatched;
    using DFieldBSMem = DFieldMem<IdxRangeBS>;
    using DFieldBS = DField<IdxRangeBS>;

    SplineBuilder2DCache<SplineBuilder2D, IdxRangeBatched>& m_builder_cache;
    SplineEvaluator2D const& m_evaluator;

public:
    /**
    * @brief Construct an instance of the class Spline2DPartialDerivative.
    *
    * @param builder_cache A 2D spline builder cache.
    * @param evaluator A 2D spline evaluator.
    * @param field The field to be differentiated.
    */
    explicit Spline2DPartialDerivative(
            SplineBuilder2DCache<SplineBuilder2D, IdxRangeBatched>& builder_cache,
            SplineEvaluator2D const& evaluator,
            DConstFieldType const field)
        : m_builder_cache(builder_cache)
        , m_evaluator(evaluator)
    {
        m_builder_cache.template compute_coeffs<DerivativeDimension>(field);
    }

    /**
    * @brief Compute the partial derivative of a field in the direction 
    * where the field is represented using 2d splines.
    *
    * @param[out] differentiated_field Contains on output the value of the differentiated field.
    */
    void operator()(DFieldType differentiated_field) const final
    {
        if constexpr (std::is_same_v<
                              DerivativeDimension,
                              typename SplineEvaluator2D::continuous_dimension_type1>) {
            m_evaluator.deriv_dim_1(differentiated_field, m_builder_cache());
        } else {
            m_evaluator.deriv_dim_2(differentiated_field, m_builder_cache());
        }
    }
};

/**
 * @brief A class which stores information necessary to create a pointer to 
 * an instance of the Spline2DPartialDerivative class.
 *
 * This class allows an instance of the Spline2DPartialDerivative class to be 
 * instantiated where necessary. Typically, the Spline2DPartialDerivativeCreator
 * is instantiated in the initialisation of the simulation, and the corresponding 
 * Spline2DPartialDerivative object is instantiated where computing partial 
 * derivatives is required. 
 *
 * @tparam SplineBuilder2D A 2D spline builder.
 * @tparam SplineEvaluator2D A 2D spline evaluator.
 * @tparam DerivativeDimension The dimension along which we are differentiating.
 * @tparam IdxRangeBatched The type af the index range over which this operator will operate.
 */
template <
        class SplineBuilder2D,
        class SplineEvaluator2D,
        class DerivativeDimension,
        class IdxRangeBatched>
class Spline2DPartialDerivativeCreator
    : public IPartialDerivativeCreator<IdxRangeBatched, DerivativeDimension>
{
private:
    using DConstFieldType = DConstField<IdxRangeBatched>;

    SplineBuilder2DCache<SplineBuilder2D, IdxRangeBatched>& m_builder_cache;
    SplineEvaluator2D const& m_evaluator;

public:
    /**
     * @brief Construct an instance of the Spline2DPartialDerivativeCreator class.
     * @param[in] builder_cache A 2d spline builder cache.
     * @param[in] evaluator A 2d spline evaluator.
     */
    Spline2DPartialDerivativeCreator(
            SplineBuilder2DCache<SplineBuilder2D, IdxRangeBatched>& builder_cache,
            SplineEvaluator2D const& evaluator)
        : m_builder_cache(builder_cache)
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
    std::unique_ptr<IPartialDerivative<IdxRangeBatched, DerivativeDimension>> create_instance(
            DConstFieldType field) const
    {
        return std::make_unique<Spline2DPartialDerivative<
                SplineBuilder2D,
                SplineEvaluator2D,
                DerivativeDimension,
                IdxRangeBatched>>(m_builder_cache, m_evaluator, field);
    }
};
