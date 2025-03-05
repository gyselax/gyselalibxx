// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolator_2d.hpp"

/**
 * @brief A class for interpolating a function using splines in polar coordinates.
 *
 * @tparam RadialExtrapolationRule The extrapolation rule applied at the outer radial bound.
 */
template <class Spline2DBuilder, class Spline2DEvaluator>
class SplineInterpolator2D
    : public IInterpolator2D<
              typename Spline2DBuilder::interpolation_domain_type,
              typename Spline2DBuilder::batched_interpolation_domain_type>
{
    using base_type = IInterpolator2D<
            typename Spline2DBuilder::interpolation_domain_type,
            typename Spline2DBuilder::batched_interpolation_domain_type>;

public:
    using typename base_type::CConstFieldType;
    using typename base_type::CoordType;
    using typename base_type::DFieldType;

private:
    Spline2DBuilder const& m_builder;

    Spline2DEvaluator const& m_evaluator;

    using IdxRangeBSRTheta = typename Spline2DBuilder::batched_spline_domain_type;

    mutable DFieldMem<IdxRangeBSRTheta> m_coefs;

    using r_deriv_type = DConstField<typename Spline2DBuilder::batched_derivs_domain_type1>;
    using theta_deriv_type = DConstField<typename Spline2DBuilder::batched_derivs_domain_type2>;
    using mixed_deriv_type = DConstField<typename Spline2DBuilder::batched_derivs_domain_type>;

public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    SplineInterpolator2D(Spline2DBuilder const& builder, Spline2DEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(get_spline_idx_range(builder))
    {
    }

    /**
     * @brief Approximate the value of a function at a set of polar coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data
     * 					On input: an array containing the value of the function at the interpolation points.
     * 					On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates
     * 					The polar coordinates where the function should be evaluated.
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    DFieldType operator()(DFieldType const inout_data, CConstFieldType const coordinates)
            const override
    {
        m_builder(get_field(m_coefs), get_const_field(inout_data));
        m_evaluator(get_field(inout_data), coordinates, get_const_field(m_coefs));

        return inout_data;
    }
};



/**
 * @brief A class which stores information necessary to create a pointer to an instance of the SplineInterpolator2D class.
 *
 * This class allows an instance of the SplineInterpolator2D class where necessary. This allows the
 * memory allocated in the private members of the SplineInterpolator2D to be freed when the object is not in use.
 * These objects are: m_coefs.
 */
template <class Spline2DBuilder, class Spline2DEvaluator>
class PreallocatableSplineInterpolator2D
    : public IPreallocatableInterpolator2D<
              typename Spline2DBuilder::interpolation_domain_type,
              typename Spline2DBuilder::batched_interpolation_domain_type>
{
private:
    Spline2DBuilder const& m_builder;

    Spline2DEvaluator const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolator2D objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSplineInterpolator2D(
            Spline2DBuilder const& builder,
            Spline2DEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    /**
     * Create a pointer to an instance of the SplineInterpolator2D class.
     *
     * @return A pointer to an instance of the SplineInterpolator2D class.
     */
    std::unique_ptr<IInterpolator2D<
            typename Spline2DBuilder::interpolation_domain_type,
            typename Spline2DBuilder::batched_interpolation_domain_type>>
    preallocate() const override
    {
        return std::make_unique<
                SplineInterpolator2D<Spline2DBuilder, Spline2DEvaluator>>(m_builder, m_evaluator);
    }
};
