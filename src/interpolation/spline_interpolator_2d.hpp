// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolator_2d.hpp"

/**
 * @brief A class for interpolating a function using a 2D tensor product of splines.
 *
 * The class is parametrised by multiple template parameters. Please note that CTAD will deduce all these
 * template parameters from the Builder and Evaluator passed as constructor arguments.
 *
 * @tparam Spline2DBuilder The type of the 2D spline builder.
 * @tparam Spline2DEvaluator The type of the 2D spline evaluator.
 * @tparam IdxRangeBatched The type af the index range over which this operator will operate. This is
 *              necessary to define the internal spline representation.
 */
template <class Spline2DBuilder, class Spline2DEvaluator, class IdxRangeBatched>
class SplineInterpolator2D
    : public IInterpolator2D<typename Spline2DBuilder::interpolation_domain_type, IdxRangeBatched>
{
    using base_type
            = IInterpolator2D<typename Spline2DBuilder::interpolation_domain_type, IdxRangeBatched>;

public:
    using typename base_type::CConstFieldType;
    using typename base_type::CoordType;
    using typename base_type::DFieldType;

private:
    Spline2DBuilder const& m_builder;

    Spline2DEvaluator const& m_evaluator;

    using IdxRangeBSRTheta =
            typename Spline2DBuilder::template batched_spline_domain_type<IdxRangeBatched>;

    mutable DFieldMem<IdxRangeBSRTheta> m_coefs;

    using r_deriv_type = DConstField<
            typename Spline2DBuilder::template batched_derivs_domain_type1<IdxRangeBatched>>;
    using theta_deriv_type = DConstField<
            typename Spline2DBuilder::template batched_derivs_domain_type2<IdxRangeBatched>>;
    using mixed_deriv_type = DConstField<
            typename Spline2DBuilder::template batched_derivs_domain_type<IdxRangeBatched>>;

public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known
     *                  interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     * @param[in] idx_range_batched The index range on which this operator operates.
     */
    SplineInterpolator2D(
            Spline2DBuilder const& builder,
            Spline2DEvaluator const& evaluator,
            IdxRangeBatched idx_range_batched)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(builder.batched_spline_domain(idx_range_batched))
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
 * This class allows an instance of the SplineInterpolator2D class to be instantiated where necessary. This allows the
 * memory allocated in the private members of the SplineInterpolator2D to be freed when the object is not in use.
 * These objects are: m_coefs.
 *
 * The class is parametrised by multiple template parameters. Please note that CTAD will deduce all these
 * template parameters from the Builder and Evaluator passed as constructor arguments.
 *
 * @tparam Spline2DBuilder The type of the 2D spline builder.
 * @tparam Spline2DEvaluator The type of the 2D spline evaluator.
 * @tparam IdxRangeBatched The type af the index range over which this operator will operate. This is
 *              necessary to define the internal spline representation.
 */
template <class Spline2DBuilder, class Spline2DEvaluator, class IdxRangeBatched>
class PreallocatableSplineInterpolator2D
    : public IPreallocatableInterpolator2D<
              typename Spline2DBuilder::interpolation_domain_type,
              IdxRangeBatched>
{
private:
    Spline2DBuilder const& m_builder;

    Spline2DEvaluator const& m_evaluator;

    IdxRangeBatched m_idx_range_batched;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolator2D objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at
     *                      known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     * @param[in] idx_range_batched The index range on which this operator operates.
     */
    PreallocatableSplineInterpolator2D(
            Spline2DBuilder const& builder,
            Spline2DEvaluator const& evaluator,
            IdxRangeBatched idx_range_batched)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_idx_range_batched(idx_range_batched)
    {
    }

    /**
     * Create a pointer to an instance of the SplineInterpolator2D class.
     *
     * @return A pointer to an instance of the SplineInterpolator2D class.
     */
    std::unique_ptr<
            IInterpolator2D<typename Spline2DBuilder::interpolation_domain_type, IdxRangeBatched>>
    preallocate() const override
    {
        return std::make_unique<SplineInterpolator2D<
                Spline2DBuilder,
                Spline2DEvaluator,
                IdxRangeBatched>>(m_builder, m_evaluator, m_idx_range_batched);
    }
};
