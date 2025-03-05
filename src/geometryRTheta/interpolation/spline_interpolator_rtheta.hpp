// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "i_interpolator_rtheta.hpp"

/*
 *ddc::SplineEvaluator2D<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            RadialExtrapolationRule,
            RadialExtrapolationRule,
            ddc::PeriodicExtrapolationRule<Theta>,
            ddc::PeriodicExtrapolationRule<Theta>,
            GridR,
            GridTheta>;
*/

/**
 * @brief A class for interpolating a function using splines in polar coordinates.
 *
 * @tparam RadialExtrapolationRule The extrapolation rule applied at the outer radial bound.
 */
template <class SplineRThetaBuilder, class SplineRThetaEvaluator>
class SplineInterpolatorRTheta
    : public IInterpolator2D<
              typename SplineRThetaBuilder::interpolation_domain_type,
              typename SplineRThetaBuilder::batched_interpolation_domain_type>
{
private:
    SplineRThetaBuilder const& m_builder;

    SplineRThetaEvaluator const& m_evaluator;

    using IdxRangeBSRTheta = typename SplineRThetaBuilder::batched_spline_domain_type;

    mutable DFieldMem<IdxRangeBSRTheta> m_coefs;

    using r_deriv_type = DConstField<typename SplineRThetaBuilder::batched_derivs_domain_type1>;
    using theta_deriv_type = DConstField<typename SplineRThetaBuilder::batched_derivs_domain_type2>;
    using mixed_deriv_type = DConstField<typename SplineRThetaBuilder::batched_derivs_domain_type>;

public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    SplineInterpolatorRTheta(
            SplineRThetaBuilder const& builder,
            SplineRThetaEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(get_spline_idx_range(builder))
    {
    }

    ~SplineInterpolatorRTheta() override = default;

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
    DFieldRTheta operator()(
            DFieldRTheta const inout_data,
            ConstField<CoordRTheta, IdxRangeRTheta> const coordinates) const override
    {
#ifndef NDEBUG
        // To ensure that the interpolator is C0, we ensure that
        // the value at (r=0,theta) is the same for all theta.
        IdxRangeR r_idx_range = get_idx_range<GridR>(inout_data);
        IdxRangeTheta theta_idx_range = get_idx_range<GridTheta>(inout_data);
        if (ddc::coordinate(r_idx_range.front()) == 0) {
            ddc::parallel_for_each(
                    theta_idx_range,
                    KOKKOS_LAMBDA(IdxTheta const itheta) {
                        bool const unicity_centre_point
                                = inout_data(r_idx_range.front(), itheta)
                                  == inout_data(r_idx_range.front(), theta_idx_range.front());
                        if (!unicity_centre_point) {
                            Kokkos::abort(
                                    "Unicity of the value at the centre point is not verified.");
                        }
                    });
        }
#endif

        m_builder(get_field(m_coefs), get_const_field(inout_data));
        m_evaluator(get_field(inout_data), coordinates, get_const_field(m_coefs));

        return inout_data;
    }
};



/**
 * @brief A class which stores information necessary to create a pointer to an instance of the SplineInterpolatorRTheta class.
 *
 * This class allows an instance of the SplineInterpolatorRTheta class where necessary. This allows the
 * memory allocated in the private members of the SplineInterpolatorRTheta to be freed when the object is not in use.
 * These objects are: m_coefs.
 */
template <class RadialExtrapolationRule>
class PreallocatableSplineInterpolatorRTheta
    : public IPreallocatableInterpolator2D<IdxRangeRTheta, IdxRangeRTheta>
{
public:
    /// The type of the 2D Spline Evaluator used by this class
    using evaluator_type = ddc::SplineEvaluator2D<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            RadialExtrapolationRule,
            RadialExtrapolationRule,
            ddc::PeriodicExtrapolationRule<Theta>,
            ddc::PeriodicExtrapolationRule<Theta>,
            GridR,
            GridTheta>;

private:
    SplineRThetaBuilder const& m_builder;

    evaluator_type const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolatorRTheta objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSplineInterpolatorRTheta(
            SplineRThetaBuilder const& builder,
            evaluator_type const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    ~PreallocatableSplineInterpolatorRTheta() override = default;

    /**
     * Create a pointer to an instance of the SplineInterpolatorRTheta class.
     *
     * @return A pointer to an instance of the SplineInterpolatorRTheta class.
     */
    std::unique_ptr<IInterpolator2D<IdxRangeRTheta, IdxRangeRTheta>> preallocate() const override
    {
        return std::make_unique<SplineInterpolatorRTheta<
                SplineRThetaBuilder,
                evaluator_type>>(m_builder, m_evaluator);
    }
};
