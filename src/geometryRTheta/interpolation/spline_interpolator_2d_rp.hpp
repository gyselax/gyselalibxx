// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "i_interpolator_2d_rp.hpp"


/**
 * @brief A class for interpolating a function using splines in polar coordinates.
 *
 * @tparam RadialExtrapolationRule The extrapolation rule applied at the outer radial bound.
 */
template <class RadialExtrapolationRule>
class SplineInterpolatorRTheta : public IInterpolatorRTheta
{
public:
    /// The type of the 2D Spline Evaluator used by this class
    using evaluator_type = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
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
    SplineRThetaBuilder_host const& m_builder;

    evaluator_type const& m_evaluator;

    mutable host_t<DFieldMem<IdxRangeBSRTheta>> m_coefs;

    using r_deriv_type = host_t<DConstField<SplineRThetaBuilder_host::batched_derivs_domain_type1>>;
    using p_deriv_type = host_t<DConstField<SplineRThetaBuilder_host::batched_derivs_domain_type2>>;
    using mixed_deriv_type
            = host_t<DConstField<SplineRThetaBuilder_host::batched_derivs_domain_type>>;

public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    SplineInterpolatorRTheta(
            SplineRThetaBuilder_host const& builder,
            evaluator_type const& evaluator)
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
    host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> const inout_data,
            host_t<Field<CoordRTheta const, IdxRangeRTheta>> const coordinates) const override
    {
#ifndef NDEBUG
        // To ensure that the interpolator is C0, we ensure that
        // the value at (r=0,theta) is the same for all theta.
        IdxRangeR r_idx_range = get_idx_range<GridR>(inout_data);
        IdxRangeTheta theta_idx_range = get_idx_range<GridTheta>(inout_data);
        if (ddc::coordinate(r_idx_range.front()) == 0) {
            ddc::for_each(theta_idx_range, [&](IdxTheta const ip) {
                bool const unicity_center_point
                        = inout_data(r_idx_range.front(), ip)
                          == inout_data(r_idx_range.front(), theta_idx_range.front());
                if (!unicity_center_point) {
                    std::printf("Unicity of the value at the center point is not verified.");
                    assert(unicity_center_point);
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
class PreallocatableSplineInterpolatorRTheta : public IPreallocatableInterpolatorRTheta
{
public:
    /// The type of the 2D Spline Evaluator used by this class
    using evaluator_type = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
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
    SplineRThetaBuilder_host const& m_builder;

    evaluator_type const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolatorRTheta objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSplineInterpolatorRTheta(
            SplineRThetaBuilder_host const& builder,
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
    std::unique_ptr<IInterpolatorRTheta> preallocate() const override
    {
        return std::make_unique<
                SplineInterpolatorRTheta<RadialExtrapolationRule>>(m_builder, m_evaluator);
    }
};
