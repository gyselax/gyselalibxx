#pragma once
#include "geometry.hpp"
#include "i_interpolator_2d_rp.hpp"


/**
 * @brief A class for interpolating a function using splines in polar coordinates.
 *
 * @tparam RadialExtrapolationRule The extrapolation rule applied at the outer radial bound.
 */
template <class RadialExtrapolationRule>
class SplineInterpolatorRP : public IInterpolatorRP
{
public:
    /// The type of the 2D Spline Evaluator used by this class
    using evaluator_type = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesP,
            IDimR,
            IDimP,
            RadialExtrapolationRule,
            RadialExtrapolationRule,
            ddc::PeriodicExtrapolationRule<RDimP>,
            ddc::PeriodicExtrapolationRule<RDimP>,
            IDimR,
            IDimP>;

private:
    SplineRPBuilder const& m_builder;

    evaluator_type const& m_evaluator;

    mutable ddc::Chunk<double, BSDomainRP> m_coefs;

    using r_deriv_type = ddc::ChunkSpan<double const, SplineRPBuilder::batched_derivs_domain_type1>;
    using p_deriv_type = ddc::ChunkSpan<double const, SplineRPBuilder::batched_derivs_domain_type2>;
    using mixed_deriv_type
            = ddc::ChunkSpan<double const, SplineRPBuilder::batched_derivs_domain_type>;

public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    SplineInterpolatorRP(SplineRPBuilder const& builder, evaluator_type const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(builder.spline_domain())
    {
    }

    ~SplineInterpolatorRP() override = default;

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
    DSpanRP operator()(
            DSpanRP const inout_data,
            ddc::ChunkSpan<CoordRP const, IDomainRP> const coordinates) const override
    {
#ifndef NDEBUG
        // To ensure that the interpolator is C0, we ensure that
        // the value at (r=0,theta) is the same for all theta.
        auto r_domain = ddc::get_domain<IDimR>(inout_data);
        auto theta_domain = ddc::get_domain<IDimP>(inout_data);
        if (ddc::coordinate(r_domain.front()) == 0) {
            ddc::for_each(theta_domain, [&](IndexP const ip) {
                assert(("Unicity of the value at the center point:",
                        inout_data(r_domain.front(), ip)
                                == inout_data(r_domain.front(), theta_domain.front())));
            });
        }
#endif

        m_builder(m_coefs.span_view(), inout_data.span_cview());
        m_evaluator(inout_data.span_view(), coordinates, m_coefs.span_cview());
        return inout_data;
    }
};



/**
 * @brief A class which stores information necessary to create a pointer to an instance of the SplineInterpolatorRP class.
 *
 * This class allows an instance of the SplineInterpolatorRP class where necessary. This allows the
 * memory allocated in the private members of the SplineInterpolatorRP to be freed when the object is not in use.
 * These objects are: m_coefs.
 */
template <class RadialExtrapolationRule>
class PreallocatableSplineInterpolatorRP : public IPreallocatableInterpolatorRP
{
public:
    /// The type of the 2D Spline Evaluator used by this class
    using evaluator_type = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesP,
            IDimR,
            IDimP,
            RadialExtrapolationRule,
            RadialExtrapolationRule,
            ddc::PeriodicExtrapolationRule<RDimP>,
            ddc::PeriodicExtrapolationRule<RDimP>,
            IDimR,
            IDimP>;

private:
    SplineRPBuilder const& m_builder;

    evaluator_type const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolatorRP objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSplineInterpolatorRP(
            SplineRPBuilder const& builder,
            evaluator_type const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    ~PreallocatableSplineInterpolatorRP() override = default;

    /**
     * Create a pointer to an instance of the SplineInterpolatorRP class.
     *
     * @return A pointer to an instance of the SplineInterpolatorRP class.
     */
    std::unique_ptr<IInterpolatorRP> preallocate() const override
    {
        return std::make_unique<
                SplineInterpolatorRP<RadialExtrapolationRule>>(m_builder, m_evaluator);
    }
};
