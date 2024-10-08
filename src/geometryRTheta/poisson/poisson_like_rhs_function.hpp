#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "spline_interpolator_2d_rp.hpp"



/**
 * @brief Type of right-hand side (rhs) function of the Poisson equation.
 * @tparam RadialExtrapolationRule The extrapolation rule applied at the outer radial bound.
 */
template <class RadialExtrapolationRule>
class PoissonLikeRHSFunction
{
public:
    /// The type of the 2D Spline Evaluator used by this class
    using evaluator_type = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
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
    host_t<Spline2DConstField> const m_coefs;
    evaluator_type const& m_evaluator;

public:
    /**
	 * @brief Instantiate a PoissonLikeRHSFunction.
	 *
	 * @param[in] coefs
	 *      The bsplines coefficients of the right-hand side function.
	 * @param[in] evaluator
	 *      Evaluator on bsplines.
	 */
    PoissonLikeRHSFunction(host_t<Spline2DConstField> coefs, evaluator_type const& evaluator)
        : m_coefs(coefs)
        , m_evaluator(evaluator)
    {
    }

    ~PoissonLikeRHSFunction() {};

    /**
	 * @brief Get the value of the function at a given coordinate.
	 *
	 * @param[in] coord_rp
	 *      Polar coordinate where we want to evaluate the rhs function.
	 *
	 * @return A double with the value of the rhs at the given coordinate.
	 */
    double operator()(CoordRTheta const& coord_rp) const
    {
        return m_evaluator(coord_rp, m_coefs);
    }
};
