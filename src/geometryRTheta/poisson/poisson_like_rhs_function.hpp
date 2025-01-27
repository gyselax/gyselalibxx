// SPDX-License-Identifier: MIT
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
    host_t<ConstSpline2D> const m_coefs;
    evaluator_type const& m_evaluator;

public:
    /**
	 * @brief Instantiate a PoissonLikeRHSFunction.
	 *
	 * @param[in] coefs
	 *      The B-splines coefficients of the right-hand side function.
	 * @param[in] evaluator
	 *      Evaluator on B-splines.
	 */
    PoissonLikeRHSFunction(host_t<ConstSpline2D> coefs, evaluator_type const& evaluator)
        : m_coefs(coefs)
        , m_evaluator(evaluator)
    {
    }

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
