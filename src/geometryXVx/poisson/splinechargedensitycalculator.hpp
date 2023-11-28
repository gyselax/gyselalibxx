// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "geometry.hpp"
#include "ichargedensitycalculator.hpp"

/**
 * A class which calculates the charge density.
 *
 * A class which calculates the charge density using an intermediate spline representation.
 * The charge density is then used as the right hand side of the Poisson equation.
 *
 * When using a spline representation where no boundary conditions need to be provided
 * (e.g. with periodic or Greville boundary conditions) then this class is equivalent
 * to the ChargeDensityCalculator when initialised with a quadrature based on splines.
 * The version using the quadrature should always be preferred over this version as in
 * this implementation the matrix equation is solved at each call of the operator. This
 * makes it significantly more costly than the quadrature method.
 *
 * When using a spline representation where boundary conditions need to be provided (e.g
 * Hermite) this class provides zeros (derivative = 0.)
 */
class SplineChargeDensityCalculator : public IChargeDensityCalculator
{
    SplineVxBuilder const& m_spline_vx_builder;

    SplineEvaluator<BSplinesVx> const& m_spline_vx_evaluator;

    std::vector<double> m_derivs_vxmin_data;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_data;

    Span1D<double> m_derivs_vxmax;

public:
    /**
     * Constructor of SplineChargeDensityCalculator
     *
     * @param spline_vx_builder A spline builder which calculates the coefficients of a spline representation.
     * @param spline_vx_evaluator A spline evaluator which provides the value of a spline representation from its coefficients.
     */
    SplineChargeDensityCalculator(
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    /**
     * Calculate the charge density rho from the distribution function.
     *
     * Calculate the charge density by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is then integrated and multiplied by the charge to find the
     * charge density.
     *
     * @param[out] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     *
     * @return rho The charge density.
     */
    DSpanX operator()(DSpanX rho, DViewSpXVx allfdistribu) const final;
};
