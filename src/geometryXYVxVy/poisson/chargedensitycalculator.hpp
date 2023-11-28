// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "ichargedensitycalculator.hpp"

/**
 * @brief A class which computes charges density.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \int_{v} q_s f_s(x,v) dv @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,v) @f$ is the distribution function.
 *
 * This equation is solved using an intermediate spline representation.
 */
class ChargeDensityCalculator : public IChargeDensityCalculator
{
    SplineVxVyBuilder const& m_spline_vxvy_builder;

    SplineVxVyEvaluator m_spline_vxvy_evaluator;

    int m_nbc_Vx;
    int m_nbc_Vy;
    int m_interp_dom_size_Vx;
    int m_interp_dom_size_Vy;

public:
    /**
     * Constructor of SplineChargeDensityCalculator
     *
     * @param spline_vxvy_builder A spline builder which calculates the coefficients of a spline representation.
     * @param spline_vxvy_evaluator A spline evaluator which provides the value of a spline representation from its coefficients.
     */
    ChargeDensityCalculator(
            SplineVxVyBuilder const& spline_vxvy_builder,
            SplineVxVyEvaluator const& spline_vxvy_evaluator);

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
     */
    void operator()(DSpanXY rho, DViewSpXYVxVy allfdistribu) const final;
};
