// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

#include "ichargedensitycalculator.hpp"
#include "quadrature.hpp"
#include "simpson_quadrature.hpp"

/**
 * @brief A class which computes charges density.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \int_{v} q_s f_s(x,v) dv @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,v) @f$ is the distribution function.
 */
class ChargeDensityCalculator : public IChargeDensityCalculator
{
    const Quadrature<IDimVx>& m_quad;

public:
    /**
     * @brief Create a ChargeDensityCalculator object.
     *
     * @param quad The quadrature method which should be used.
     */
    ChargeDensityCalculator(const Quadrature<IDimVx>& quad);

    /**
     * @brief Computes the charge density rho from the distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     *
     * @return rho The charge density.
     */
    DSpanX operator()(DSpanX rho, DViewSpXVx allfdistribu) const final;
};
