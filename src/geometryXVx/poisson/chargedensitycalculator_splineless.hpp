// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

#include "quadrature.hpp"
#include "simpson_quadrature.hpp"

/**
 * @brief A class which computes charges density .      
 */
class ChargeDensityCalculatorSplineless
{
    const Quadrature<IDimVx>& m_quad;

public:
    /**
 * @brief Create a ChargeDensityCalculatorSplineless object, integration is done using composite Simpson rule.
 */
    ChargeDensityCalculatorSplineless(const Quadrature<IDimVx>& quad);
    /**
 * @brief Computes the charge density rho from the distribution function.
 * @param[in, out] rho
 * @param[in] allfdistribu 
 */
    void operator()(DSpanX rho, DViewSpXVx allfdistribu) const;
};
