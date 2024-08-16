// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <quadrature.hpp>

#include "ichargedensitycalculator.hpp"

/**
 * @brief A class which computes charges density with Kokkos.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \sum_s \int_{v} q_s f_s(x,v) dv @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,v) @f$ is the distribution function.
 */
class ChargeDensityCalculator : public IChargeDensityCalculator
{
private:
    Quadrature<IdxRangeVx, IdxRangeXVx> m_quadrature;

public:
    /**
     * @brief Create a ChargeDensityCalculator object.
     * @param[in] coeffs
     *            The coefficients of the quadrature.
     */
    explicit ChargeDensityCalculator(DConstFieldVx coeffs);

    /**
     * @brief Computes the charge density rho from the distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     *
     * @return rho The charge density.
     */
    DFieldX operator()(DFieldX rho, DConstFieldSpXVx allfdistribu) const final;
};
