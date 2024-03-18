// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>
#include <geometry.hpp>

#include "ichargedensitycalculator.hpp"

/**
 * @brief A class which computes charges density with Kokkos.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \int_{vx} \int_{vy} q_s f_s(x,y,vx,vy) dvx dvy @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,y,vx,vy) @f$ is the distribution function.
 */
class ChargeDensityCalculator : public IChargeDensityCalculator
{
private:
    using ChunkViewType = device_t<DViewVxVy>;
    ChunkViewType m_coefficients;

public:
    /**
     * @brief Create a ChargeDensityCalculator object.
     * @param[in] coeffs
     *            The coefficients of the quadrature.
     */
    explicit ChargeDensityCalculator(const ChunkViewType& coeffs);

    /**
     * @brief Computes the charge density rho from the distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DSpanXY rho, DViewSpXYVxVy allfdistribu) const final;
};
