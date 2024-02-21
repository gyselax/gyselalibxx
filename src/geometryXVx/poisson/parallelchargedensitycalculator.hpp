// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

#include "ichargedensitycalculator.hpp"

/**
 * @brief A class which computes charges density with Kokkos.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \int_{v} q_s f_s(x,v) dv @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,v) @f$ is the distribution function.
 */
class ParallelChargeDensityCalculator : public IChargeDensityCalculator
{
private:
    using ChunkViewType = device_t<DViewVx>;
    ChunkViewType m_coefficients;

public:
    /**
     * @brief Create a ParallelChargeDensityCalculator object.
     * @param[in] coeffs
     *            The coefficients of the quadrature.
     */
    explicit ParallelChargeDensityCalculator(const ChunkViewType& coeffs);

    /**
     * @brief Computes the charge density rho from the distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     *
     * @return rho The charge density.
     */
    device_t<DSpanX> operator()(device_t<DSpanX> rho, device_t<DViewSpXVx> allfdistribu)
            const final;
};
