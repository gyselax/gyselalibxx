// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "geometry.hpp"
#include "ichargedensitycalculator.hpp"
#include "quadrature.hpp"

/**
 * @brief A class which computes charges density with Kokkos.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \int_{vx} \int_{vy} q_s f_s(x,y,vx,vy) dvx dvy @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,y,vx,vy) @f$ is the distribution function.
 */
class MpiChargeDensityCalculator : public IChargeDensityCalculator
{
private:
    IChargeDensityCalculator const& m_local_charge_density_calculator;
    MPI_Comm m_comm;

public:
    /**
     * @brief Create a MpiChargeDensityCalculator object.
     * @param[in] coeffs
     *            The coefficients of the quadrature.
     */
    explicit MpiChargeDensityCalculator(
            MPI_Comm comm,
            IChargeDensityCalculator const& local_charge_density_calculator);

    /**
     * @brief Computes the charge density rho from the distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const final;
};
