// SPDX-License-Identifier: MIT

#pragma once

#include <mpi.h>

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "geometry.hpp"
#include "imomentscalculator.hpp"
#include "quadrature.hpp"

/**
 * @brief A class which computes charges density with Kokkos.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \int_{vx} \int_{vy} q_s f_s(x,y,vx,vy) dvx dvy @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,y,vx,vy) @f$ is the distribution function.
 */
class MpiMomentsCalculator : public IMomentsCalculator
{
private:
    IMomentsCalculator const& m_local_moments_calculator;
    MPI_Comm m_comm;

public:
    /**
     * @brief Create a MpiMomentsCalculator object.
     * @param[in] comm The MPI communicator across which the calculation is carried out.
     * @param[in] local_moments_calculator
     *                 An operator which calculates the density locally
     *                 on a given MPI node. The results from this operator
     *                 will then be combined using MPI.
     */
    explicit MpiMomentsCalculator(
            MPI_Comm comm,
            IMomentsCalculator const& local_moments_calculator);

    /**
     * @brief Computes the charge density rho from the distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldX rho, DConstFieldSpVxVyVzX allfdistribu) const final;

    /**
     * @brief Computes the kinetic energy density rho from the distribution function.
     * @param[in, out] kinetic
     * @param[in] allfdistribu 
     */
    void operator()(DFieldX kinetic, DConstFieldSpVxVyVzX allfdistribu, char ki) const final;

    /**
     * @brief Computes the mean velocities from the distribution function.
     * @param[out] mean_current_x, mean_current_y
     * @param[in] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldX mean_current_x, DFieldX mean_current_y, DFieldX mean_current_z, DFieldX rho, 
                    DConstFieldSpVxVyVzX allfdistribu) const final;

    void operator()(DFieldX momentum_x, DFieldX momentum_y, DFieldX momentum_z, 
                    DConstFieldSpVxVyVzX allfdistribu) const final;

    /**
     * @brief Computes the charge density rho for each distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldSpX rho, DConstFieldSpVxVyVzX allfdistribu) const final;

    /**
     * @brief Computes the mean velocities for each distribution function.
     * @param[out] mean_current_x, mean_current_y
     * @param[in] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldSpX mean_current_x, DFieldSpX mean_current_y, DFieldSpX mean_current_z, DFieldSpX rho, 
                    DConstFieldSpVxVyVzX allfdistribu) const final;
};
