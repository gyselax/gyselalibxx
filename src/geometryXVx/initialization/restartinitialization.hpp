// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <species_info.hpp>

#include "iinitialization.hpp"

/**
 * @brief A class that initializes the distribution function from a previous simulation.
 *
 * A class that triggers a PDI event to read the values of 
 * a distribution function saved in a hdf5 file. These
 * values are copied to the field that represents the 
 * distribution function. 
 */
class RestartInitialization : public IInitialization
{
private:
    int m_iter_start; /* iteration number to perform the restart from */
    double& m_time_start; /* corresponding simulation time */

public:
    /**
     * @brief Create an initialization object.
     * @param[in] iter_start An integer representing the number of iteration already performed 
     *                       to produce the distribution function used to initialize the current simulation.
     * @param[in] time_start The physical time corresponding to iter_start.
     */
    RestartInitialization(int iter_start, double& time_start);

    ~RestartInitialization() override = default;

    /**
     * @brief Triggers a PDI event to fill the distribution function with values from a hdf5 file.
     * @param[out] allfdistribu The distribution function initialized with the values 
     *                          read from an external file.
     * @return The initialized distribution function.
     */
    DFieldSpXVx operator()(DFieldSpXVx allfdistribu) const override;
};
