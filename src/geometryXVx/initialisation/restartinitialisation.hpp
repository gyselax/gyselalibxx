// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "iinitialisation.hpp"
#include "species_info.hpp"

/**
 * @brief A class that initialises the distribution function from a previous simulation.
 *
 * A class that triggers a PDI event to read the values of 
 * a distribution function saved in a hdf5 file. These
 * values are copied to the field that represents the 
 * distribution function. 
 */
class RestartInitialisation : public IInitialisation
{
private:
    int m_iter_start; /* iteration number to perform the restart from */
    double& m_time_start; /* corresponding simulation time */

public:
    /**
     * @brief Create an initialisation object.
     * @param[in] iter_start An integer representing the number of iteration already performed 
     *                       to produce the distribution function used to initialise the current simulation.
     * @param[in] time_start The physical time corresponding to iter_start.
     */
    RestartInitialisation(int iter_start, double& time_start);

    ~RestartInitialisation() override = default;

    /**
     * @brief Triggers a PDI event to fill the distribution function with values from a hdf5 file.
     * @param[out] allfdistribu The distribution function initialised with the values 
     *                          read from an external file.
     * @return The initialised distribution function.
     */
    DFieldSpXVx operator()(DFieldSpXVx allfdistribu) const override;
};
