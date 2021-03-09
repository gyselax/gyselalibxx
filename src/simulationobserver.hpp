#pragma once

// local headers
#include "distributedfield.hpp"

/** An interface for any class interested in observing a Simulation.
 */
class SimulationObserver {
public:
    /// The destructor
    ~SimulationObserver() = default;

    /** Notification when the state of the simulation changes
   * @param data the new state of the simulation
   */
    virtual void simulation_updated(const Data2D& data) = 0;
};
