#pragma once

// local headers
#include "distributed2dfield.hpp"

/** An interface for any class interested in observing a Simulation.
 */
class SimulationObserver
{
public:
	/// The destructor
	~SimulationObserver() = default;
	
	/** Notification when the state of the simulation changes
	 * @param data the new state of the simulation
	 */
	virtual void simulation_updated( const Distributed2DField& data ) = 0;
	
};
