#pragma once

// standard C++ library headers
#include <set>

// local headers
#include "configuration.hpp"
#include "distributed2dfield.hpp"
#include "initialconditionner.hpp"
#include "simulationobserver.hpp"
#include "timestep.hpp"

/** The main simulation class that implements time-stepping.
 * 
 * The simulation is executed with the run() function.
 * The class relies on external services to provide most of its behaviour.
 * A Configuration instance provide the user-specified parameters, a 
 * InitialConditionner instance provides the initial and boundary conditions
 * while a TimeStep instance implements the operator to apply at each time-step.
 * 
 * The simulation is observable by SimulationObserver instances thanks to the
 * observe() and unobserve() functions.
 */
class Simulation
{
private:
	/// The configuration of the simulation
	const Configuration& m_config;

	/// The time-step operator to apply
	const TimeStep& m_time_step;
	
	/// The initial-conditions to apply
	const InitialConditionner& m_init;

	/// The MPI communicator to use
	MPI_Comm m_comm;

	/// The observers interested to be notified when the simulation state changes
	std::set<SimulationObserver*> m_observers;

public:
	/** Constructs a new simulation
	 * @param comm the MPI communicator to use
	 * @param config the configuration of the simulation
	 * @param time_step the time-step operator to apply
	 * @param init the initial-conditions to apply
	 */
	Simulation( MPI_Comm comm, const Configuration& config, const TimeStep& time_step, const InitialConditionner& init );

	/** Run the simulation for the number of time-steps specified in the config
	 */
	void run() const;

	/** Add an observer interested by the simulation
	 * @param observer the observer to add
	 */
	void observe( SimulationObserver& observer ) { m_observers.insert( &observer ); }

	/** Remove an observer no longer interested by the simulation
	 * @param observer the observer to remove
	 */
	void unobserve( SimulationObserver& observer ) { m_observers.erase( &observer ); }

};
