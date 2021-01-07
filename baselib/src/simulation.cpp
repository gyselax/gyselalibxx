// the implemented class (last)
#include "simulation.hpp"

Simulation::Simulation( MPI_Comm comm, const Configuration& config, const TimeStep& time_step, const InitialConditionner& init )
	: m_config( config )
	, m_time_step( time_step )
	, m_init( init)
	, m_comm( comm )
{
}

void Simulation::run() const
{
	// allocate data for the current iteration
	Distributed2DField current( m_comm, m_config.dist_extents(), m_config.global_shape(), m_time_step.required_ghosts(), m_config.delta_space() );

	// allocate data for the next iteration
	Distributed2DField next( m_comm, m_config.dist_extents(), m_config.global_shape(), m_time_step.required_ghosts(), m_config.delta_space() );

	// initialize data at t=0
	m_init.initial_condition( current );
	m_init.initial_condition( next );

	for ( auto&& observer : m_observers ) {
		observer->simulation_updated( current );
	}

	// the main (time) iteration
	for ( int ii = 0; ii < m_config.nb_iter(); ++ii ) {

		// compute the temperature at the next iteration
		m_time_step.iter( current, next );

		for ( auto&& observer : m_observers ) {
			observer->simulation_updated( next );
		}

		// swap the current and next buffers
		current.swap( next );
	}

}
