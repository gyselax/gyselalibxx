#pragma once

// library headers
#include <configuration.hpp>
#include <timestep.hpp>

class FinitediffHeatSolver:
	public TimeStep
{
	// The time-step to use
	double m_delta_t;

public:
	/** Constructs a new FinitediffHeatSolver
	 * @param config the configuration where to look for the time-step to use
	 */
	FinitediffHeatSolver( const Configuration& config );

	// see overridden function
	void iter( const Distributed2DField& cur, Distributed2DField& next ) const override;

	// see overridden function
	Shape2D required_ghosts() const override { return {1, 1}; }

};
