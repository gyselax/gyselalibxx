#pragma once

// local headers
#include "distributed2dfield.hpp"

/** An interface for a class that implements the operator applied at each and
 * every time-step of the Simulation.
 */
class TimeStep
{
public:
	/// The destructor
	virtual ~TimeStep() = default;

	/** Compute the temperature at t+delta_t based on the temperature at t
	* @param cur  the current value (t) of the local data block
	* @param next the next value (t+delta_t) of the local data block
	*/
	virtual void iter( const Distributed2DField& cur, Distributed2DField& next ) const = 0;
	
	/** Access the required ghost size to support this solver (size of the stencil)
	 * @return the required ghost size to support this solver (size of the stencil)
	 */
	virtual Shape2D required_ghosts() const = 0;
	
};
