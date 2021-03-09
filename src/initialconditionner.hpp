#pragma once

// local headers
#include "distributedfield.hpp"

/** An interface to set the initial and boundary conditions of the simulation
 */
class InitialConditionner
{
public:
	/// The destructor
	virtual ~InitialConditionner() = default;
	
	/** Initializes the temperature at t=0
	* @param data the local data block to initialize
	*/
	virtual void initial_condition( Data2D& data ) const = 0;
	
};

