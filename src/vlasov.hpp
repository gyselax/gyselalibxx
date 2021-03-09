#pragma once

// local headers
#include "distributedfield.hpp"

/** An interface for a class that implements the operator applied at each and
 * every time-step of the Simulation.
 */
class Vlasov
{
public:
	/// The destructor
	virtual ~Vlasov() = default;

	virtual void operator()( const Field2D& cur, Field2D& next ) const = 0;
	
};

