#pragma once

// library headers
#include <initialconditionner.hpp>

/** An implementation of InitialConditionner that simply sets predefined initial
 * and Dirichlet boundary conditions:
 * - at t=0, 0 is set on the whole domain,
 * - at the boundary on top, right, and bottom, the value 0 is set also
 * - at the boundary on the left, 2097.152 is used.
 */
class FixedConditions:
	public InitialConditionner
{
	// see overridden function
	virtual void initial_condition( Distributed2DField& data ) const override;

};
