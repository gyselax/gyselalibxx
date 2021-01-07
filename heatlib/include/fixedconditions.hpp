#pragma once

// library headers
#include <initialconditionner.hpp>

class FixedConditions:
	public InitialConditionner
{
	double m_delta_t;

public:
	// see overridden function
	virtual void initial_condition( Distributed2DField& data ) const override;

};
