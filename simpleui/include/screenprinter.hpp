#pragma once

// library headers
#include <simulationobserver.hpp>

class ScreenPrinter
		: public SimulationObserver
{
public:
	// see overriden function
	void simulation_updated( const Distributed2DField& data ) override;
	
};
