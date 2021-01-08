#pragma once

// library headers
#include <simulationobserver.hpp>

/** An implementation of SimulationObserver that prints the simulation data to
 * screen.
 * 
 * \warning This class is inefficient, brittle and does not guarantee the correct
 * ordering of the values printed on screen.
 */
class ScreenPrinter
		: public SimulationObserver
{
public:
	// see overriden function
	void simulation_updated( const Distributed2DField& data ) override;
	
};
