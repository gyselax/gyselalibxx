#include <cstdlib>
#include <iostream>
#include <string>

#include "commandlineconfig.hpp"

using std::cerr;
using std::endl;
using std::exit;
using std::stoi;
using std::stod;

CommandLineConfig::CommandLineConfig( const int argc, const char* const argv[] )
{
	if ( argc != 9 ) {
		cerr << "Usage: " << argv[0] << " <Nb_iter> <height> <width> <process_height> <process_width> <delta_t> <delta_y> <delta_x>" << endl;
		exit( 1 );
	}

	m_nb_iter = stoi(argv[1]);
	m_global_shape[DY] = stoi(argv[2]);
	m_global_shape[DX] = stoi(argv[3]);
	m_dist_extents[DY] = stoi(argv[4]);
	m_dist_extents[DX] = stoi(argv[5]);
	m_delta_t = stod(argv[6]);
	m_delta_space[DY] = stod(argv[7]);
	m_delta_space[DX] = stod(argv[8]);
}
