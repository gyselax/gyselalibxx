#include <mpi.h>

#include "distributed2dfield.hpp"

int main( int argc, char* argv[] )
{
	// initialize the MPI library
	MPI_Init( &argc, &argv );

	DataND<2> f2d({1,1}, MPI_COMM_WORLD, Dimension::DX, {1.,1.});

	DataND<1> ex({1}, MPI_COMM_WORLD, Dimension::DX, {1.});

	// finalize MPI
	MPI_Finalize();

	return 0;
}
