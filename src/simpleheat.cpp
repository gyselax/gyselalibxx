#include <mpi.h>

#include "distributed2dfield.hpp"
#include "advection1d.hpp"
#include "bsplines_uniform.h"
#include "spline_interpolator_1d.h"

int main( int argc, char* argv[] )
{
	// initialize the MPI library
	MPI_Init( &argc, &argv );

	DataND<2> f2d({1,1}, MPI_COMM_WORLD, Dimension::DX, {1.,1.});

	DataND<1> ex({1}, MPI_COMM_WORLD, Dimension::DX, {1.});

    int degree(3);
    bool periodic(false);
    double xmin(0.0), xmax(1.0);
    int ncells(10);
    BSplines_uniform bspl(degree, periodic, xmin, xmax, ncells);
    Spline_interpolator_1D spl_interp(bspl, BoundaryCondition::sll_p_greville, BoundaryCondition::sll_p_greville);

    double dt(0.1);

    Advection1D adv(bspl, spl_interp, dt,
                [](double x) { return 0.0; }, [](double x) { return 0.0; });

	// finalize MPI
	MPI_Finalize();

	return 0;
}
