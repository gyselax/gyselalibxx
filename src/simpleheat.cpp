#include <mpi.h>

#include "advection1d.hpp"
#include "bsplines_uniform.h"
#include "distributedfield.hpp"
#include "spline_interpolator_1d.h"

int main(int argc, char* argv[])
{
    // initialize the MPI library
    MPI_Init(&argc, &argv);

    Field2D f2d;

    Field2D ex;

    int degree(3);
    bool periodic(false);
    double xmin(0.0), xmax(1.0);
    int ncells(10);
    BSplines_uniform bspl(degree, periodic, xmin, xmax, ncells);
    Spline_interpolator_1D spl_interp(
        bspl, BoundaryCondition::sll_p_greville, BoundaryCondition::sll_p_greville);

    NullBoundaryValue& bv = NullBoundaryValue::value;

    double dt(0.1);

    Advection1D adv(bspl, spl_interp, bv, bv);

    // finalize MPI
    MPI_Finalize();

    return 0;
}
