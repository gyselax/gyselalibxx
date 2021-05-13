#include <vector>

#include "bsplines_uniform.h"
#include "nulladvectionvx.h"
#include "nullefieldsolver.h"
#include "predcorr.h"
#include "spline_interpolator_1d.h"
#include "splineadvectionx.h"
#include "splitvlasovsolver.h"

// The meshed domain:
// * origin: (0,0)
// * unit vector: (1,1)
// * domain-start: (0,0)
// * domain-bound: (100, 100)
MDomainXVx const dom2d(RCoordXVx(0., 0.), RCoordXVx(1., 1.), MCoordXVx(0, 0), MCoordXVx(100, 100));

BSplines_uniform const bsplines_x = {3, true, 4, 5, 6};

BSplines_uniform const bsplines_vx = {3, true, 4, 5, 6};

Spline_interpolator_1D const interp_x(bsplines_x, BoundCond::PERIODIC, BoundCond::PERIODIC);

Spline_interpolator_1D const interp_vx(bsplines_vx, BoundCond::GREVILLE, BoundCond::GREVILLE);

SplineAdvectionX const advection_x(bsplines_x, interp_x);

NullAdvectionVx const advection_vx;

SplitVlasovSolver const vlasov(advection_x, advection_vx);

NullEfieldSolver const efield;

PredCorr const predcorr(vlasov, efield, 1.);

int main()
{
    DBlockXVx fdistribu(dom2d);

    predcorr(fdistribu, 1, 100);
}
