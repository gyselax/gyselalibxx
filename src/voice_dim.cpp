#include <vector>

#include "advectionx.h"
#include "bsplines_uniform.h"
#include "nullefieldsolver.h"
#include "predcorr2.h"
#include "spline_interpolator_1d.h"
#include "splitvlasovsolver.h"

int main()
{
    // The meshed domain:
    // * origin: (0,0)
    // * unit vector: (1,1)
    // * domain-start: (0,0)
    // * domain-bound: (100, 100)
    MDomainXVx dom2d(RCoordXVx(0., 0.), RCoordXVx(1., 1.), MCoordXVx(0, 0), MCoordXVx(100, 100));

    const BSplines_uniform bsplines_x = {3, true, 4, 5, 6};

    const BSplines_uniform bsplines_vx = {3, true, 4, 5, 6};

    const Spline_interpolator_1D interp_x(bsplines_x, BoundCond::PERIODIC, BoundCond::PERIODIC);

    const Spline_interpolator_1D interp_vx(bsplines_vx, BoundCond::GREVILLE, BoundCond::GREVILLE);

    //TODO: const AdvectionX advection_x = {bsplines_x, interp_x};

    //TODO: const AdvectionX advection_v = {bsplines_vx, interp_vx};

    //TODO: const SplitVlasovSolver vlasov {advection_x, advection_v};

    const NullEfieldSolver efield;

    //TODO: const PredCorr2 predcorr(vlasov, efield, 1.);

    DBlockXVx fdistribu(dom2d);

    //TODO: predcorr(fdistribu, 1, 100);
}
