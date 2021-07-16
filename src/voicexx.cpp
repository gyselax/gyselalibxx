#include "deprecated/bsplines_uniform.h"

#include "boundary_conditions.h"
#include "geometry.h"
#include "nulladvectionvx.h"
#include "nullefieldsolver.h"
#include "predcorr.h"
#include "spline_builder_1d.h"
#include "splineadvectionx.h"
#include "splitvlasovsolver.h"

// The meshed domain:
// * origin: (0,0)
// * unit vector: (1,1)
// * domain-start: (0,0)
// * domain-bound: (100, 100)
MeshX const mesh_x(RCoordX(0.), RCoordX(1.));

MeshVx const mesh_vx(RCoordVx(0.), RCoordVx(1.));

MeshXVx const mesh_x_vx(mesh_x, mesh_vx);

MDomainXVx const dom2d(mesh_x_vx, MCoordXVx(0, 0), MCoordXVx(100, 100));

deprecated::UniformBSplines const bsplines_x = {3, true, 4, 5, 6};

deprecated::UniformBSplines const bsplines_vx = {3, true, 4, 5, 6};

deprecated::SplineBuilder1D const interp_x(bsplines_x, BoundCond::PERIODIC, BoundCond::PERIODIC);

deprecated::SplineBuilder1D const interp_vx(bsplines_vx, BoundCond::GREVILLE, BoundCond::GREVILLE);

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
