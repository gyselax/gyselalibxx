#include <vector>

#include "bsplines_uniform.h"
#include "efieldsolver.h"
#include "predcorr.h"
#include "vlasovsolver.h"

int main()
{
    // The mesh dx:0.01, dv:0.2, dt:.1
    Mesh3D mesh = {Mesh(0, .01), Mesh(-10, .2), Mesh(0, .1)};

    // The meshed domain x:[0,1) v:[-10,10) t:[0,100)
    MDomain3D dom3d = {MDomain(mesh[0], 100), MDomain(mesh[1], 100), MDomain(mesh[3], 1000)};

    const BSplines_uniform bsplines_x = {3, true, 4, 5, 6};

    const BSplines_uniform bsplines_vx = {3, true, 4, 5, 6};

    const Spline_interpolator_1D interp_x(
            bsplines_x,
            BoundaryCondition::sll_p_periodic,
            BoundaryCondition::sll_p_periodic);

    const Spline_interpolator_1D interp_vx(
            bsplines_vx,
            BoundaryCondition::sll_p_greville,
            BoundaryCondition::sll_p_greville);

    const Advection1D advection_x = {bsplines_x, interp_x};

    const Advection1D advection_v = {bsplines_vx, interp_vx};

    const VlasovSolver vlasov {advection_x, advection_v};

    const EfieldSolver efield;

    const PredCorr predcorr(vlasov, efield, {dom3d[3]});

    DBlock2D fdistribu({dom3d[0], dom3d[1]});

    predcorr(fdistribu, 1);
}
