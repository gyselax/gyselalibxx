#include <vector>

#include "bsplines_uniform.h"
#include "efieldsolver.h"
#include "mesh.h"
#include "predcorr.h"
#include "space.h"
#include "vlasovsolver.h"

int main()
{
    // The mesh, d_x=.01, d_v=.2, d_t=.1
    Mesh3D mesh = {Mesher(0, .01), Mesher(-10, .2), Mesher(0, .1)};

    // The meshed domain
    MDomain3D dom3d
            = {MDomain(mesh[0], 0, 100), MDomain(mesh[1], 0, 100), MDomain(mesh[3], 0, 1000)};


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
