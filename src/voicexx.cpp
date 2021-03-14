#include <vector>

#include "bsplines_uniform.h"
#include "mesh.h"
#include "predcorr.h"
#include "space.h"
#include "vlasov.h"

int main()
{
    // a 2D geometry,
    RDomain2D geometry = {{0, 100}, {0, 100}};

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

    const Vlasov vlasov {advection_x, advection_v};

    const EfieldSolver* pefield = nullptr;

    const PredCorr predcorr = {vlasov, *pefield};

    Mesh2D mesh {Mesher {0, 1}, Mesher {0, 1}};

    MDomain2D dom2d(mesh, {0, 0}, ExtentsND<2> {10, 10});

    DBlock2D fdistribu = {dom2d};

    Mesh1D time_mesh {Mesher {0, .1}};

    predcorr(fdistribu, MDomain1D {time_mesh, {0}, Extents1D {100}});
}
