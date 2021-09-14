#include <cstdlib>

#include <ddc/MCoord>
#include <ddc/ProductMDomain>
#include <ddc/ProductMesh>

#include <paraconf.h>

#include "geometry.h"
#include "nulladvectionvx.h"
#include "nullefieldsolver.h"
#include "predcorr.h"
#include "splineadvectionx.h"
#include "splitvlasovsolver.h"

int main(int argc, char** argv)
{
    PC_tree_t conf;
    if (argc > 1) {
        conf = PC_parse_path(argv[1]);
    } else {
        return EXIT_SUCCESS;
    }

    long steps;
    PC_int(PC_get(conf, ".steps"), &steps);

    double time_step;
    PC_double(PC_get(conf, ".time_step"), &time_step);

    double mass_ratio;
    PC_double(PC_get(conf, ".mass_ratio"), &mass_ratio);

    double x_min;
    PC_double(PC_get(conf, ".MeshX.min"), &x_min);
    double x_max;
    PC_double(PC_get(conf, ".MeshX.max"), &x_max);
    long x_size;
    PC_int(PC_get(conf, ".MeshX.size"), &x_size);

    double vx_min;
    PC_double(PC_get(conf, ".MeshVx.min"), &vx_min);
    double vx_max;
    PC_double(PC_get(conf, ".MeshVx.max"), &vx_max);
    long vx_size;
    PC_int(PC_get(conf, ".MeshVx.size"), &vx_size);

    KnotsX const knots_x((RCoordX(x_min)), RCoordX(x_max), x_size);

    MDomainX dom_knots_x((ProductMesh<KnotsX>(knots_x)), MCoord<KnotsX>(x_size - 1));

    BSplinesX bsplines_x(dom_knots_x);

    SplineXBuilder builder_x(bsplines_x);

    MeshVx const mesh_vx((RCoordVx(vx_min)), RCoordVx(vx_max), vx_size);

    MeshXVx const mesh_x_vx(builder_x.interpolation_domain().mesh().get<MeshX>(), mesh_vx);

    MDomainXVx const
            dom2d(mesh_x_vx,
                  MCoordXVx(builder_x.interpolation_domain().extents() - 1, vx_size - 1));

    SplineAdvectionX const advection_x(bsplines_x, builder_x);

    NullAdvectionVx const advection_vx;

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    NullEfieldSolver const efield;

    PredCorr const predcorr(vlasov, efield, time_step);

    DBlockXVx fdistribu(dom2d);

    predcorr(fdistribu, mass_ratio, steps);

    return EXIT_SUCCESS;
}
