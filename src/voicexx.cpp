#include <cstdlib>
#include <iosfwd>

#include <paraconf.h>

#include "deprecated/bsplines_uniform.h"

#include "boundary_conditions.h"
#include "geometry.h"
#include "nulladvectionvx.h"
#include "nullefieldsolver.h"
#include "predcorr.h"
#include "spline_builder_1d.h"
#include "splineadvectionx.h"
#include "splitvlasovsolver.h"

static constexpr std::size_t spline_degree = 3;

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

    MeshX const mesh_x((RCoordX(x_min)), RCoordX(x_max), x_size + 1);

    MeshVx const mesh_vx((RCoordVx(vx_min)), RCoordVx(vx_max), vx_size + 1);

    MeshXVx const mesh_x_vx(mesh_x, mesh_vx);

    MDomainXVx const dom2d(mesh_x_vx, MCoordXVx(x_size, vx_size));

    deprecated::UniformBSplines const bsplines_x(spline_degree, true, x_min, x_max, x_size);

    deprecated::UniformBSplines const bsplines_vx(spline_degree, true, vx_min, vx_max, vx_size);

    deprecated::SplineBuilder1D const
            interp_x(bsplines_x, BoundCond::PERIODIC, BoundCond::PERIODIC);

    deprecated::SplineBuilder1D const
            interp_vx(bsplines_vx, BoundCond::PERIODIC, BoundCond::PERIODIC);

    SplineAdvectionX const advection_x(bsplines_x, interp_x);

    NullAdvectionVx const advection_vx;

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    NullEfieldSolver const efield;

    PredCorr const predcorr(vlasov, efield, time_step);

    DBlockXVx fdistribu(dom2d);

    predcorr(fdistribu, mass_ratio, steps);

    return EXIT_SUCCESS;
}
