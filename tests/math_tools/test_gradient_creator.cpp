// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>

#include "../coord_transformations/coord_transformations_testing_tools.hpp"
#include "../coord_transformations/geometry_coord_transformations_tests.hpp"

#include "circular_to_cartesian.hpp"
#include "constant_partial_derivatives.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "gradient_creator.hpp"
#include "mesh_builder.hpp"

TEST(GradientCreator, Call)
{
    IdxStepR nbcells_r(4);
    std::vector<CoordR> point_sampling_r
            = build_uniform_break_points(CoordR(1e-5), CoordR(1.), nbcells_r);

    ddc::init_discrete_space<GridR>(point_sampling_r);

    IdxStepTheta nbcells_th(4);
    std::vector<CoordTheta> point_sampling_th
            = build_uniform_break_points(CoordTheta(0.), CoordTheta(2. * M_PI), nbcells_th);
    ddc::init_discrete_space<GridTheta>(point_sampling_th);

    IdxRangeR idxrange_r(IdxR(0), nbcells_r + 1);
    IdxRangeTheta idxrange_th(IdxTheta(0), nbcells_th);
    IdxRangeRTheta idxrange_rtheta(idxrange_r, idxrange_th);

    ConstantPartialDerivativeCreator<IdxRange<GridR, GridTheta>, R> r_deriv(0);
    ConstantPartialDerivativeCreator<IdxRange<GridR, GridTheta>, Theta> theta_deriv(1);
    GradientCreator<IdxRange<GridR, GridTheta>, R, Theta> get_gradient(r_deriv, theta_deriv);

    DVectorFieldMem<IdxRange<GridR, GridTheta>, VectorIndexSet<R_cov, Theta_cov>>
            grad_func_cov_alloc(idxrange_rtheta);
    DVectorField<IdxRange<GridR, GridTheta>, VectorIndexSet<R_cov, Theta_cov>> grad_func_cov(
            grad_func_cov_alloc);

    DFieldMem<IdxRange<GridR, GridTheta>> func_alloc(idxrange_rtheta);

    get_gradient(grad_func_cov, get_const_field(func_alloc));

    auto gradient_cov_host = ddcHelper::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), grad_func_cov);

    ddc::host_for_each(get_idx_range(gradient_cov_host), [&](IdxRTheta const irtheta) {
        EXPECT_EQ(ddcHelper::get<R_cov>(gradient_cov_host)(irtheta), 0.0);
        EXPECT_EQ(ddcHelper::get<Theta_cov>(gradient_cov_host)(irtheta), 1.0);
    });
}
