// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <collisions_intra.hpp>
#include <collisions_utils.hpp>
#include <geometry.hpp>
#include <irighthandside.hpp>
#include <pdi.h>
#include <quadrature.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

/**
 * Tests the construction of the gridvx and gridvx_staggered 
 */
TEST(CollisionsIntraGridvx, CollisionsIntraGridvx)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(10);

    CoordVx const vx_min(-7);
    CoordVx const vx_max(7);
    IdxStepVx const vx_size(10);

    IdxStepSp const nb_kinspecies(2);

    IdxRangeSp const dom_sp(IdxSp(0), nb_kinspecies);
    IdxSp const my_iion = dom_sp.front();
    IdxSp const my_ielec = dom_sp.back();

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx gridvx(SplineInterpPointsVx::get_domain<GridVx>());

    SplineXBuilder_1d const builder_x(gridx);
    SplineVxBuilder_1d const builder_vx(gridvx);

    IdxRangeSpXVx const mesh(dom_sp, gridx, gridvx);

    host_t<DFieldMemSp> charges(dom_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    host_t<DFieldMemSp> masses(dom_sp);
    double const mass_ion(400.);
    double const mass_elec(1.);
    masses(my_ielec) = mass_elec;
    masses(my_iion) = mass_ion;
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    // collision operator
    double const nustar0(1.);
    CollisionsIntra collisions(mesh, nustar0);
    IdxRange<CollisionsIntra::GhostedVx> gridvx_ghosted = collisions.get_gridvx_ghosted();
    IdxRange<CollisionsIntra::GhostedVxStaggered> gridvx_ghosted_staggered
            = collisions.get_gridvx_ghosted_staggered();

    double const npoints(gridvx.size());
    std::vector<double> gridvx_ghosted_pred(npoints + 2);
    std::vector<double> gridvx_ghosted_staggered_pred(npoints + 1);

    double const dv((vx_max - vx_min) / vx_size);
    gridvx_ghosted_pred[0] = vx_min - dv;
    gridvx_ghosted_pred[npoints] = vx_max;
    gridvx_ghosted_pred[npoints + 1] = vx_max + dv;

    gridvx_ghosted_staggered_pred[0] = vx_min - dv / 2.;
    gridvx_ghosted_staggered_pred[npoints] = vx_max + dv / 2.;
    for (int i(1); i < npoints; i++) {
        gridvx_ghosted_pred[i] = vx_min + dv * (i - 1);
        gridvx_ghosted_staggered_pred[i] = vx_min - dv / 2. + dv * i;
    }

    ddc::for_each(gridvx_ghosted, [&](auto const ivx_gh) {
        EXPECT_LE(std::fabs(ddc::coordinate(ivx_gh) - gridvx_ghosted_pred[ivx_gh.uid()]), 1.e-12);
    });

    ddc::for_each(gridvx_ghosted_staggered, [&](auto const ivx_ghs) {
        EXPECT_LE(
                std::fabs(ddc::coordinate(ivx_ghs) - gridvx_ghosted_staggered_pred[ivx_ghs.uid()]),
                1.e-12);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
