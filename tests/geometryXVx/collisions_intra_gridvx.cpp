// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <collisions_intra.hpp>
#include <collisions_utils.hpp>
#include <geometry.hpp>
#include <irighthandside.hpp>
#include <maxwellianequilibrium.hpp>
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
    IVectX const x_size(10);

    CoordVx const vx_min(-7);
    CoordVx const vx_max(7);
    IVectVx const vx_size(10);

    IVectSp const nb_kinspecies(2);

    IDomainSp const dom_sp(IndexSp(0), nb_kinspecies);
    IndexSp const my_iion = dom_sp.front();
    IndexSp const my_ielec = dom_sp.back();

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(InterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(InterpPointsVx::get_sampling());

    IDomainX interpolation_domain_x(InterpPointsX::get_domain());
    IDomainVx interpolation_domain_vx(InterpPointsVx::get_domain());

    SplineXBuilder const builder_x(interpolation_domain_x);

    SplineVxBuilder const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSpXVx const mesh(dom_sp, gridx, gridvx);

    FieldSp<int> charges(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    DFieldSp masses(dom_sp);
    double const mass_ion(400);
    double const mass_elec(1);
    masses(my_ielec) = mass_elec;
    masses(my_iion) = mass_ion;
    FieldSp<int> init_perturb_mode(dom_sp);
    ddc::fill(init_perturb_mode, 0);
    DFieldSp init_perturb_amplitude(dom_sp);
    ddc::fill(init_perturb_amplitude, 0);
    ddc::init_discrete_space<IDimSp>(
            std::move(charges),
            std::move(masses),
            std::move(init_perturb_amplitude),
            std::move(init_perturb_mode));

    // collision operator
    double const nustar0(1.);
    CollisionsIntra collisions(mesh, nustar0);
    ddc::DiscreteDomain<CollisionsIntra::ghosted_vx_point_sampling> gridvx_ghosted
            = collisions.get_gridvx_ghosted();
    ddc::DiscreteDomain<CollisionsIntra::ghosted_vx_staggered_point_sampling>
            gridvx_ghosted_staggered = collisions.get_gridvx_ghosted_staggered();

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

    ddc::for_each(ddc::policies::parallel_host, gridvx_ghosted, [&](auto const ivx_gh) {
        EXPECT_LE(std::fabs(ddc::coordinate(ivx_gh) - gridvx_ghosted_pred[ivx_gh.uid()]), 1.e-12);
    });

    ddc::for_each(ddc::policies::parallel_host, gridvx_ghosted_staggered, [&](auto const ivx_ghs) {
        EXPECT_LE(
                std::fabs(ddc::coordinate(ivx_ghs) - gridvx_ghosted_staggered_pred[ivx_ghs.uid()]),
                1.e-12);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
