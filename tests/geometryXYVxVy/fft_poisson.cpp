// SPDX-License-Identifier: MIT

#include <sll/null_boundary_value.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "chargedensitycalculator.hpp"
#include "fftpoissonsolver.hpp"
#include "geometry.hpp"
#include "neumann_spline_quadrature.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"

TEST(FftPoissonSolver, CosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IVectX const x_size(10);

    CoordY const y_min(0.0);
    CoordY const y_max(2.0 * M_PI);
    IVectY const y_size(10);

    CoordVx const vx_min(-0.5);
    CoordVx const vx_max(0.5);
    IVectVx const vx_size(11);

    CoordVy const vy_min(-0.5);
    CoordVy const vy_max(0.5);
    IVectVy const vy_size(11);

    IVectSp const nb_species(2);
    IDomainSp const dom_sp(IndexSp(0), nb_species);
    IndexSp const my_ielec = dom_sp.front();
    IndexSp const my_iion = dom_sp.back();

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);
    ddc::init_discrete_space<BSplinesVy>(vy_min, vy_max, vy_size);
    ddc::init_discrete_space<DDCBSplinesVx>(vx_min, vx_max, vx_size);
    ddc::init_discrete_space<DDCBSplinesVy>(vy_min, vy_max, vy_size);

    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());
    ddc::init_discrete_space<IDimVy>(SplineInterpPointsVy::get_sampling());

    ddc::init_discrete_space(IDimX::init(x_min, x_max, x_size + 1));
    ddc::init_discrete_space(IDimY::init(y_min, y_max, y_size + 1));

    IDomainSp const gridsp = IDomainSp(my_iion, IVectSp(1));

    IDomainX gridx = IDomainX(IndexX(0), x_size);
    IDomainY gridy = IDomainY(IndexY(0), y_size);

    IDomainXY gridxy(gridx, gridy);

    IDomainVx gridvx = SplineInterpPointsVx::get_domain();
    IDomainVy gridvy = SplineInterpPointsVy::get_domain();

    IDomainVxVy gridvxvy(gridvx, gridvy);

    ddc::init_fourier_space<RDimX>(gridx);
    ddc::init_fourier_space<RDimY>(gridy);

    SplineVxVyBuilder const builder_vx_vy(gridvxvy);
    SplineVxVyEvaluator const evaluator_vx_vy(
            g_null_boundary_2d<BSplinesVx, BSplinesVy>,
            g_null_boundary_2d<BSplinesVx, BSplinesVy>,
            g_null_boundary_2d<BSplinesVx, BSplinesVy>,
            g_null_boundary_2d<BSplinesVx, BSplinesVy>);

    SplineVxBuilder_1d const builder_vx_1d(gridvx);
    SplineVyBuilder_1d const builder_vy_1d(gridvy);

    IDomainSpXYVxVy const mesh(gridsp, gridxy, gridvxvy);

    // Initialise infomation about species
    FieldSp<int> charges(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    DFieldSp masses(dom_sp);
    ddc::parallel_fill(masses, 1);

    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));

    host_t<DFieldVxVy> const quadrature_coeffs_host
            = neumann_spline_quadrature_coefficients(gridvxvy, builder_vx_1d, builder_vy_1d);
    auto quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    ChargeDensityCalculator rhs(quadrature_coeffs);
    FftPoissonSolver poisson(rhs);

    DFieldXY electrostatic_potential(gridxy);
    DFieldXY electric_field_x(gridxy);
    DFieldXY electric_field_y(gridxy);
    DFieldSpXYVxVy allfdistribu(mesh);

    // Initialization of the distribution function --> fill values
    auto c_dom = ddc::remove_dims_of(mesh, gridxy);
    ddc::for_each(gridxy, [&](IndexXY const ixy) {
        IndexX ix = ddc::select<IDimX>(ixy);
        IndexY iy = ddc::select<IDimY>(ixy);
        double x = ddc::coordinate(ix);
        double y = ddc::coordinate(iy);
        double fdistribu_val = cos(x) + cos(y);
        ddc::for_each(c_dom, [&](auto const ispvxvy) {
            allfdistribu(ixy, ispvxvy) = fdistribu_val;
        });
    });

    poisson(electrostatic_potential, electric_field_x, electric_field_y, allfdistribu);

    double error_pot = 0.0;
    double error_field_x = 0.0;
    double error_field_y = 0.0;

    ddc::for_each(gridxy, [&](IndexXY const ixy) {
        IndexX ix = ddc::select<IDimX>(ixy);
        IndexY iy = ddc::select<IDimY>(ixy);
        double const exact_pot = cos(ddc::coordinate(ix)) + cos(ddc::coordinate(iy));
        error_pot = fmax(fabs(electrostatic_potential(ixy) - exact_pot), error_pot);
        double const exact_field_x = sin(ddc::coordinate(ix));
        double const exact_field_y = sin(ddc::coordinate(iy));
        error_field_x = fmax(fabs(electric_field_x(ixy) - exact_field_x), error_field_x);
        error_field_y = fmax(fabs(electric_field_y(ixy) - exact_field_y), error_field_y);
    });
    EXPECT_LE(error_pot, 1e-12);
    EXPECT_LE(error_field_x, 1e-10);
    EXPECT_LE(error_field_y, 1e-10);
}
