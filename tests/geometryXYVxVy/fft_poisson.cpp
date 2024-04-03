// SPDX-License-Identifier: MIT

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "chargedensitycalculator.hpp"
#include "fftpoissonsolver.hpp"
#include "geometry.hpp"
#include "neumann_spline_quadrature.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"

namespace {

static void TestFftPoissonSolverCosineSource()
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

    IDomainXYVxVy const meshxyvxvy(gridxy, gridvxvy);

    SplineVxBuilder_1d const builder_vx_1d(gridvx);
    SplineVyBuilder_1d const builder_vy_1d(gridvy);

    IDomainSpXYVxVy const mesh(gridsp, meshxyvxvy);

    ddc::init_discrete_space<IDimFx>(ddc::init_fourier_space<RDimX>(gridx));
    ddc::init_discrete_space<IDimFy>(ddc::init_fourier_space<RDimY>(gridy));

    // Initialise infomation about species
    host_t<FieldSp<int>> charges(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    host_t<DFieldSp> masses(dom_sp);
    ddc::parallel_fill(masses, 1);

    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));

    host_t<DFieldVxVy> const quadrature_coeffs_host
            = neumann_spline_quadrature_coefficients(gridvxvy, builder_vx_1d, builder_vy_1d);
    auto quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    ChargeDensityCalculator rhs(quadrature_coeffs);
    FftPoissonSolver poisson(rhs);

    DFieldXY electrostatic_potential_alloc(gridxy);
    DFieldXY electric_field_x_alloc(gridxy);
    DFieldXY electric_field_y_alloc(gridxy);
    DFieldSpXYVxVy allfdistribu_alloc(mesh);

    DSpanXY electrostatic_potential = electrostatic_potential_alloc.span_view();
    DSpanXY electric_field_x = electric_field_x_alloc.span_view();
    DSpanXY electric_field_y = electric_field_y_alloc.span_view();
    DSpanSpXYVxVy allfdistribu = allfdistribu_alloc.span_view();

    // Initialization of the distribution function --> fill values
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            mesh,
            KOKKOS_LAMBDA(IndexSpXYVxVy const ispxyvxvy) {
                IndexX ix = ddc::select<IDimX>(ispxyvxvy);
                IndexY iy = ddc::select<IDimY>(ispxyvxvy);
                double x = ddc::coordinate(ix);
                double y = ddc::coordinate(iy);
                double fdistribu_val = Kokkos::cos(x) + Kokkos::cos(y);
                allfdistribu(ispxyvxvy) = fdistribu_val;
            });

    poisson(electrostatic_potential, electric_field_x, electric_field_y, allfdistribu);

    double error_pot = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IndexXY const ixy) {
                IndexX ix = ddc::select<IDimX>(ixy);
                IndexY iy = ddc::select<IDimY>(ixy);
                double const exact_pot
                        = Kokkos::cos(ddc::coordinate(ix)) + Kokkos::cos(ddc::coordinate(iy));
                return Kokkos::abs(electrostatic_potential(ixy) - exact_pot);
            });
    double error_field_x = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IndexXY const ixy) {
                IndexX ix = ddc::select<IDimX>(ixy);
                double const exact_field_x = Kokkos::sin(ddc::coordinate(ix));
                return Kokkos::abs(electric_field_x(ixy) - exact_field_x);
            });
    double error_field_y = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IndexXY const ixy) {
                IndexY iy = ddc::select<IDimY>(ixy);
                double const exact_field_y = Kokkos::sin(ddc::coordinate(iy));
                return Kokkos::abs(electric_field_y(ixy) - exact_field_y);
            });

    EXPECT_LE(error_pot, 1e-12);
    EXPECT_LE(error_field_x, 1e-10);
    EXPECT_LE(error_field_y, 1e-10);
}

} // namespace

TEST(FftPoissonSolver, CosineSource)
{
    TestFftPoissonSolverCosineSource();
}
