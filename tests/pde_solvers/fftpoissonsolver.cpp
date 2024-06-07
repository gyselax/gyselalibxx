// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <sll/math_tools.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "fft_poisson_solver.hpp"

namespace FFTPoissonSolverTest {

struct RDimX
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

struct RDimY
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;

struct IDimX : ddc::UniformPointSampling<RDimX>
{
};
struct IDimY : ddc::UniformPointSampling<RDimY>
{
};

using IndexX = ddc::DiscreteElement<IDimX>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;

using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;

using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;

using DFieldX = device_t<ddc::Chunk<double, IDomainX>>;
using DFieldY = device_t<ddc::Chunk<double, IDomainY>>;
using DFieldXY = device_t<ddc::Chunk<double, IDomainXY>>;

using DSpanX = device_t<ddc::ChunkSpan<double, IDomainX>>;
using DSpanY = device_t<ddc::ChunkSpan<double, IDomainY>>;
using DSpanXY = device_t<ddc::ChunkSpan<double, IDomainXY>>;


TEST(FftPoissonSolver, CosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IVectX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_size + 1));
    IDomainX gridx = IDomainX(IndexX(0), x_size);

    // Creating operators
    FFTPoissonSolver<IDomainX, IDomainX, Kokkos::DefaultExecutionSpace> poisson(gridx);

    host_t<DFieldX> electrostatic_potential_host(gridx);
    host_t<DFieldX> electric_field_host(gridx);
    host_t<DFieldX> rhs_host(gridx);

    // Initialization of the distribution function --> fill values
    for (IndexX const ix : gridx) {
        rhs_host(ix) = cos(ddc::coordinate(ix));
    }
    DFieldX electrostatic_potential(gridx);
    DFieldX electric_field(gridx);
    DFieldX rhs(gridx);

    ddc::parallel_deepcopy(rhs, rhs_host);
    poisson(electrostatic_potential, electric_field, rhs);
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);

    double error_pot = 0.0;
    double error_field = 0.0;

    for (IndexX const ix : gridx) {
        double const exact_pot = cos(ddc::coordinate(ix));
        error_pot = fmax(fabs(electrostatic_potential_host(ix) - exact_pot), error_pot);
        double const exact_field = sin(ddc::coordinate(ix));
        error_field = fmax(fabs(electric_field_host(ix) - exact_field), error_field);
    }
    EXPECT_LE(error_pot, 1e-8);
    EXPECT_LE(error_field, 1e-6);
}

TEST(FftPoissonSolver, BatchedCosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IVectX const x_size(10);

    CoordY const y_min(0.0);
    CoordY const y_max(2.0 * M_PI);
    IVectY const y_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_size + 1));
    ddc::init_discrete_space<IDimY>(IDimY::init<IDimY>(y_min, y_max, y_size + 1));
    IDomainX gridx = IDomainX(IndexX(0), x_size);
    IDomainY gridy = IDomainY(IndexY(0), y_size);
    IDomainXY gridxy(gridx, gridy);

    // Creating operators
    FFTPoissonSolver<IDomainY, IDomainXY, Kokkos::DefaultExecutionSpace> poisson(gridy);

    host_t<DFieldXY> electrostatic_potential_host(gridxy);
    host_t<DFieldXY> electric_field_host(gridxy);
    host_t<DFieldXY> rhs_host(gridxy);

    // Initialization of the distribution function --> fill values
    for (IndexX const ix : gridx) {
        double const c = ix.uid() + 1;
        for (IndexY const iy : gridy) {
            rhs_host(ix, iy) = ipow(c, 2) * cos(c * ddc::coordinate(iy));
        }
    }
    DFieldXY electrostatic_potential(gridxy);
    DFieldXY electric_field(gridxy);
    DFieldXY rhs(gridxy);

    ddc::parallel_deepcopy(rhs, rhs_host);
    poisson(electrostatic_potential, electric_field, rhs);
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);

    double error_pot = 0.0;
    double error_field = 0.0;

    for (IndexX const ix : gridx) {
        double const c = ix.uid() + 1;
        for (IndexY const iy : gridy) {
            double const exact_pot = cos(c * ddc::coordinate(iy));
            error_pot = fmax(fabs(electrostatic_potential_host(ix, iy) - exact_pot), error_pot);
            double const exact_field = c * sin(c * ddc::coordinate(iy));
            error_field = fmax(fabs(electric_field_host(ix, iy) - exact_field), error_field);
        }
    }
    EXPECT_LE(error_pot, 1e-8);
    EXPECT_LE(error_field, 1e-6);
}

static void TestFftPoissonSolver2DCosineSource()
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IVectX const x_size(10);

    CoordY const y_min(0.0);
    CoordY const y_max(2.0 * M_PI);
    IVectY const y_size(10);

    // Creating mesh & supports
    ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_size + 1));
    ddc::init_discrete_space<IDimY>(IDimY::init<IDimY>(y_min, y_max, y_size + 1));

    IDomainX gridx = IDomainX(IndexX(0), x_size);
    IDomainY gridy = IDomainY(IndexY(0), y_size);

    IDomainXY gridxy(gridx, gridy);

    FFTPoissonSolver<IDomainXY, IDomainXY, Kokkos::DefaultExecutionSpace> poisson(gridxy);

    DFieldXY electrostatic_potential_alloc(gridxy);
    VectorField<double, IDomainXY, NDTag<RDimX, RDimY>, ddc::DeviceAllocator<double>>
            electric_field_alloc(gridxy);
    DFieldXY rhs_alloc(gridxy);

    DSpanXY electrostatic_potential = electrostatic_potential_alloc.span_view();
    VectorFieldSpan electric_field = electric_field_alloc.span_view();
    DSpanXY rhs = rhs_alloc.span_view();

    // Initialization of the distribution function --> fill values
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(IndexXY const ixy) {
                IndexX ix = ddc::select<IDimX>(ixy);
                IndexY iy = ddc::select<IDimY>(ixy);
                double x = ddc::coordinate(ix);
                double y = ddc::coordinate(iy);
                rhs(ixy) = Kokkos::cos(x) + Kokkos::cos(y);
            });

    poisson(electrostatic_potential_alloc, electric_field_alloc, rhs_alloc);

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
    ddc::ChunkSpan electric_field_x = electric_field.get<RDimX>();
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
    ddc::ChunkSpan electric_field_y = electric_field.get<RDimY>();
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

TEST(FftPoissonSolver2D, CosineSource)
{
    TestFftPoissonSolver2DCosineSource();
}

} // namespace FFTPoissonSolverTest
