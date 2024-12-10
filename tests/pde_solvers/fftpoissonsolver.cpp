// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <sll/math_tools.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "fft_poisson_solver.hpp"

namespace FFTPoissonSolverTest {

struct X
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

struct Y
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

using CoordX = Coord<X>;
using CoordY = Coord<Y>;

struct GridX : UniformGridBase<X>
{
};
struct GridY : UniformGridBase<Y>
{
};

using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxXY = Idx<GridX, GridY>;

using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;

using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;

using DFieldMemX = DFieldMem<IdxRangeX>;
using DFieldMemY = DFieldMem<IdxRangeY>;
using DFieldMemXY = DFieldMem<IdxRangeXY>;

using DFieldX = DField<IdxRangeX>;
using DFieldY = DField<IdxRangeY>;
using DFieldXY = DField<IdxRangeXY>;


TEST(FftPoissonSolver, CosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IdxStepX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_size + 1));
    IdxRangeX gridx = IdxRangeX(IdxX(0), x_size);

    // Creating operators
    FFTPoissonSolver<IdxRangeX, IdxRangeX, Kokkos::DefaultExecutionSpace> poisson(gridx);

    host_t<DFieldMemX> electrostatic_potential_host(gridx);
    host_t<DFieldMemX> electric_field_host(gridx);
    host_t<DFieldMemX> rhs_host(gridx);

    // Initialization of the distribution function --> fill values
    for (IdxX const ix : gridx) {
        rhs_host(ix) = cos(ddc::coordinate(ix));
    }
    DFieldMemX electrostatic_potential(gridx);
    DFieldMemX electric_field(gridx);
    DFieldMemX rhs(gridx);

    ddc::parallel_deepcopy(rhs, rhs_host);
    poisson(get_field(electrostatic_potential), get_field(electric_field), get_field(rhs));
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);

    double error_pot = 0.0;
    double error_field = 0.0;

    for (IdxX const ix : gridx) {
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
    IdxStepX const x_size(10);

    CoordY const y_min(0.0);
    CoordY const y_max(2.0 * M_PI);
    IdxStepY const y_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_size + 1));
    ddc::init_discrete_space<GridY>(GridY::init<GridY>(y_min, y_max, y_size + 1));
    IdxRangeX gridx = IdxRangeX(IdxX(0), x_size);
    IdxRangeY gridy = IdxRangeY(IdxY(0), y_size);
    IdxRangeXY gridxy(gridx, gridy);

    // Creating operators
    FFTPoissonSolver<IdxRangeY, IdxRangeXY, Kokkos::DefaultExecutionSpace> poisson(gridy);

    host_t<DFieldMemXY> electrostatic_potential_host(gridxy);
    host_t<DFieldMemXY> electric_field_host(gridxy);
    host_t<DFieldMemXY> rhs_host(gridxy);

    // Initialization of the distribution function --> fill values
    for (IdxX const ix : gridx) {
        double const c = (ix - gridx.front()) + 1;
        for (IdxY const iy : gridy) {
            rhs_host(ix, iy) = ipow(c, 2) * cos(c * ddc::coordinate(iy));
        }
    }
    DFieldMemXY electrostatic_potential(gridxy);
    DFieldMemXY electric_field(gridxy);
    DFieldMemXY rhs(gridxy);

    ddc::parallel_deepcopy(rhs, rhs_host);
    poisson(get_field(electrostatic_potential), get_field(electric_field), get_field(rhs));
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);

    double error_pot = 0.0;
    double error_field = 0.0;

    for (IdxX const ix : gridx) {
        double const c = (ix - gridx.front()) + 1;
        for (IdxY const iy : gridy) {
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
    IdxStepX const x_size(10);

    CoordY const y_min(0.0);
    CoordY const y_max(2.0 * M_PI);
    IdxStepY const y_size(10);

    // Creating mesh & supports
    ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_size + 1));
    ddc::init_discrete_space<GridY>(GridY::init<GridY>(y_min, y_max, y_size + 1));

    IdxRangeX gridx = IdxRangeX(IdxX(0), x_size);
    IdxRangeY gridy = IdxRangeY(IdxY(0), y_size);

    IdxRangeXY gridxy(gridx, gridy);

    FFTPoissonSolver<IdxRangeXY, IdxRangeXY, Kokkos::DefaultExecutionSpace> poisson(gridxy);

    DFieldMemXY electrostatic_potential_alloc(gridxy);
    VectorFieldMem<double, IdxRangeXY, NDTag<X, Y>> electric_field_alloc(gridxy);
    DFieldMemXY rhs_alloc(gridxy);

    DFieldXY electrostatic_potential = get_field(electrostatic_potential_alloc);
    VectorField electric_field = get_field(electric_field_alloc);
    DFieldXY rhs = get_field(rhs_alloc);

    // Initialization of the distribution function --> fill values
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                IdxX ix = ddc::select<GridX>(ixy);
                IdxY iy = ddc::select<GridY>(ixy);
                double x = ddc::coordinate(ix);
                double y = ddc::coordinate(iy);
                rhs(ixy) = Kokkos::cos(x) + Kokkos::cos(y);
            });

    poisson(get_field(electrostatic_potential_alloc),
            get_field(electric_field_alloc),
            get_field(rhs_alloc));

    double error_pot = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxXY const ixy) {
                IdxX ix = ddc::select<GridX>(ixy);
                IdxY iy = ddc::select<GridY>(ixy);
                double const exact_pot
                        = Kokkos::cos(ddc::coordinate(ix)) + Kokkos::cos(ddc::coordinate(iy));
                return Kokkos::abs(electrostatic_potential(ixy) - exact_pot);
            });
    DFieldXY electric_field_x = electric_field.get<X>();
    double error_field_x = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxXY const ixy) {
                IdxX ix = ddc::select<GridX>(ixy);
                double const exact_field_x = Kokkos::sin(ddc::coordinate(ix));
                return Kokkos::abs(electric_field_x(ixy) - exact_field_x);
            });
    DFieldXY electric_field_y = electric_field.get<Y>();
    double error_field_y = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxXY const ixy) {
                IdxY iy = ddc::select<GridY>(ixy);
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
