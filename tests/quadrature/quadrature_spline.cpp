// SPDX-License-Identifier: MIT
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc/ddc.hpp"
#include "ddc/kernels/splines.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "quadrature.hpp"
#include "spline_quadrature.hpp"

namespace {

struct X
{
#ifdef PERIODIC_RDIMX
    static bool constexpr PERIODIC = true;
#else
    static bool constexpr PERIODIC = false;
#endif
};

using CoordX = Coord<X>;

struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};

auto constexpr SplineXBoundary = X::PERIODIC ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};

using IdxStepX = IdxStep<X>;
using IdxRangeX = IdxRange<GridX>;
using DFieldMemX = DFieldMem<IdxRangeX>;


TEST(SplineUniformQuadrature, ExactForConstantFunc)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IdxStepX const x_size(10);

    using SplineXBuilder = ddc::SplineBuilder<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesX,
            GridX,
            SplineXBoundary,
            SplineXBoundary,
            ddc::SplineSolver::LAPACK,
            GridX>;

    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());

    SplineXBuilder const builder_x(gridx);

    DFieldMemX const quadrature_coeffs(
            spline_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx, builder_x));
    Quadrature const integrate(get_const_field(quadrature_coeffs));

    DFieldMemX values(gridx);

    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), values, 1.0);
    double integral = integrate(Kokkos::DefaultExecutionSpace(), get_const_field(values));
    double expected_val = x_max - x_min;
    EXPECT_LE(abs(integral - expected_val), 1e-15);
}

template <std::size_t N>
struct ComputeErrorTraits
{
    struct Y
    {
        static bool constexpr PERIODIC = false;
    };
    struct BSplinesY : ddc::UniformBSplines<Y, 3>
    {
    };
    using GrevillePointsY = ddc::GrevilleInterpolationPoints<
            BSplinesY,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>;
    struct GridY : GrevillePointsY::interpolation_discrete_dimension_type
    {
    };
};

template <std::size_t N>
double compute_error(int n_elems)
{
    using Y = typename ComputeErrorTraits<N>::Y;
    using BSplinesY = typename ComputeErrorTraits<N>::BSplinesY;
    using GrevillePointsY = typename ComputeErrorTraits<N>::GrevillePointsY;
    using GridY = typename ComputeErrorTraits<N>::GridY;
    auto constexpr SplineYBoundary = ddc::BoundCond::GREVILLE;
    using SplineYBuilder = ddc::SplineBuilder<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesY,
            GridY,
            SplineYBoundary,
            SplineYBoundary,
            ddc::SplineSolver::LAPACK,
            GridY>;
    using IdxRangeY = IdxRange<GridY>;
    using DFieldMemY = DFieldMem<IdxRangeY>;
    using DFieldY = DField<IdxRangeY>;

    Coord<Y> const y_min(0.0);
    Coord<Y> const y_max(M_PI);

    ddc::init_discrete_space<BSplinesY>(y_min, y_max, n_elems);

    ddc::init_discrete_space<GridY>(GrevillePointsY::template get_sampling<GridY>());
    IdxRangeY const gridy(GrevillePointsY::template get_domain<GridY>());

    SplineYBuilder const builder_y(gridy);

    DFieldMemY quadrature_coeffs(
            spline_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridy, builder_y));
    Quadrature<IdxRangeY> const integrate(get_const_field(quadrature_coeffs));

    DFieldMemY values_alloc(gridy);
    DFieldY values = get_field(values_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridy,
            KOKKOS_LAMBDA(Idx<GridY> const idx) {
                values(idx) = Kokkos::sin(ddc::coordinate(idx));
            });
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values);
    return std::abs(2 - integral);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2)...};
}

TEST(SplineUniformQuadrature, Convergence)
{
    constexpr int NTESTS(6);

    std::array<double, NTESTS> error = compute_errors(std::make_index_sequence<NTESTS>(), 10);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        double order_error = abs(4 - order);
        EXPECT_LE(order_error, 5e-2);
    }
}

} // namespace
