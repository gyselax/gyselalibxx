// SPDX-License-Identifier: MIT

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_helper.hpp"
#include "neumann_spline_quadrature.hpp"
#include "quadrature.hpp"

namespace {

struct X
{
    static bool constexpr PERIODIC = false;
};

using CoordX = Coord<X>;

struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};

auto constexpr SplineXBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsX
        = ddc::KnotsAsInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};

using SplineXBuilder_1d_h = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        GridX>;
using SplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        GridX>;

using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

using DFieldMemX = FieldMem<double, IdxRangeX>;

TEST(NeumannSplineUniformQuadrature1D, ExactForConstantFuncHost)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IdxStepX const x_size(10);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());

    SplineXBuilder_1d_h const builder_x(gridx);

    host_t<DFieldMemX> const quadrature_coeffs(
            neumann_spline_quadrature_coefficients<
                    Kokkos::DefaultHostExecutionSpace>(gridx, builder_x));
    Quadrature const integrate(get_const_field(quadrature_coeffs));

    host_t<DFieldMemX> values(gridx);
    ddc::parallel_fill(values, 1.0);
    double integral = integrate(Kokkos::DefaultHostExecutionSpace(), get_const_field(values));
    double expected_val = x_max - x_min;
    EXPECT_LE(abs(integral - expected_val), 1e-15);
}


TEST(NeumannSplineUniformQuadrature1D, ExactForConstantFuncDefaultExecSpace)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IdxStepX const x_size(10);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());

    SplineXBuilder_1d const builder_x(gridx);

    DFieldMemX const quadrature_coeffs(neumann_spline_quadrature_coefficients<
                                       Kokkos::DefaultExecutionSpace>(gridx, builder_x));
    /* Quadrature const integrate(get_const_field(quadrature_coeffs));

    host_t<DFieldMemX> values(gridx);
    ddc::parallel_fill(values, 1.0);
    double integral = integrate(Kokkos::DefaultHostExecutionSpace(), get_const_field(values));
    double expected_val = x_max - x_min;
    EXPECT_LE(abs(integral - expected_val), 1e-15);*/
}

/*
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
    using GrevillePointsY = ddc::
            KnotsAsInterpolationPoints<BSplinesY, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;
    struct GridY : GrevillePointsY::interpolation_discrete_dimension_type
    {
    };
};

template <std::size_t N>
double compute_error(int n_elems)
{
    using DimY = typename ComputeErrorTraits<N>::Y;
    using BSplinesY = typename ComputeErrorTraits<N>::BSplinesY;
    using GrevillePointsY = typename ComputeErrorTraits<N>::GrevillePointsY;
    using GridY = typename ComputeErrorTraits<N>::GridY;
    using SplineYBuilder = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesY,
            GridY,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::SplineSolver::LAPACK,
            GridY>;
    using IdxRangeY = IdxRange<GridY>;
    using DFieldMemY = FieldMem<double, IdxRangeY>;
    using DFieldY = Field<double, IdxRangeY>;

    Coord<DimY> const y_min(-1.0);
    Coord<DimY> const y_max(1.0);

    ddc::init_discrete_space<BSplinesY>(y_min, y_max, n_elems);

    ddc::init_discrete_space<GridY>(GrevillePointsY::template get_sampling<GridY>());
    IdxRangeY const gridy(GrevillePointsY::template get_domain<GridY>());

    SplineYBuilder const builder_y(gridy);

    DFieldMemY const quadrature_coeffs(neumann_spline_quadrature_coefficients<
                                       Kokkos::DefaultExecutionSpace>(gridy, builder_y));
    Quadrature const integrate(get_field(quadrature_coeffs));

    DFieldMemY values_alloc(gridy);
    DFieldY values = get_field(values_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridy,
            KOKKOS_LAMBDA(Idx<GridY> const idx) {
                double x = ddc::coordinate(idx);
                values(idx) = (x + 1) * (x + 1) * (x + 1) * (x - 1) * (x - 1);
            });
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values);
    return std::abs(16.0 / 15.0 - integral);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2)...};
}

TEST(NeumannSplineUniformQuadrature1D, Convergence)
{
    constexpr int NTESTS(3);

    std::array<double, NTESTS> error = compute_errors(std::make_index_sequence<NTESTS>(), 10);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        double order_error = abs(4 - order);
        EXPECT_LE(order_error, 5e-2);
    }
}
*/
} // namespace
