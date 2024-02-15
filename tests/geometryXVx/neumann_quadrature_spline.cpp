// SPDX-License-Identifier: MIT

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "geometry.hpp"
#include "neumann_spline_quadrature.hpp"
#include "quadrature.hpp"

namespace {

TEST(NeumannSplineQuadratureTest, ExactForConstantFunc)
{
    CoordVx const vx_min(0.0);
    CoordVx const vx_max(M_PI);
    IVectVx const vx_size(10);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());
    ddc::DiscreteDomain<IDimVx> interpolation_domain_vx(SplineInterpPointsVx::get_domain());

    SplineVxBuilder_1d const builder_vx(interpolation_domain_vx);

    IDomainVx const gridvx = builder_vx.interpolation_domain();

    DFieldVx const quadrature_coeffs = neumann_spline_quadrature_coefficients(gridvx, builder_vx);
    Quadrature const integrate(quadrature_coeffs.span_cview());

    DFieldVx values(gridvx);

    ddc::for_each(gridvx, [&](ddc::DiscreteElement<IDimVx> const idx) { values(idx) = 1.0; });
    double integral = integrate(values);
    double expected_val = vx_max - vx_min;
    EXPECT_LE(abs(integral - expected_val), 1e-15);
}


template <std::size_t N>
struct Y
{
    static bool constexpr PERIODIC = false;
};

template <std::size_t N>
double compute_error(int n_elems)
{
    using DimY = Y<N>;
    using BSplinesY = ddc::UniformBSplines<DimY, 3>;
    using GrevillePointsY = ddc::
            KnotsAsInterpolationPoints<BSplinesY, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;
    using IDimY = typename GrevillePointsY::interpolation_mesh_type;
    using SplineYBuilder = ddc::SplineBuilder<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesY,
            IDimY,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::SplineSolver::GINKGO,
            IDimY>;
    using IDomainY = ddc::DiscreteDomain<IDimY>;
    using DFieldY = ddc::Chunk<double, IDomainY>;

    ddc::Coordinate<Y<N>> const y_min(-1.0);
    ddc::Coordinate<Y<N>> const y_max(1.0);

    ddc::init_discrete_space<BSplinesY>(y_min, y_max, n_elems);

    ddc::init_discrete_space<IDimY>(GrevillePointsY::get_sampling());
    IDomainY const gridy(GrevillePointsY::get_domain());

    SplineYBuilder const builder_y(gridy);

    DFieldY const quadrature_coeffs = neumann_spline_quadrature_coefficients(gridy, builder_y);
    Quadrature const integrate(quadrature_coeffs.span_cview());

    DFieldY values(gridy);

    ddc::for_each(gridy, [&](ddc::DiscreteElement<IDimY> const idx) {
        double x = ddc::coordinate(idx);
        values(idx) = (x + 1) * (x + 1) * (x + 1) * (x - 1) * (x - 1);
    });
    double integral = integrate(values);
    return std::abs(16.0 / 15.0 - integral);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2)...};
}

TEST(NeumannSplineQuadratureTest, UniformConverge)
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

} // namespace
