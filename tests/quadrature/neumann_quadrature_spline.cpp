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

using CoordX = ddc::Coordinate<X>;

struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};

auto constexpr SplineXBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsX
        = ddc::KnotsAsInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;

struct IDimX : SplineInterpPointsX::interpolation_mesh_type
{
};

using SplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX>;

using IVectX = ddc::DiscreteVector<IDimX>;
using IDomainX = ddc::DiscreteDomain<IDimX>;

using DFieldX = ddc::Chunk<double, IDomainX>;

TEST(NeumannSplineUniformQuadrature1D, ExactForConstantFunc)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IVectX const x_size(10);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::DiscreteDomain<IDimX> interpolation_domain_x(SplineInterpPointsX::get_domain<IDimX>());

    SplineXBuilder_1d const builder_x(interpolation_domain_x);

    IDomainX const gridx = builder_x.interpolation_domain();

    DFieldX const quadrature_coeffs = neumann_spline_quadrature_coefficients(gridx, builder_x);
    Quadrature const integrate(quadrature_coeffs.span_cview());

    DFieldX values(gridx);

    ddc::for_each(gridx, [&](ddc::DiscreteElement<IDimX> const idx) { values(idx) = 1.0; });
    double integral = integrate(values);
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
    using GrevillePointsY = ddc::
            KnotsAsInterpolationPoints<BSplinesY, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;
    struct IDimY : GrevillePointsY::interpolation_mesh_type
    {
    };
};

template <std::size_t N>
double compute_error(int n_elems)
{
    using DimY = typename ComputeErrorTraits<N>::Y;
    using BSplinesY = typename ComputeErrorTraits<N>::BSplinesY;
    using GrevillePointsY = typename ComputeErrorTraits<N>::GrevillePointsY;
    using IDimY = typename ComputeErrorTraits<N>::IDimY;
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
    using DFieldY = device_t<ddc::Chunk<double, IDomainY>>;

    ddc::Coordinate<DimY> const y_min(-1.0);
    ddc::Coordinate<DimY> const y_max(1.0);

    ddc::init_discrete_space<BSplinesY>(y_min, y_max, n_elems);

    ddc::init_discrete_space<IDimY>(GrevillePointsY::template get_sampling<IDimY>());
    IDomainY const gridy(GrevillePointsY::template get_domain<IDimY>());

    SplineYBuilder const builder_y(gridy);

    host_t<DFieldY> const quadrature_coeffs
            = neumann_spline_quadrature_coefficients(gridy, builder_y);
    Quadrature const integrate(quadrature_coeffs.span_cview());

    host_t<DFieldY> values(gridy);

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

} // namespace
