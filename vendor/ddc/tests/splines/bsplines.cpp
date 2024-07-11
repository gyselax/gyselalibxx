// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <array>
#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "test_utils.hpp"

template <class T>
struct BSplinesFixture;

template <std::size_t D, std::size_t Nc, bool periodic>
struct BSplinesFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        std::integral_constant<std::size_t, Nc>,
        std::integral_constant<bool, periodic>>> : public testing::Test
{
    struct DimX
    {
        static constexpr bool PERIODIC = periodic;
    };
    struct UBSplinesX : ddc::UniformBSplines<DimX, D>
    {
    };
    struct NUBSplinesX : ddc::NonUniformBSplines<DimX, D>
    {
    };
    static constexpr std::size_t spline_degree = D;
    static constexpr std::size_t ncells = Nc;
};

using degrees = std::integer_sequence<std::size_t, 1, 2, 3, 4, 5, 6>;
using ncells = std::integer_sequence<std::size_t, 10, 20, 100>;
using periodicity = std::integer_sequence<bool, true, false>;

using Cases = tuple_to_types_t<cartesian_product_t<degrees, ncells, periodicity>>;

// Trailing comma is needed to avoid spurious `gnu-zero-variadic-macro-arguments` warning with clang
TYPED_TEST_SUITE(BSplinesFixture, Cases, );

TYPED_TEST(BSplinesFixture, PartitionOfUnity_Uniform)
{
    std::size_t constexpr degree = TestFixture::spline_degree;
    using DimX = typename TestFixture::DimX;
    using BSplinesX = typename TestFixture::UBSplinesX;
    using CoordX = ddc::Coordinate<DimX>;
    static constexpr CoordX xmin = CoordX(0.0);
    static constexpr CoordX xmax = CoordX(0.2);
    static constexpr std::size_t ncells = TestFixture::ncells;
    ddc::init_discrete_space<BSplinesX>(xmin, xmax, ncells);

    std::array<double, degree + 1> values_ptr;
    std::experimental::mdspan<double, std::experimental::extents<std::size_t, degree + 1>> const
            values(values_ptr.data());

    std::size_t const n_test_points = ncells * 30;
    double const dx = (xmax - xmin) / (n_test_points - 1);

    for (std::size_t i(0); i < n_test_points; ++i) {
        CoordX test_point(xmin + dx * i);
        ddc::discrete_space<BSplinesX>().eval_basis(values, test_point);
        double sum = 0.0;
        for (std::size_t j(0); j < degree + 1; ++j) {
            sum += values(j);
        }
        EXPECT_LE(fabs(sum - 1.0), 1.0e-15);
    }
}

TYPED_TEST(BSplinesFixture, PartitionOfUnity_NonUniform)
{
    std::size_t constexpr degree = TestFixture::spline_degree;
    using DimX = typename TestFixture::DimX;
    using BSplinesX = typename TestFixture::NUBSplinesX;
    using CoordX = ddc::Coordinate<DimX>;
    static constexpr CoordX xmin = CoordX(0.0);
    static constexpr CoordX xmax = CoordX(0.2);
    static constexpr std::size_t ncells = TestFixture::ncells;
    std::vector<CoordX> breaks(ncells + 1);
    double dx = (xmax - xmin) / ncells;
    for (std::size_t i(0); i < ncells + 1; ++i) {
        breaks[i] = CoordX(xmin + i * dx);
    }
    ddc::init_discrete_space<BSplinesX>(breaks);

    std::array<double, degree + 1> values_ptr;
    std::experimental::mdspan<double, std::experimental::extents<std::size_t, degree + 1>> const
            values(values_ptr.data());


    std::size_t n_test_points = ncells * 30;
    dx = (xmax - xmin) / (n_test_points - 1);

    for (std::size_t i(0); i < n_test_points; ++i) {
        CoordX test_point(xmin + dx * i);
        ddc::discrete_space<BSplinesX>().eval_basis(values, test_point);
        double sum = 0.0;
        for (std::size_t j(0); j < degree + 1; ++j) {
            sum += values(j);
        }
        EXPECT_LE(fabs(sum - 1.0), 1.0e-15);
    }
}
