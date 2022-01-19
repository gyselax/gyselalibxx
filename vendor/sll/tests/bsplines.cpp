#include <array>
#include <cmath>

#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/UniformDiscretization>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/view.hpp>

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
    static constexpr std::size_t spline_degree = D;
    static constexpr std::size_t ncells = Nc;
};

using degrees = std::integer_sequence<std::size_t, 1, 2, 3, 4, 5, 6>;
using ncells = std::integer_sequence<std::size_t, 10, 20, 100>;
using periodicity = std::integer_sequence<bool, true, false>;

using Cases = tuple_to_types_t<cartesian_product_t<degrees, ncells, periodicity>>;

TYPED_TEST_SUITE(BSplinesFixture, Cases);

TYPED_TEST(BSplinesFixture, PartitionOfUnity_Uniform)
{
    std::size_t constexpr degree = TestFixture::spline_degree;
    using DimX = typename TestFixture::DimX;
    using CoordX = Coordinate<DimX>;
    static constexpr CoordX xmin = CoordX(0.0);
    static constexpr CoordX xmax = CoordX(0.2);
    static constexpr std::size_t ncells = TestFixture::ncells;
    UniformBSplines<DimX, degree> const bsplines {xmin, xmax, ncells};

    std::array<double, degree + 1> vals_data;
    DSpan1D values(vals_data.data(), degree + 1);

    std::size_t const n_test_points = ncells * 30;
    double const dx = (xmax - xmin) / (n_test_points - 1);

    int jmin;
    for (std::size_t i(0); i < n_test_points; ++i) {
        double const test_point = xmin + dx * i;
        bsplines.eval_basis(values, jmin, test_point);
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
    using CoordX = Coordinate<DimX>;
    static constexpr CoordX xmin = CoordX(0.0);
    static constexpr CoordX xmax = CoordX(0.2);
    static constexpr std::size_t ncells = TestFixture::ncells;
    std::vector<double> breaks(ncells + 1);
    double dx = (xmax - xmin) / ncells;
    for (std::size_t i(0); i < ncells + 1; ++i) {
        breaks[i] = xmin + i * dx;
    }
    NonUniformBSplines<DimX, degree> const bsplines {breaks};

    std::array<double, degree + 1> vals_data;
    DSpan1D values(vals_data.data(), degree + 1);

    std::size_t n_test_points = ncells * 30;
    dx = (xmax - xmin) / (n_test_points - 1);

    int jmin;
    for (std::size_t i(0); i < n_test_points; ++i) {
        double test_point = xmin + dx * i;
        bsplines.eval_basis(values, jmin, test_point);
        double sum = 0.0;
        for (std::size_t j(0); j < degree + 1; ++j) {
            sum += values(j);
        }
        EXPECT_LE(fabs(sum - 1.0), 1.0e-15);
    }
}
