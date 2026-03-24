// SPDX-License-Identifier: MIT

#include <array>
#include <random>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "mesh_builder.hpp"
#include "nd_lagrange_evaluator.hpp"
#include "test_utils.hpp"
#include "view.hpp"

namespace {

struct X
{
    static constexpr bool PERIODIC = false;
};

struct Y
{
    static constexpr bool PERIODIC = false;
};

template <class T>
struct NDLagrangeNonPeriodicFixture;

template <std::size_t D, class T, bool Uniform>
struct NDLagrangeNonPeriodicFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        T,
        std::integral_constant<bool, Uniform>>> : public testing::Test
{
    struct GridX : public std::conditional_t<Uniform, UniformGridBase<X>, NonUniformGridBase<X>>
    {
    };
    struct GridY : public std::conditional_t<Uniform, UniformGridBase<Y>, NonUniformGridBase<Y>>
    {
    };
    struct LagBasisX
        : public std::conditional_t<
                  Uniform,
                  UniformLagrangeBasis<X, D, T>,
                  NonUniformLagrangeBasis<X, D, T>>
    {
    };
    struct LagBasisY
        : public std::conditional_t<
                  Uniform,
                  UniformLagrangeBasis<Y, D, T>,
                  NonUniformLagrangeBasis<Y, D, T>>
    {
    };
    // Uniform evaluation grids at different points from the coefficient grids.
    struct TestGridX : public UniformGridBase<X>
    {
    };
    struct TestGridY : public UniformGridBase<Y>
    {
    };

    using DataType = T;

    static constexpr std::size_t degree = D;
    static constexpr bool UNIFORM = Uniform;
    static constexpr double TOL = std::is_same_v<T, float> ? 5e-6 : 1e-12;
};

using degrees = std::integer_sequence<std::size_t, 2, 3, 4>;
using uniformity = std::integer_sequence<bool, true, false>;
using Cases = tuple_to_types_t<cartesian_product_t<degrees, std::tuple<double, float>, uniformity>>;

template <class Dim, class DataType, std::size_t N>
KOKKOS_FUNCTION DataType polynomial(Coord<Dim> coord, std::array<DataType, N> const& coeffs)
{
    DataType result = 0;
    for (std::size_t j = 0; j < N; ++j) {
        result += coeffs[j] * static_cast<DataType>(Kokkos::pow(static_cast<DataType>(coord), j));
    }
    return result;
}

} // namespace

TYPED_TEST_SUITE(NDLagrangeNonPeriodicFixture, Cases);

/**
 * @brief Test that NDLagrangeEvaluator reproduces a separable polynomial exactly.
 *
 * A Lagrange evaluator of degree D reproduces polynomials of degree <= D exactly.
 * The tensor-product NDLagrangeEvaluator therefore reproduces f(x,y) = p(x)*q(y)
 * when deg(p) <= D and deg(q) <= D.
 * Coefficients are provided at the knot points and evaluation is performed on a
 * separate uniform test grid to avoid trivial pass-through.
 */
TYPED_TEST(NDLagrangeNonPeriodicFixture, ExactPolynomialInterpolation)
{
    using DataType = typename TestFixture::DataType;
    using GridX = typename TestFixture::GridX;
    using GridY = typename TestFixture::GridY;
    using LagBasisX = typename TestFixture::LagBasisX;
    using LagBasisY = typename TestFixture::LagBasisY;
    using TestGridX = typename TestFixture::TestGridX;
    using TestGridY = typename TestFixture::TestGridY;

    constexpr std::size_t degree = TestFixture::degree;
    static constexpr double TOL = TestFixture::TOL;

    using EvalX = LagrangeEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            LagBasisX,
            TestGridX,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>;
    using EvalY = LagrangeEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            LagBasisY,
            TestGridY,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>;
    using Eval2D = NDLagrangeEvaluator<EvalX, EvalY>;

    using KnotGridX = typename EvalX::coeff_grid_type;
    using KnotGridY = typename EvalY::coeff_grid_type;

    // Set up the domains
    Coord<X> xmin(0.0), xmax(2.0);
    Coord<Y> ymin(0.0), ymax(3.0);
    std::size_t const ncells = 10;

    if constexpr (TestFixture::UNIFORM) {
        ddc::init_discrete_space<GridX>(GridX::init(xmin, xmax, IdxStep<GridX>(ncells + 1)));
        ddc::init_discrete_space<GridY>(GridY::init(ymin, ymax, IdxStep<GridY>(ncells + 1)));
    } else {
        ddc::init_discrete_space<GridX>(
                build_random_non_uniform_break_points(xmin, xmax, IdxStep<GridX>(ncells), 0.5));
        ddc::init_discrete_space<GridY>(
                build_random_non_uniform_break_points(ymin, ymax, IdxStep<GridY>(ncells), 0.5));
    }

    IdxRange<GridX> const x_range(Idx<GridX>(0), IdxStep<GridX>(ncells + 1));
    IdxRange<GridY> const y_range(Idx<GridY>(0), IdxStep<GridY>(ncells + 1));
    IdxRange<GridX, GridY> const idx_range(x_range, y_range);
    ddc::init_discrete_space<LagBasisX>(x_range);
    ddc::init_discrete_space<LagBasisY>(y_range);

    // Uniform evaluation grid — distinct from the coefficient knot grid
    std::size_t const ntest = 7;
    ddc::init_discrete_space<TestGridX>(TestGridX::init(xmin, xmax, IdxStep<TestGridX>(ntest)));
    ddc::init_discrete_space<TestGridY>(TestGridY::init(ymin, ymax, IdxStep<TestGridY>(ntest)));
    IdxRange<TestGridX> const test_x_range(Idx<TestGridX>(0), IdxStep<TestGridX>(ntest));
    IdxRange<TestGridY> const test_y_range(Idx<TestGridY>(0), IdxStep<TestGridY>(ntest));
    IdxRange<TestGridX, TestGridY> const test_range(test_x_range, test_y_range);

    // Random polynomial coefficients, degree D in each dimension
    std::array<DataType, degree + 1> coeffs_x, coeffs_y;
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(0.5, 1.5);
    for (std::size_t i(0); i <= degree; ++i) {
        coeffs_x[i] = static_cast<DataType>(dis(gen));
        coeffs_y[i] = static_cast<DataType>(dis(gen));
    }

    // Build 2D function f(x_i, y_j) = poly_x(x_i) * poly_y(y_j)
    using LagrangeKnotsX = std::conditional_t<
            TestFixture::UNIFORM,
            UniformLagrangeKnots<LagBasisX>,
            NonUniformLagrangeKnots<LagBasisX>>;
    using LagrangeKnotsY = std::conditional_t<
            TestFixture::UNIFORM,
            UniformLagrangeKnots<LagBasisY>,
            NonUniformLagrangeKnots<LagBasisY>>;
    IdxRange<LagrangeKnotsX, LagrangeKnotsY> lagrange_coeff_idx_range(
            ddc::discrete_space<LagBasisX>().full_domain(),
            ddc::discrete_space<LagBasisY>().full_domain());

    FieldMem<DataType, IdxRange<LagrangeKnotsX, LagrangeKnotsY>> poly_coeffs_alloc(
            lagrange_coeff_idx_range);
    Field<DataType, IdxRange<LagrangeKnotsX, LagrangeKnotsY>> poly_coeffs(poly_coeffs_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            lagrange_coeff_idx_range,
            KOKKOS_LAMBDA(Idx<LagrangeKnotsX, LagrangeKnotsY> idx) {
                poly_coeffs(idx) = polynomial(ddc::coordinate(Idx<KnotGridX>(idx)), coeffs_x)
                                   * polynomial(ddc::coordinate(Idx<KnotGridY>(idx)), coeffs_y);
            });

    // Evaluate on the test grid
    ddc::NullExtrapolationRule const null_extrap;
    EvalX lagrange_evaluator_x(null_extrap, null_extrap);
    EvalY lagrange_evaluator_y(null_extrap, null_extrap);
    Eval2D const eval2d(lagrange_evaluator_x, lagrange_evaluator_y);

    FieldMem<DataType, IdxRange<TestGridX, TestGridY>> result_alloc(test_range);
    eval2d(get_field(result_alloc), get_const_field(poly_coeffs));

    // Check against the known polynomial values
    auto const result_host = ddc::create_mirror_view_and_copy(get_const_field(result_alloc));
    ddc::host_for_each(test_range, [&](Idx<TestGridX, TestGridY> idx) {
        DataType const expected = polynomial(ddc::coordinate(Idx<TestGridX>(idx)), coeffs_x)
                                  * polynomial(ddc::coordinate(Idx<TestGridY>(idx)), coeffs_y);
        EXPECT_NEAR(
                static_cast<double>(result_host(idx)),
                static_cast<double>(expected),
                TOL * expected);
    });
}
