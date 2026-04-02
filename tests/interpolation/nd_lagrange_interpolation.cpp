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
#include "nd_identity_interpolation_builder.hpp"
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

struct XPeriodic
{
    static constexpr bool PERIODIC = true;
};

template <class T>
struct NDLagrangeNonPeriodicFixture;

template <class T>
struct NDLagrangePeriodicFixture;

template <std::size_t D, class T, bool Uniform, class X>
struct NDLagrangeFixture : public testing::Test
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

template <std::size_t D, class T, bool Uniform>
struct NDLagrangeNonPeriodicFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        T,
        std::integral_constant<bool, Uniform>>> : public NDLagrangeFixture<D, T, Uniform, X>
{
};

template <std::size_t D, class T, bool Uniform>
struct NDLagrangePeriodicFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        T,
        std::integral_constant<bool, Uniform>>> : public NDLagrangeFixture<D, T, Uniform, XPeriodic>
{
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

template <class Dim, class DataType, std::size_t N>
KOKKOS_FUNCTION DataType polynomial_deriv(Coord<Dim> coord, std::array<DataType, N> const& coeffs)
{
    DataType result = 0;
    for (std::size_t j = 1; j < N; ++j) {
        result += static_cast<DataType>(j) * coeffs[j]
                  * static_cast<DataType>(Kokkos::pow(static_cast<DataType>(coord), j - 1));
    }
    return result;
}

template <class DataType, class Dim1, class Dim2, std::size_t N>
void fill_polynomial_2d(
        Field<DataType, IdxRange<Dim1, Dim2>> poly_coeffs,
        std::array<DataType, N> const& coeffs1,
        std::array<DataType, N> const& coeffs2)
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(poly_coeffs),
            KOKKOS_LAMBDA(Idx<Dim1, Dim2> idx) {
                poly_coeffs(idx) = polynomial(ddc::coordinate(Idx<Dim1>(idx)), coeffs1)
                                   * polynomial(ddc::coordinate(Idx<Dim2>(idx)), coeffs2);
            });
}

} // namespace

//TYPED_TEST_SUITE(NDLagrangeNonPeriodicFixture, Cases);
TYPED_TEST_SUITE(NDLagrangePeriodicFixture, Cases);

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
    using Builder = NDIdentityInterpolationBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            IdxRange<GridX, GridY>,
            IdxRange<LagBasisX, LagBasisY>>;

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
    FieldMem<DataType, IdxRange<GridX, GridY>> vals_alloc("vals", idx_range);
    fill_polynomial_2d(get_field(vals_alloc), coeffs_x, coeffs_y);
    Builder const builder;
    using CoeffIdxRange = InterpolationBuilderTraits<
            Builder>::template batched_basis_idx_range_type<IdxRange<GridX, GridY>>;
    FieldMem<DataType, CoeffIdxRange>
            poly_coeffs_alloc("coeffs", batched_basis_idx_range(builder, idx_range));
    builder(get_field(poly_coeffs_alloc), get_const_field(vals_alloc));

    // Evaluate on the test grid
    ddc::NullExtrapolationRule const null_extrap;
    EvalX lagrange_evaluator_x(null_extrap, null_extrap);
    EvalY lagrange_evaluator_y(null_extrap, null_extrap);
    Eval2D const eval2d(lagrange_evaluator_x, lagrange_evaluator_y);

    FieldMem<DataType, IdxRange<TestGridX, TestGridY>> result_alloc(test_range);
    eval2d(get_field(result_alloc), get_const_field(poly_coeffs_alloc));

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

/**
 * @brief Test that NDLagrangeEvaluator reproduces partial and cross-derivatives of a separable
 * polynomial exactly.
 *
 * For f(x,y) = p(x)*q(y):
 *   - d/dx f = p'(x)*q(y)
 *   - d/dy f = p(x)*q'(y)
 *   - d²/dxdy f = p'(x)*q'(y)
 *
 * Each of these is a polynomial of lower degree that the degree-D Lagrange evaluator reproduces
 * exactly. We therefore expect only floating-point roundoff error.
 */
TYPED_TEST(NDLagrangeNonPeriodicFixture, ExactDerivatives)
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

    using LagrangeKnotsX = std::conditional_t<
            TestFixture::UNIFORM,
            UniformLagrangeKnots<LagBasisX>,
            NonUniformLagrangeKnots<LagBasisX>>;
    using LagrangeKnotsY = std::conditional_t<
            TestFixture::UNIFORM,
            UniformLagrangeKnots<LagBasisY>,
            NonUniformLagrangeKnots<LagBasisY>>;

    Coord<X> xmin(0.0), xmax(1.0);
    Coord<Y> ymin(0.0), ymax(1.0);
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
    ddc::init_discrete_space<LagBasisX>(x_range);
    ddc::init_discrete_space<LagBasisY>(y_range);

    std::size_t const ntest = 7;
    ddc::init_discrete_space<TestGridX>(TestGridX::init(xmin, xmax, IdxStep<TestGridX>(ntest)));
    ddc::init_discrete_space<TestGridY>(TestGridY::init(ymin, ymax, IdxStep<TestGridY>(ntest)));
    IdxRange<TestGridX> const test_x_range(Idx<TestGridX>(0), IdxStep<TestGridX>(ntest));
    IdxRange<TestGridY> const test_y_range(Idx<TestGridY>(0), IdxStep<TestGridY>(ntest));
    IdxRange<TestGridX, TestGridY> const test_range(test_x_range, test_y_range);

    std::array<DataType, degree + 1> coeffs_x, coeffs_y;
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (std::size_t i(0); i < degree + 1; ++i) {
        coeffs_x[i] = static_cast<DataType>(dis(gen));
        coeffs_y[i] = static_cast<DataType>(dis(gen));
    }

    IdxRange<LagrangeKnotsX, LagrangeKnotsY> const lagrange_coeff_idx_range(
            ddc::discrete_space<LagBasisX>().full_domain(),
            ddc::discrete_space<LagBasisY>().full_domain());
    FieldMem<DataType, IdxRange<LagrangeKnotsX, LagrangeKnotsY>> poly_coeffs_alloc(
            lagrange_coeff_idx_range);
    fill_polynomial_2d(get_field(poly_coeffs_alloc), coeffs_x, coeffs_y);

    ddc::NullExtrapolationRule const null_extrap;
    Eval2D const eval2d(EvalX(null_extrap, null_extrap), EvalY(null_extrap, null_extrap));

    // Derivative tolerances scale with cell size (derivative basis values scale as 1/dx^n)
    double const dx_max = ddcHelper::maximum_distance_between_adjacent_points(x_range);
    double const dy_max = ddcHelper::maximum_distance_between_adjacent_points(y_range);
    double const tol_dx = TOL * (degree + 1) / dx_max;
    double const tol_dy = TOL * (degree + 1) / dy_max;
    double const tol_dxdy = TOL * (degree + 1) * (degree + 1) / (dx_max * dy_max);

    FieldMem<DataType, IdxRange<TestGridX, TestGridY>> deriv_alloc(test_range);
    auto deriv_host = ddc::create_mirror_view(get_field(deriv_alloc));

    // d/dx f = p'(x) * q(y)
    Idx<ddc::Deriv<X>> dx(1);
    eval2d.deriv(dx, get_field(deriv_alloc), get_const_field(poly_coeffs_alloc));
    ddc::parallel_deepcopy(get_field(deriv_host), get_const_field(deriv_alloc));
    ddc::host_for_each(test_range, [&](Idx<TestGridX, TestGridY> idx) {
        double const expected = polynomial_deriv(ddc::coordinate(Idx<TestGridX>(idx)), coeffs_x)
                                * polynomial(ddc::coordinate(Idx<TestGridY>(idx)), coeffs_y);
        EXPECT_NEAR(static_cast<double>(deriv_host(idx)), expected, tol_dx);
    });

    // d/dy f = p(x) * q'(y)
    Idx<ddc::Deriv<Y>> dy(1);
    eval2d.deriv(dy, get_field(deriv_alloc), get_const_field(poly_coeffs_alloc));
    ddc::parallel_deepcopy(get_field(deriv_host), get_const_field(deriv_alloc));
    ddc::host_for_each(test_range, [&](Idx<TestGridX, TestGridY> idx) {
        double const expected = polynomial(ddc::coordinate(Idx<TestGridX>(idx)), coeffs_x)
                                * polynomial_deriv(ddc::coordinate(Idx<TestGridY>(idx)), coeffs_y);
        EXPECT_NEAR(static_cast<double>(deriv_host(idx)), expected, tol_dy);
    });

    // d²/dxdy f = p'(x) * q'(y)
    Idx<ddc::Deriv<X>, ddc::Deriv<Y>> dxdy(1, 1);
    eval2d.deriv(dxdy, get_field(deriv_alloc), get_const_field(poly_coeffs_alloc));
    ddc::parallel_deepcopy(get_field(deriv_host), get_const_field(deriv_alloc));
    ddc::host_for_each(test_range, [&](Idx<TestGridX, TestGridY> idx) {
        double const expected = polynomial_deriv(ddc::coordinate(Idx<TestGridX>(idx)), coeffs_x)
                                * polynomial_deriv(ddc::coordinate(Idx<TestGridY>(idx)), coeffs_y);
        EXPECT_NEAR(static_cast<double>(deriv_host(idx)), expected, tol_dxdy);
    });
}

/**
 * @brief Test that NDIdentityInterpolationBuilder correctly fills wrap-around coefficients
 * for a periodic dimension.
 *
 * The test uses a degree-3 uniform Lagrange basis on XPeriodic × Y (non-periodic).
 * The periodic dimension runs over [0, 2π] with ncells+1 nodes; f(x, y) = cos(x) * (1 + y).
 *
 * The builder fills:
 *   - coeffs[0..ncells-1, j]            from vals (main body, avoiding the duplicate at 2π)
 *   - coeffs[ncells..ncells+degree-1, j] as wrap-around = coeffs[0..degree-1, j]
 *
 * Two checks:
 *  1. Wrap-around: each wrap-around coefficient equals the corresponding coefficient from
 *     the start of the periodic domain.
 *  2. Round-trip: evaluating the built coefficients on a test grid inside (0, 2π) × (0, 1)
 *     reproduces cos(x)*(1+y) within the Lagrange interpolation error.
 */
TYPED_TEST(NDLagrangePeriodicFixture, PeriodicWraparound)
{
    using DataType = typename TestFixture::DataType;
    using GridX = typename TestFixture::GridX;
    using GridY = typename TestFixture::GridY;
    using LagBasisX = typename TestFixture::LagBasisX;
    using LagBasisY = typename TestFixture::LagBasisY;
    using LagKnotsX = std::conditional_t<
            TestFixture::UNIFORM,
            UniformLagrangeKnots<LagBasisX>,
            NonUniformLagrangeKnots<LagBasisX>>;
    using LagKnotsY = std::conditional_t<
            TestFixture::UNIFORM,
            UniformLagrangeKnots<LagBasisY>,
            NonUniformLagrangeKnots<LagBasisY>>;
    using TestGridX = typename TestFixture::TestGridX;
    using TestGridY = typename TestFixture::TestGridY;
    using Builder = NDIdentityInterpolationBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            IdxRange<GridX, GridY>,
            IdxRange<LagBasisX, LagBasisY>>;
    using EvalX = LagrangeEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            LagBasisX,
            TestGridX,
            ddc::PeriodicExtrapolationRule<XPeriodic>,
            ddc::PeriodicExtrapolationRule<XPeriodic>>;
    using EvalY = LagrangeEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            LagBasisY,
            TestGridY,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>;
    using Eval2D = NDLagrangeEvaluator<EvalX, EvalY>;

    Coord<XPeriodic> const xmin(0.0), xmax(2.0 * M_PI);
    Coord<Y> const ymin(0.0), ymax(1.0);
    std::size_t const ncells = 10;
    constexpr std::size_t degree = TestFixture::degree;

    if constexpr (TestFixture::UNIFORM) {
        ddc::init_discrete_space<GridX>(GridX::init(xmin, xmax, IdxStep<GridX>(ncells + 1)));
        ddc::init_discrete_space<GridY>(GridY::init(ymin, ymax, IdxStep<GridY>(ncells + 1)));
    } else {
        ddc::init_discrete_space<GridX>(
                build_random_non_uniform_break_points(xmin, xmax, IdxStep<GridX>(ncells), 0.5));
        ddc::init_discrete_space<GridY>(
                build_random_non_uniform_break_points(ymin, ymax, IdxStep<GridY>(ncells), 0.5));
    }

    IdxRange<GridX> const knot_x_range(Idx<GridX>(0), IdxStep<GridX>(ncells + 1));
    IdxRange<GridY> const knot_y_range(Idx<GridY>(0), IdxStep<GridY>(ncells + 1));
    IdxRange<GridX, GridY> const
            idx_range(knot_x_range.remove_last(IdxStep<GridX>(1)), knot_y_range);

    ddc::init_discrete_space<LagBasisX>(knot_x_range);
    ddc::init_discrete_space<LagBasisY>(knot_y_range);

    std::size_t const ntest = 7;
    ddc::init_discrete_space<TestGridX>(TestGridX::init(xmin, xmax, IdxStep<TestGridX>(ntest)));
    ddc::init_discrete_space<TestGridY>(TestGridY::init(ymin, ymax, IdxStep<TestGridY>(ntest)));
    IdxRange<TestGridX> const test_x_range(Idx<TestGridX>(0), IdxStep<TestGridX>(ntest));
    IdxRange<TestGridY> const test_y_range(Idx<TestGridY>(0), IdxStep<TestGridY>(ntest));
    IdxRange<TestGridX, TestGridY> const test_range(test_x_range, test_y_range);

    // Function values on the interpolation mesh: f(x, y) = cos(x) * (1 + y)
    FieldMem<DataType, IdxRange<GridX, GridY>> vals_alloc("vals", idx_range);
    Field<DataType, IdxRange<GridX, GridY>> vals(vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(Idx<GridX, GridY> idx) {
                DataType const x = ddc::coordinate(Idx<GridX>(idx));
                DataType const y = ddc::coordinate(Idx<GridY>(idx));
                vals(idx) = Kokkos::cos(x) * (1.0 + y);
            });

    // Build coefficients via the identity builder
    Builder const builder;
    using CoeffIdxRange = InterpolationBuilderTraits<
            Builder>::template batched_basis_idx_range_type<IdxRange<GridX, GridY>>;
    static_assert(std::is_same_v<IdxRange<LagKnotsX, LagKnotsY>, CoeffIdxRange>);
    CoeffIdxRange knot_idx_range(batched_basis_idx_range(builder, idx_range));
    FieldMem<DataType, CoeffIdxRange> coeffs_alloc("coeffs", knot_idx_range);
    builder(get_field(coeffs_alloc), get_const_field(vals_alloc));

    // Check 1: wrap-around coefficients.
    // For each j in Y and each i in [0, degree], coeffs[ncells + i, j] must equal
    // coeffs[i, j] because the periodic wrap-around copies the first (degree+1) values.
    auto const coeffs_host = ddc::create_mirror_view_and_copy(get_const_field(coeffs_alloc));
    IdxRange<LagKnotsX> periodic_knot_idx_range(knot_idx_range);
    IdxRange<LagKnotsY> non_periodic_knot_idx_range(knot_idx_range);
    for (std::size_t i = 0; i <= degree; ++i) {
        Idx<LagKnotsX> const main_idx(periodic_knot_idx_range.front() + i);
        Idx<LagKnotsX> const wrap_idx(periodic_knot_idx_range.back() - degree + i);
        ddc::host_for_each(non_periodic_knot_idx_range, [&](Idx<LagKnotsY> idx_y) {
            EXPECT_DOUBLE_EQ(coeffs_host(main_idx, idx_y), coeffs_host(wrap_idx, idx_y));
        });
    }

    // Check 2: round-trip evaluation inside the domain matches cos(x)*(1+y).
    // Degree-3 Lagrange on 10 cells: error is O(h^4) with h = 2π/10 ≈ 0.63, so ~0.16.
    // Use a loose tolerance to test correctness of the coefficient layout, not accuracy.
    double const tol = 0.2;
    FieldMem<DataType, IdxRange<TestGridX, TestGridY>> result_alloc(test_range);
    ddc::PeriodicExtrapolationRule<XPeriodic> const periodic_extrap;
    ddc::NullExtrapolationRule const null_extrap;
    EvalX lagrange_evaluator_x(periodic_extrap, periodic_extrap);
    EvalY lagrange_evaluator_y(null_extrap, null_extrap);
    Eval2D const eval2d(lagrange_evaluator_x, lagrange_evaluator_y);
    eval2d(get_field(result_alloc), get_const_field(coeffs_alloc));

    auto const result_host = ddc::create_mirror_view_and_copy(get_const_field(result_alloc));
    ddc::host_for_each(test_range, [&](Idx<TestGridX, TestGridY> idx) {
        double const x = ddc::coordinate(Idx<TestGridX>(idx));
        double const y = ddc::coordinate(Idx<TestGridY>(idx));
        double const expected = std::cos(x) * (1.0 + y);
        EXPECT_NEAR(result_host(idx), expected, tol);
    });
}
