#include <array>
#include <random>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "identity_interpolation_builder.hpp"
#include "l_norm_tools.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "lagrange_evaluator.hpp"
#include "mesh_builder.hpp"
#include "test_utils.hpp"
#include "view.hpp"

namespace {

struct X
{
    static constexpr bool PERIODIC = false;
};

struct Y
{
    static constexpr bool PERIODIC = true;
};

template <class T>
struct LagrangeNonPeriodicEvaluatorFixture;

template <std::size_t D, class T, bool Uniform>
struct LagrangeNonPeriodicEvaluatorFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        T,
        std::integral_constant<bool, Uniform>>> : public testing::Test
{
    struct GridX : public std::conditional_t<Uniform, UniformGridBase<X>, NonUniformGridBase<X>>
    {
    };
    struct LagBasis
        : public std::conditional_t<
                  Uniform,
                  UniformLagrangeBasis<X, D, T>,
                  NonUniformLagrangeBasis<X, D, T>>
    {
    };
    using DataType = T;

    static constexpr std::size_t degree = D;
    static constexpr bool UNIFORM = Uniform;

    // Replace with your actual tolerance policy
    static constexpr double TOL = std::is_same_v<T, float> ? 5e-6 : 1e-12;
};

template <class T>
struct LagrangePeriodicEvaluatorFixture;

template <std::size_t D, class T, bool Uniform>
struct LagrangePeriodicEvaluatorFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        T,
        std::integral_constant<bool, Uniform>>> : public testing::Test
{
    struct GridY : public std::conditional_t<Uniform, UniformGridBase<Y>, NonUniformGridBase<Y>>
    {
    };
    struct LagBasis
        : public std::conditional_t<
                  Uniform,
                  UniformLagrangeBasis<Y, D, T>,
                  NonUniformLagrangeBasis<Y, D, T>>
    {
    };
    using DataType = T;

    struct RefinedGridY
        : public std::conditional_t<Uniform, UniformGridBase<Y>, NonUniformGridBase<Y>>
    {
    };
    struct RefinedLagBasis
        : public std::conditional_t<
                  Uniform,
                  UniformLagrangeBasis<Y, D, T>,
                  NonUniformLagrangeBasis<Y, D, T>>
    {
    };

    struct TestGrid : public UniformGridBase<Y>
    {
    };

    static constexpr std::size_t degree = D;
    static constexpr bool UNIFORM = Uniform;
};

using degrees = std::integer_sequence<std::size_t, 2, 3, 4>;
using uniformity = std::integer_sequence<bool, true, false>;
using Cases = tuple_to_types_t<cartesian_product_t<degrees, std::tuple<double, float>, uniformity>>;

template <class X, class DataType, std::size_t N>
DataType polynomial(Coord<X> coord, std::array<DataType, N> coeffs)
{
    DataType result = 0.0;
    for (int j(0); j < N; ++j) {
        result += coeffs[j] * std::pow(coord, j);
    }
    return result;
}

} // namespace

TYPED_TEST_SUITE(LagrangeNonPeriodicEvaluatorFixture, Cases);

TYPED_TEST_SUITE(LagrangePeriodicEvaluatorFixture, Cases);

TYPED_TEST(LagrangeNonPeriodicEvaluatorFixture, ExactPolynomialInterpolation)
{
    using DataType = typename TestFixture::DataType;
    using GridX = TestFixture::GridX;
    using LagBasis = TestFixture::LagBasis;
    using Builder = IdentityInterpolationBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            GridX,
            LagBasis>;
    using Evaluator = LagrangeEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            LagBasis,
            GridX,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>;
    using IdxRangeCoeff = typename Builder::template batched_basis_domain_type<IdxRange<GridX>>;

    constexpr std::size_t degree = TestFixture::degree;
    static constexpr double TOL = TestFixture::TOL;

    Coord<X> xmin(0);
    Coord<X> xmax(2);
    std::size_t ncells(10);
    if constexpr (TestFixture::UNIFORM) {
        ddc::init_discrete_space<GridX>(GridX::init(xmin, xmax, IdxStep<GridX>(ncells + 1)));
    } else {
        std::vector<Coord<X>> points
                = build_random_non_uniform_break_points(xmin, xmax, IdxStep<GridX>(ncells), 0.5);
        ddc::init_discrete_space<GridX>(points);
    }
    IdxRange<GridX> idx_range(Idx<GridX>(0), IdxStep<GridX>(ncells + 1));
    ddc::init_discrete_space<LagBasis>(idx_range);

    std::array<DataType, degree + 1> coeffs;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution dis(0., 1.);
    for (int i(0); i < degree + 1; ++i) {
        coeffs[i] = dis(gen);
    }

    host_t<FieldMem<DataType, IdxRange<GridX>>> function_values_host(idx_range);
    ddc::host_for_each(idx_range, [&](Idx<GridX> i) {
        function_values_host(i) = polynomial(ddc::coordinate(i), coeffs);
    });

    auto function_values_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(function_values_host));
    Field<DataType, IdxRange<GridX>> function_values(function_values_alloc);

    Builder builder;
    FieldMem<DataType, IdxRangeCoeff> lagrange_coeffs_alloc(
            ddc::discrete_space<LagBasis>().full_domain());
    Field<DataType, IdxRangeCoeff> lagrange_coeffs(lagrange_coeffs_alloc);

    ddc::NullExtrapolationRule extrapol;
    Evaluator evaluator(extrapol, extrapol);

    host_t<FieldMem<Coord<X>, IdxRange<GridX>>> test_coords_host_alloc(idx_range);

    ddc::host_for_each(get_idx_range(test_coords_host_alloc), [&](Idx<GridX> idx) {
        test_coords_host_alloc(idx) = xmin + (xmax - xmin) * dis(gen);
    });

    auto test_coords_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(test_coords_host_alloc));
    ConstField<Coord<X>, IdxRange<GridX>> test_coords = get_const_field(test_coords_alloc);

    builder(lagrange_coeffs, get_const_field(function_values));
    evaluator(function_values, get_const_field(test_coords), get_const_field(lagrange_coeffs));

    ddc::parallel_deepcopy(get_field(function_values_host), get_const_field(function_values));

    ddc::host_for_each(get_idx_range(test_coords_host_alloc), [&](Idx<GridX> idx) {
        EXPECT_NEAR(
                function_values_host(idx),
                polynomial(test_coords_host_alloc(idx), coeffs),
                TOL);
    });
}

template <class DataType, class LagBasis, class GridType, class TestGridType>
DataType get_cosine_error(IdxRange<GridType> idx_range, IdxRange<TestGridType> test_idx_range)
{
    using Dim = typename GridType::continuous_dimension_type;
    using Builder = IdentityInterpolationBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            GridType,
            LagBasis>;
    using Evaluator = LagrangeEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            LagBasis,
            TestGridType,
            ddc::PeriodicExtrapolationRule<Dim>,
            ddc::PeriodicExtrapolationRule<Dim>>;
    using IdxRangeCoeff = typename Builder::template batched_basis_domain_type<IdxRange<GridType>>;

    FieldMem<DataType, IdxRange<GridType>> function_values_alloc(idx_range);
    Field<DataType, IdxRange<GridType>> function_values(function_values_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(Idx<GridType> i) {
                function_values(i) = Kokkos::cos(ddc::coordinate(i));
            });

    Builder builder;
    FieldMem<DataType, IdxRangeCoeff> lagrange_coeffs_alloc(
            ddc::discrete_space<LagBasis>().full_domain());
    Field<DataType, IdxRangeCoeff> lagrange_coeffs(lagrange_coeffs_alloc);

    ddc::PeriodicExtrapolationRule<Dim> extrapol;
    Evaluator evaluator(extrapol, extrapol);

    builder(lagrange_coeffs, get_const_field(function_values));

    FieldMem<DataType, IdxRange<TestGridType>> values_at_test_points_alloc(test_idx_range);
    Field<DataType, IdxRange<TestGridType>> values_at_test_points(values_at_test_points_alloc);
    evaluator(values_at_test_points, get_const_field(lagrange_coeffs));

    return error_norm_inf(
            Kokkos::DefaultExecutionSpace(),
            get_const_field(values_at_test_points),
            KOKKOS_LAMBDA(Idx<TestGridType> idx) {
                return static_cast<DataType>(Kokkos::cos(ddc::coordinate(idx)));
            });
}

TYPED_TEST(LagrangePeriodicEvaluatorFixture, Convergence)
{
    using DataType = TestFixture::DataType;
    using GridY = TestFixture::GridY;
    using LagBasis = TestFixture::LagBasis;
    using RefinedGridY = TestFixture::RefinedGridY;
    using RefinedLagBasis = TestFixture::RefinedLagBasis;
    using TestGrid = TestFixture::TestGrid;

    constexpr std::size_t degree = TestFixture::degree;

    Coord<Y> ymin(0);
    Coord<Y> ymax(2 * Kokkos::numbers::pi);
    std::size_t ncells(10);
    if constexpr (TestFixture::UNIFORM) {
        ddc::init_discrete_space<GridY>(GridY::init(ymin, ymax, IdxStep<GridY>(ncells + 1)));
        ddc::init_discrete_space<RefinedGridY>(
                RefinedGridY::init(ymin, ymax, IdxStep<RefinedGridY>(2 * ncells + 1)));
    } else {
        std::vector<Coord<Y>> points
                = build_random_non_uniform_break_points(ymin, ymax, IdxStep<GridY>(ncells), 0.1);
        ddc::init_discrete_space<GridY>(points);
        std::vector<Coord<Y>> refined_points(ncells * 2 + 1);
        for (int i(0); i < ncells; ++i) {
            refined_points[2 * i] = points[i];
            refined_points[2 * i + 1] = 0.5 * (points[i] + points[i + 1]);
        }
        refined_points.back() = points.back();
        ddc::init_discrete_space<RefinedGridY>(refined_points);
    }
    IdxRange<GridY> idx_range(Idx<GridY>(0), IdxStep<GridY>(ncells + 1));
    ddc::init_discrete_space<LagBasis>(idx_range);

    IdxRange<RefinedGridY>
            refined_idx_range(Idx<RefinedGridY>(0), IdxStep<RefinedGridY>(2 * ncells + 1));
    ddc::init_discrete_space<RefinedLagBasis>(refined_idx_range);

    ddc::init_discrete_space<TestGrid>(TestGrid::init(ymin, ymax, IdxStep<TestGrid>(6 * ncells)));

    IdxRange<TestGrid> test_idx_range(Idx<TestGrid>(0), IdxStep<TestGrid>(6 * ncells));

    DataType cosine_error = get_cosine_error<DataType, LagBasis>(
            idx_range.remove_last(IdxStep<GridY>(1)), // Remove repeat periodic point
            test_idx_range);
    DataType refined_cosine_error = get_cosine_error<DataType, RefinedLagBasis>(
            refined_idx_range.remove_last(IdxStep<RefinedGridY>(1)), // Remove repeat periodic point
            test_idx_range);

    DataType order = std::log(cosine_error / refined_cosine_error) / std::log(2.0);

    std::cout << cosine_error << " " << refined_cosine_error << " " << order << std::endl;
    if constexpr (TestFixture::UNIFORM) {
        EXPECT_NEAR(order, degree + 1, 0.5);
    } else {
        EXPECT_GT(order, degree);
    }
}
