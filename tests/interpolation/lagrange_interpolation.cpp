#include <array>
#include <random>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "identity_interpolation_builder.hpp"
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

template <class T>
struct LagrangeEvaluatorFixture;

template <std::size_t D, class T, bool Uniform>
struct LagrangeEvaluatorFixture<std::tuple<
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
                  UniformLagrangeBasis<GridX, D, T>,
                  NonUniformLagrangeBasis<GridX, D, T>>
    {
    };
    using DataType = T;

    static constexpr std::size_t degree = D;
    static constexpr bool UNIFORM = Uniform;

    // Replace with your actual tolerance policy
    static constexpr double TOL = std::is_same_v<T, float> ? 1e-6 : 1e-12;
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

TYPED_TEST_SUITE(LagrangeEvaluatorFixture, Cases);

TYPED_TEST(LagrangeEvaluatorFixture, ExactPolynomialInterpolation)
{
    static constexpr bool Periodic = X::PERIODIC;
    using DataType = typename TestFixture::DataType;
    using GridX = TestFixture::GridX;
    using LagBasis = TestFixture::LagBasis;
    using Builder = IdentityInterpolationBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            GridX,
            LagBasis>;
    using Evaluator = LagrangeEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DataType,
            LagBasis,
            GridX,
            std::conditional_t<
                    Periodic,
                    ddc::PeriodicExtrapolationRule<X>,
                    ddc::NullExtrapolationRule>,
            std::conditional_t<
                    Periodic,
                    ddc::PeriodicExtrapolationRule<X>,
                    ddc::NullExtrapolationRule>>;
    using CoeffIdxRange = typename Builder::template batched_basis_domain_type<IdxRange<GridX>>;

    constexpr std::size_t degree = TestFixture::degree;
    static constexpr double TOL = TestFixture::TOL;

    Coord<X> xmin(0);
    Coord<X> xmax(2);
    std::size_t ncells(10);
    if constexpr (TestFixture::UNIFORM) {
        ddc::init_discrete_space<GridX>(GridX::init(xmin, xmax, IdxStep<GridX>(ncells + 1)));
    } else {
        std::vector<Coord<X>> points = build_random_non_uniform_break_points(
                xmin,
                xmax,
                IdxStep<GridX>(ncells),
                0.0); //5);
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
    FieldMem<DataType, CoeffIdxRange> lagrange_coeffs_alloc(
            ddc::discrete_space<LagBasis>().full_domain());
    Field<DataType, CoeffIdxRange> lagrange_coeffs(lagrange_coeffs_alloc);

    std::conditional_t<Periodic, ddc::PeriodicExtrapolationRule<X>, ddc::NullExtrapolationRule>
            extrapol;
    Evaluator evaluator(extrapol, extrapol);

    host_t<FieldMem<Coord<X>, IdxRange<GridX>>> test_coords_host_alloc(idx_range);

    ddc::host_for_each(get_idx_range(test_coords_host_alloc), [&](Idx<GridX> idx) {
        test_coords_host_alloc(idx) = xmin + (xmax - xmin) * dis(gen);
    });

    auto test_coords_alloc = ddc::create_mirror_view_and_copy(get_field(test_coords_host_alloc));
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

