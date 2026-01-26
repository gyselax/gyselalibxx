#include <array>
#include <random>

#include <gtest/gtest.h>

#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "mesh_builder.hpp"
#include "test_utils.hpp"
#include "view.hpp"

namespace {

template <class T>
struct LagrangeBasisFixture;

template <std::size_t D, class T, bool Uniform, bool Periodic>
struct LagrangeBasisFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        T,
        std::integral_constant<bool, Uniform>,
        std::integral_constant<bool, Periodic>>> : public testing::Test
{
    struct X
    {
        static constexpr bool PERIODIC = Periodic;
    };
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
    static constexpr bool PERIODIC = Periodic;

    // Replace with your actual tolerance policy
    static constexpr double TOL = std::is_same_v<T, float> ? 1e-6 : 1e-12;
};

using degrees = std::integer_sequence<std::size_t, 2, 3, 4>;
using uniformity = std::integer_sequence<bool, true, false>;
using periodicity = std::integer_sequence<bool, true, false>;
using Cases = tuple_to_types_t<
        cartesian_product_t<degrees, std::tuple<double, float>, uniformity, periodicity>>;

} // namespace

TYPED_TEST_SUITE(LagrangeBasisFixture, Cases);

TYPED_TEST(LagrangeEvaluatorFixture, ExactPolynomialInterpolation)
{
    static constexpr bool Periodic = TestFixture::Periodic;
    using DataType = typename TestFixture::DataType;
    using Builder = IdentityInterpolationBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            GridX,
            LagBasis>;
    using Evaluator = LagrangeEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            LagBasis,
            GridX,
            std::conditional_t<
                    Periodic,
                    ddc::PeriodicExtrapolationRule<X>,
                    ddc::NullExtrapolationRule<X>>,
            std::conditional_t<
                    Periodic,
                    ddc::PeriodicExtrapolationRule<X>,
                    ddc::NullExtrapolationRule<X>>>;
    using CoeffIdxRange = Builder::batched_basis_domain_type<GridX>;

    constexpr std::size_t degree = TestFixture::degree;

    Coord<X> xmin(0);
    Coord<X> xmax(2);
    std::size_t ncells(20);
    if constexpr (TestFixture::UNIFORM) {
        ddc::init_discrete_space<GridX>(GridX::init(xmin, xmax, IdxStep<GridX>(ncells)));
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
    for (int i(0); i < Degree + 1; ++i) {
        coeffs[i] = dis(gen);
    }

    host_t<FieldMem<IdxRange<GridX>>> function_values_host(idx_range);
    ddc::for_each(idx_range, [&](Idx<GridX> i) {
        function_values(i) = 0.0;
        DataType x = ddc::coordinate(i);
        for (int i(0); i < Degree + 1; ++i) {
            values(i) += coeffs[i] * std::pow(x, i);
        }
    });

    FieldMem<DataType, IdxRange<GridX>> function_values_alloc(idx_range);
    Field<DataType, IdxRange<GridX>> function_values(function_values_alloc);

    ddc::parallel_deepcopy(function_values, get_const_field(function_values_host));

    FieldMem<DataType, CoeffIdxRange> lagrange_coeffs_alloc(idx_range);
    Field<DataType, CoeffIdxRange> lagrange_coeffs(lagrange_coeffs_alloc);
    Builder builder;

    std::conditional_t<Periodic, ddc::PeriodicExtrapolationRule<X>, ddc::NullExtrapolationRule<X>>
            extrapol;
    Evaluator evaluator(extrapol, extrapol);

    FieldMem<Coord<X>, IdxRange<GridX>> test_coords_alloc(idx_range);
    Field<Coord<X>, IdxRange<GridX>> test_coords(test_coords_alloc);

    choose_coords(test_coords);

    builder(lagrange_coeffs, get_const_field(function_values));
    evaluator(function_values, test_coords, get_const_field(lagrange_coeffs));

    for (int k = 0; k < 20; ++k) {
        DataType x = TestFixture::random_point_inside_domain();
        EXPECT_NEAR(evaluator(x), f(x), kTol);
    }
}

