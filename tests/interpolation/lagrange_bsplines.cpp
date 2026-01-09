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
        static constexpr bool PERIODIC = false;
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

TYPED_TEST(LagrangeBasisFixture, KroneckerDeltaAtKnots)
{
    using X = typename TestFixture::X;
    using GridX = typename TestFixture::GridX;
    using DataType = typename TestFixture::DataType;
    using LagBasis = typename TestFixture::LagBasis;
    using knot_discrete_dimension_type = std::conditional_t<
            TestFixture::UNIFORM,
            UniformLagrangeKnots<LagBasis>,
            NonUniformLagrangeKnots<LagBasis>>;
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

    constexpr Idx<knot_discrete_dimension_type> poly_start(0);

    std::array<DataType, degree + 1> basis_storage;
    Span1D<DataType> basis_values(basis_storage.data(), degree + 1);
    std::cout << basis_values.size() << " " << degree + 1 << std::endl;

    // --- Test ---
    for (std::size_t j = 0; j < degree + 1; ++j) {
        ddc::discrete_space<LagBasis>()
                .eval_basis(basis_values, ddc::coordinate(poly_start + j), poly_start);

        for (std::size_t i = 0; i < degree + 1; ++i) {
            if (i == j) {
                EXPECT_NEAR(basis_values[i], DataType(1), TestFixture::TOL)
                        << "Basis " << i << " at its own knot";
            } else {
                EXPECT_NEAR(basis_values[i], DataType(0), TestFixture::TOL)
                        << "Basis " << i << " at knot " << j;
            }
        }
    }
}

TYPED_TEST(LagrangeBasisFixture, PartitionOfUnity)
{
    using X = typename TestFixture::X;
    using GridX = typename TestFixture::GridX;
    using DataType = typename TestFixture::DataType;
    using LagBasis = typename TestFixture::LagBasis;
    using knot_discrete_dimension_type = std::conditional_t<
            TestFixture::UNIFORM,
            UniformLagrangeKnots<LagBasis>,
            NonUniformLagrangeKnots<LagBasis>>;

    constexpr std::size_t degree = LagBasis::degree();

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

    constexpr Idx<knot_discrete_dimension_type> poly_start(0);

    std::array<DataType, degree + 1> basis_storage;
    Span1D<DataType> basis_values(basis_storage.data(), degree + 1);

    // Test some points inside the cell
    for (int k = 1; k < 4; ++k) {
        Coord<X> x(
                DataType(ddc::coordinate(poly_start) + ddc::coordinate(poly_start + 1))
                * DataType(k) / DataType(8));

        ddc::discrete_space<LagBasis>().eval_basis(basis_values, x, poly_start);

        DataType sum = 0;
        for (std::size_t i = 0; i < degree + 1; ++i) {
            sum += basis_values[i];
        }

        EXPECT_NEAR(sum, DataType(1), TestFixture::TOL);
    }
}
