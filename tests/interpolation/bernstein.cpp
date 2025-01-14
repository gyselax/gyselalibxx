#include <random>

#include <ddc/ddc.hpp>

#include <sll/mapping/mapping_tools.hpp>
#include <sll/test_utils.hpp>

#include <gtest/gtest.h>

#include "bernstein.hpp"
#include "ddc_alias_inline_functions.hpp"

template <class Tag1, class Tag2>
Coord<Tag1, Tag2> generate_random_point_in_triangle(
        Coord<Tag1, Tag2> const& corner1,
        Coord<Tag1, Tag2> const& corner2,
        Coord<Tag1, Tag2> const& corner3)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double rand1 = dist(mt);
    double rand2 = dist(mt);
    if (rand1 + rand2 > 1) {
        rand1 = 1 - rand1;
        rand2 = 1 - rand2;
    }
    const double c1_x = ddc::get<Tag1>(corner1);
    const double c1_y = ddc::get<Tag2>(corner1);
    const double c2_x = ddc::get<Tag1>(corner2);
    const double c2_y = ddc::get<Tag2>(corner2);
    const double c3_x = ddc::get<Tag1>(corner3);
    const double c3_y = ddc::get<Tag2>(corner3);
    const double point_x = c1_x + (c2_x - c1_x) * rand1 + (c3_x - c1_x) * rand2;
    const double point_y = c1_y + (c2_y - c1_y) * rand1 + (c3_y - c1_y) * rand2;

    return Coord<Tag1, Tag2>(point_x, point_y);
}

template <class T>
struct BernsteinFixture;

template <std::size_t D>
struct BernsteinFixture<std::tuple<std::integral_constant<std::size_t, D>>> : public testing::Test
{
    struct X
    {
        static constexpr bool PERIODIC = false;
    };
    struct Y
    {
        static constexpr bool PERIODIC = false;
    };
    struct Corner1
    {
    };
    struct Corner2
    {
    };
    struct Corner3
    {
    };
    static constexpr std::size_t poly_degree = D;
    struct Bernstein
        : TriangularBernsteinPolynomialBasis<X, Y, Corner1, Corner2, Corner3, poly_degree>
    {
    };
};

using degrees = std::integer_sequence<std::size_t, 0, 1, 2, 3>;

using Cases = tuple_to_types_t<cartesian_product_t<degrees>>;

TYPED_TEST_SUITE(BernsteinFixture, Cases);

TYPED_TEST(BernsteinFixture, PartitionOfUnity)
{
    using X = typename TestFixture::X;
    using Y = typename TestFixture::Y;
    using Corner1 = typename TestFixture::Corner1;
    using Corner2 = typename TestFixture::Corner2;
    using Corner3 = typename TestFixture::Corner3;
    using Bernstein = typename TestFixture::Bernstein;
    using CoordXY = Coord<X, Y>;

    const CoordXY c1(-1.0, -1.0);
    const CoordXY c2(0.0, 1.0);
    const CoordXY c3(1.0, -1.0);

    CartesianToBarycentric<X, Y, Corner1, Corner2, Corner3> coordinate_converter(c1, c2, c3);
    static_assert(is_mapping_v<CartesianToBarycentric<X, Y, Corner1, Corner2, Corner3>>);
    ddc::init_discrete_space<Bernstein>(coordinate_converter);

    IdxRange<Bernstein> idx_range(Idx<Bernstein>(0), IdxStep<Bernstein>(Bernstein::nbasis()));

    host_t<DFieldMem<IdxRange<Bernstein>>> values(idx_range);

    std::size_t const n_test_points = 100;
    for (std::size_t i(0); i < n_test_points; ++i) {
        CoordXY const test_point = generate_random_point_in_triangle(c1, c2, c3);
        ddc::discrete_space<Bernstein>().eval_basis(get_field(values), test_point);
        double total = ddc::transform_reduce(
                idx_range,
                0.0,
                ddc::reducer::sum<double>(),
                [&](Idx<Bernstein> const ix) { return values(ix); });
        EXPECT_LE(fabs(total - 1.0), 1.0e-15);
    }
}
