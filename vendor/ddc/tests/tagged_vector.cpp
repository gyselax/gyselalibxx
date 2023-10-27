// SPDX-License-Identifier: MIT
#include <memory>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

TEST(TaggedVector, Constructor)
{
    [[maybe_unused]] ddc::detail::TaggedVector<int, double, float> map(1, 2);
}

TEST(TaggedVector, ReorderingConstructor)
{
    ddc::detail::TaggedVector<int, double, float> map_ref(1, 2);
    ddc::detail::TaggedVector<int, double> submap_double = ddc::select<double>(map_ref);
    ddc::detail::TaggedVector<int, float> submap_float = ddc::select<float>(map_ref);
    ddc::detail::TaggedVector<int, double, float> map_v1(submap_double, submap_float);
    ddc::detail::TaggedVector<int, double, float> map_v2(submap_float, submap_double);
    ddc::detail::TaggedVector<int, float, double> map_v3(map_ref);
    EXPECT_EQ(map_v1, map_ref);
    EXPECT_EQ(map_v2, map_ref);
    EXPECT_EQ(map_v3, map_ref);
}

TEST(TaggedVector, Accessor)
{
    ddc::detail::TaggedVector<int, double, float> map(1, 2);

    EXPECT_EQ(map.get<double>(), 1);
    EXPECT_EQ(ddc::get<float>(map), 2);
}

TEST(TaggedVector, ConstAccessor)
{
    ddc::detail::TaggedVector<int, double, float> map(1, 2);
    ddc::detail::TaggedVector<int, double, float> const& cmap = map;

    EXPECT_EQ(cmap.get<double>(), 1);
    EXPECT_EQ(ddc::get<float>(cmap), 2);
}

TEST(TaggedVector, AccessorSingleElement)
{
    ddc::detail::TaggedVector<int, double> map(1);

    EXPECT_EQ(map.get<double>(), 1);
    EXPECT_EQ(ddc::get<double>(map), 1);
    EXPECT_EQ(map.value(), 1);
}

TEST(TaggedVector, Transpose)
{
    ddc::detail::TaggedVector<int, int, double, float> coord(0, 1, 2);
    EXPECT_EQ(coord.get<int>(), 0);
    EXPECT_EQ(coord.get<double>(), 1);
    EXPECT_EQ(coord.get<float>(), 2);

    ddc::detail::TaggedVector<int, double, float, int> coord_reordered(coord);
    EXPECT_EQ(coord.get<int>(), coord_reordered.get<int>());
    EXPECT_EQ(coord.get<double>(), coord_reordered.get<double>());
    EXPECT_EQ(coord.get<float>(), coord_reordered.get<float>());
}

TEST(TaggedVector, Operators)
{
    ddc::detail::TaggedVector<int, double, float> a(1, 2);
    ddc::detail::TaggedVector<int, float, double> b(3, 4);
    ddc::detail::TaggedVector<int, double> c = ddc::select<double>(a);
    ASSERT_EQ(a + b, (ddc::detail::TaggedVector<int, double, float>(5, 5)));
    ASSERT_EQ(b - a, (ddc::detail::TaggedVector<int, double, float>(3, 1)));
    ASSERT_EQ(c + 4, (ddc::detail::TaggedVector<int, double>(5)));
    ASSERT_EQ(4 + c, (ddc::detail::TaggedVector<int, double>(5)));
    ASSERT_EQ(4 * a, (ddc::detail::TaggedVector<int, double, float>(4, 8)));
}

TEST(TaggedVector, Assignment)
{
    ddc::detail::TaggedVector<int, double, float> a(1, 2);
    ddc::detail::TaggedVector<int, double, float> b;
    b = a;
    EXPECT_EQ(a.get<double>(), b.get<double>());
    EXPECT_EQ(a.get<float>(), b.get<float>());
}

TEST(TaggedVector, ReorderingAssignment)
{
    ddc::detail::TaggedVector<int, double, float> a(1, 2);
    ddc::detail::TaggedVector<int, float, double> b;
    b = a;
    EXPECT_EQ(a.get<double>(), b.get<double>());
    EXPECT_EQ(a.get<float>(), b.get<float>());
}

TEST(TaggedVector, MoveAssignment)
{
    ddc::detail::TaggedVector<int, double, float> a(1, 2);
    ddc::detail::TaggedVector<int, double, float> b;
    b = std::move(a);
    EXPECT_EQ(1, b.get<double>());
    EXPECT_EQ(2, b.get<float>());
}

TEST(TaggedVector, ReorderingMoveAssignment)
{
    ddc::detail::TaggedVector<int, double, float> a(1, 2);
    ddc::detail::TaggedVector<int, float, double> b;
    b = std::move(a);
    EXPECT_EQ(1, b.get<double>());
    EXPECT_EQ(2, b.get<float>());
}
