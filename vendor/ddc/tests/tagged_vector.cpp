// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <utility>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

TEST(TaggedVector, Constructor)
{
    [[maybe_unused]] ddc::detail::TaggedVector<int, double, float> const map(1, 2);
}

TEST(TaggedVector, ConstructorFromTaggedVectors)
{
    ddc::detail::TaggedVector<int, float, double> const map_v1(1, 2);
    ddc::detail::TaggedVector<int, long double> const map_v2(3);
    ddc::detail::TaggedVector<int, float, long double, double> const map_v3(map_v1, map_v2);
    EXPECT_EQ(map_v3.get<float>(), 1);
    EXPECT_EQ(map_v3.get<double>(), 2);
    EXPECT_EQ(map_v3.get<long double>(), 3);
}

TEST(TaggedVector, ReorderingConstructor)
{
    ddc::detail::TaggedVector<int, double, float> const map_ref(1, 2);
    ddc::detail::TaggedVector<int, double> const submap_double = ddc::select<double>(map_ref);
    ddc::detail::TaggedVector<int, float> const submap_float = ddc::select<float>(map_ref);
    ddc::detail::TaggedVector<int, double, float> const map_v1(submap_double, submap_float);
    ddc::detail::TaggedVector<int, double, float> const map_v2(submap_float, submap_double);
    ddc::detail::TaggedVector<int, float, double> const map_v3(map_ref);
    EXPECT_EQ(map_v1, map_ref);
    EXPECT_EQ(map_v2, map_ref);
    EXPECT_EQ(map_v3, map_ref);
}

TEST(TaggedVector, Accessor)
{
    ddc::detail::TaggedVector<int, double, float> const map(1, 2);

    EXPECT_EQ(map.get<double>(), 1);
    EXPECT_EQ(ddc::get<float>(map), 2);
}

TEST(TaggedVector, AccessorSingleElement)
{
    ddc::detail::TaggedVector<int, double> const map(1);

    EXPECT_EQ(map.get<double>(), 1);
    EXPECT_EQ(ddc::get<double>(map), 1);
    EXPECT_EQ(map.value(), 1);
}

TEST(TaggedVector, Transpose)
{
    ddc::detail::TaggedVector<int, int, double, float> const coord(0, 1, 2);
    EXPECT_EQ(coord.get<int>(), 0);
    EXPECT_EQ(coord.get<double>(), 1);
    EXPECT_EQ(coord.get<float>(), 2);

    ddc::detail::TaggedVector<int, double, float, int> const coord_reordered(coord);
    EXPECT_EQ(coord.get<int>(), coord_reordered.get<int>());
    EXPECT_EQ(coord.get<double>(), coord_reordered.get<double>());
    EXPECT_EQ(coord.get<float>(), coord_reordered.get<float>());
}

TEST(TaggedVector, Operators)
{
    ddc::detail::TaggedVector<int, double, float> const a(1, 2);
    ddc::detail::TaggedVector<int, float, double> const b(3, 4);
    ddc::detail::TaggedVector<int, double> const c = ddc::select<double>(a);
    EXPECT_EQ(a + b, (ddc::detail::TaggedVector<int, double, float>(5, 5)));
    EXPECT_EQ(b - a, (ddc::detail::TaggedVector<int, double, float>(3, 1)));
    EXPECT_EQ(c + 4, (ddc::detail::TaggedVector<int, double>(5)));
    EXPECT_EQ(4 + c, (ddc::detail::TaggedVector<int, double>(5)));
    EXPECT_EQ(4 * a, (ddc::detail::TaggedVector<int, double, float>(4, 8)));
}

TEST(TaggedVector, Assignment)
{
    ddc::detail::TaggedVector<int, double, float> const a(1, 2);
    ddc::detail::TaggedVector<int, double, float> b;
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

TEST(TaggedVector, Conversion)
{
    ddc::detail::TaggedVector<int, float, double> const a(1, 2);
    ddc::detail::TaggedVector<double, float, double> const b(a);
    EXPECT_EQ(b.get<float>(), 1.0);
    EXPECT_EQ(b.get<double>(), 2.0);
}

TEST(TaggedVector, ConversionReorder)
{
    ddc::detail::TaggedVector<int, float, double> const a(1, 2);
    ddc::detail::TaggedVector<double, double, float> const b(a);
    EXPECT_EQ(b.get<float>(), 1.0);
    EXPECT_EQ(b.get<double>(), 2.0);
}
