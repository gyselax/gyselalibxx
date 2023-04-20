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
    ASSERT_EQ(map_v1, map_ref);
    ASSERT_EQ(map_v2, map_ref);
    ASSERT_EQ(map_v3, map_ref);
}

TEST(TaggedVector, Accessor)
{
    ddc::detail::TaggedVector<int, double, float> map(1, 2);

    ASSERT_EQ(map.get<double>(), 1);
    ASSERT_EQ(ddc::get<float>(map), 2);
}

TEST(TaggedVector, ConstAccessor)
{
    ddc::detail::TaggedVector<int, double, float> map(1, 2);
    ddc::detail::TaggedVector<int, double, float> const& cmap = map;

    ASSERT_EQ(cmap.get<double>(), 1);
    ASSERT_EQ(ddc::get<float>(cmap), 2);
}

TEST(TaggedVector, AccessorSingleElement)
{
    ddc::detail::TaggedVector<int, double> map(1);

    ASSERT_EQ(map.get<double>(), 1);
    ASSERT_EQ(ddc::get<double>(map), 1);
    ASSERT_EQ(map.value(), 1);
}

TEST(TaggedVector, Transpose)
{
    ddc::detail::TaggedVector<int, int, double, float> coord {0, 1, 2};
    ASSERT_EQ(coord.get<int>(), 0);
    ASSERT_EQ(coord.get<double>(), 1);
    ASSERT_EQ(coord.get<float>(), 2);

    ddc::detail::TaggedVector<int, double, float, int> coord_reordered(coord);
    ASSERT_EQ(coord.get<int>(), coord_reordered.get<int>());
    ASSERT_EQ(coord.get<double>(), coord_reordered.get<double>());
    ASSERT_EQ(coord.get<float>(), coord_reordered.get<float>());
}
