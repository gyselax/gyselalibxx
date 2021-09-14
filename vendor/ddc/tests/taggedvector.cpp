#include <memory>

#include <ddc/taggedvector.h>

#include <gtest/gtest.h>

TEST(TaggedVector, Constructor)
{
    TaggedVector<int, double, float> map(1, 2);
}

TEST(TaggedVector, Accessor)
{
    TaggedVector<int, double, float> map(1, 2);

    ASSERT_EQ(map.get<double>(), 1);
    ASSERT_EQ(get<float>(map), 2);
}

TEST(TaggedVector, ConstAccessor)
{
    TaggedVector<int, double, float> map(1, 2);
    const auto& cmap = map;

    ASSERT_EQ(cmap.get<double>(), 1);
    ASSERT_EQ(get<float>(cmap), 2);
}

TEST(TaggedVector, Transpose)
{
    TaggedVector<int, int, double, float> coord {0, 1, 2};
    ASSERT_EQ(coord.get<int>(), 0);
    ASSERT_EQ(coord.get<double>(), 1);
    ASSERT_EQ(coord.get<float>(), 2);

    TaggedVector<int, double, float, int> coord_reordered(coord);
    ASSERT_EQ(coord.get<int>(), coord_reordered.get<int>());
    ASSERT_EQ(coord.get<double>(), coord_reordered.get<double>());
    ASSERT_EQ(coord.get<float>(), coord_reordered.get<float>());
}
