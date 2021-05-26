#include <gtest/gtest.h>

#include "taggedarray.h"

TEST(TagRank, Rank)
{
    using namespace detail;
    ASSERT_EQ((RankIn<SingleType<int>, TypeSeq<int, double, float>>::val), 0);
    ASSERT_EQ((RankIn<SingleType<double>, TypeSeq<int, double, float>>::val), 1);
    ASSERT_EQ((RankIn<SingleType<float>, TypeSeq<int, double, float>>::val), 2);
}

TEST(TaggedArray, Constructor)
{
    TaggedArray<int, double, float> map(1, 2);
}

TEST(TaggedArray, Accessor)
{
    TaggedArray<int, double, float> map(1, 2);

    ASSERT_EQ(map.get<double>(), 1);
    ASSERT_EQ(get<float>(map), 2);
}

TEST(TaggedArray, ConstAccessor)
{
    TaggedArray<int, double, float> map(1, 2);
    const auto& cmap = map;

    ASSERT_EQ(cmap.get<double>(), 1);
    ASSERT_EQ(get<float>(cmap), 2);
}

TEST(TaggedArray, Transpose)
{
    TaggedArray<int, int, double, float> coord {0, 1, 2};
    ASSERT_EQ(coord.get<int>(), 0);
    ASSERT_EQ(coord.get<double>(), 1);
    ASSERT_EQ(coord.get<float>(), 2);

    TaggedArray<int, double, float, int> coord_reordered(coord);
    ASSERT_EQ(coord.get<int>(), coord_reordered.get<int>());
    ASSERT_EQ(coord.get<double>(), coord_reordered.get<double>());
    ASSERT_EQ(coord.get<float>(), coord_reordered.get<float>());
}
