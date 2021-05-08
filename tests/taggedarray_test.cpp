#include <gtest/gtest.h>

#include "taggedarray.h"

TEST(TaggedArray, Constructor)
{
    TaggedArray<int, double, float> map {1, 2};
}

TEST(TaggedArray, Accessor)
{
    TaggedArray<int, double, float> map {1, 2};
    
    ASSERT_EQ(map.get<double>(), 1);
    ASSERT_EQ(get<float>(map), 2);
}

TEST(TaggedArray, ConstAccessor)
{
    TaggedArray<int, double, float> map {1, 2};
    const auto& cmap = map;
    
    ASSERT_EQ(cmap.get<double>(), 1);
    ASSERT_EQ(get<float>(cmap), 2);
}

