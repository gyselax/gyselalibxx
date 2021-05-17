#include <gtest/gtest.h>

#include "mcoord.h"

#include <experimental/mdspan>

using namespace std;
using namespace std::experimental;

TEST(MCoord, mcoord_end)
{
    extents<dynamic_extent, dynamic_extent> ex(10, 20);
    auto&& coord = mcoord_end<Dim::X, Dim::Vx>(ex);
    ASSERT_EQ(10, coord.template get<Dim::X>());
    ASSERT_EQ(20, coord.template get<Dim::Vx>());
}
