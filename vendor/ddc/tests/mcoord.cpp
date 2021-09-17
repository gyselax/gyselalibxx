// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>

#include <experimental/mdspan>

#include <ddc/MCoord>
#include <ddc/TaggedVector>

#include <gtest/gtest.h>

using namespace std;
using namespace std::experimental;

class DimX;
class DimVx;

TEST(MCoord, mcoord_end)
{
    extents<dynamic_extent, dynamic_extent> ex(10, 20);
    auto&& coord = mcoord_end<DimX, DimVx>(ex);
    ASSERT_EQ(10, coord.template get<DimX>());
    ASSERT_EQ(20, coord.template get<DimVx>());
}
