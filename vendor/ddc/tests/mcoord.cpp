#include <iosfwd>
#include <memory>

#include <gtest/gtest.h>

#include "gtest/gtest_pred_impl.h"

#include "mcoord.h"

#include <experimental/mdspan>

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
