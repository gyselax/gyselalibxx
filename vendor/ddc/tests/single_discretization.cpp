// SPDX-License-Identifier: MIT
#include <memory>
#include <utility>

#include <ddc/experimental/single_discretization.hpp>

#include <gtest/gtest.h>

namespace {

using experimental::SingleDiscretization;

class DimX;

using CoordX = Coordinate<DimX>;

using DDimX = SingleDiscretization<DimX>;

using ElemX = DiscreteCoordinate<DDimX>;

} // namespace

TEST(SingleDiscretization, ClassSize)
{
    EXPECT_EQ(sizeof(DDimX), sizeof(double));
}

TEST(SingleDiscretization, Constructor)
{
    constexpr CoordX x(1.);

    SingleDiscretization<DimX> ddim_x(x);

    EXPECT_EQ(ddim_x.to_real(ElemX(0)), x);
}
