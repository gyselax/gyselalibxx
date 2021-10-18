// SPDX-License-Identifier: MIT
#include <memory>
#include <utility>

#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/NonUniformDiscretization>
#include <ddc/experimental/single_discretization.hpp>

#include <gtest/gtest.h>

namespace {

using experimental::SingleDiscretization;

class DimX;

using CoordX = Coordinate<DimX>;

using DDimX = SingleDiscretization<DimX>;

using IndexX = DiscreteCoordinate<DDimX>;

} // namespace

TEST(SingleDiscretization, class_size)
{
    EXPECT_EQ(sizeof(DDimX), sizeof(double));
}

TEST(SingleDiscretization, constructor)
{
    constexpr CoordX x(1.);

    SingleDiscretization<DimX> ddim_x(x);

    EXPECT_EQ(ddim_x.to_real(IndexX(0)), x);
}
