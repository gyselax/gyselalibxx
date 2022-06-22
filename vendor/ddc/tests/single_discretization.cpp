// SPDX-License-Identifier: MIT
#include <memory>
#include <utility>

#include <ddc/experimental/single_discretization.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace {

using experimental::SingleDiscretization;

class DimX;

using CoordX = Coordinate<DimX>;

using DDimX = SingleDiscretization<DimX>;

using DElemX = DiscreteElement<DDimX>;

} // namespace

TEST(SingleDiscretization, ClassSize)
{
    EXPECT_EQ(sizeof(DDimX::Impl<Kokkos::HostSpace>), sizeof(double));
}

TEST(SingleDiscretization, Constructor)
{
    constexpr CoordX x(1.);

    SingleDiscretization<DimX>::Impl<Kokkos::HostSpace> ddim_x(x);

    EXPECT_EQ(ddim_x.coordinate(DElemX(0)), x);
}
