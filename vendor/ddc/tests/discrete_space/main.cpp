// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "discrete_space.hpp"

namespace testing::internal {
// accessing gtest internals is not very clean, but gtest provides no public access...
extern bool g_help_flag;
} // namespace testing::internal

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    do_not_optimize_away_discrete_space_tests();
    if (::testing::GTEST_FLAG(list_tests) || ::testing::internal::g_help_flag) {
        // do not initialize DDC just to list tests or show help so as to be able to run
        // a Cuda/Hip/etc. enabled DDC with no device available
        return RUN_ALL_TESTS();
    } else {
        Kokkos::ScopeGuard const kokkos_scope(argc, argv);
        ddc::ScopeGuard const ddc_scope(argc, argv);
        return RUN_ALL_TESTS();
    }
}
