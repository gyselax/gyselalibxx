#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace testing::internal {
// accessing gtest internals is not very clean, but gtest provides no public access...
extern bool g_help_flag;
} // namespace testing::internal

int main(int argc, char** argv)
{
    ::testing::InitGoogleMock(&argc, argv);
    if (::testing::GTEST_FLAG(list_tests) || ::testing::internal::g_help_flag) {
        // do not initialize DDC just to list tests or show help so as to be able to run
        // a Cuda/Hip/etc. enabled DDC with no device available
        return RUN_ALL_TESTS();
    } else {
        ::ddc::ScopeGuard scope(argc, argv);
        return RUN_ALL_TESTS();
    }
}
