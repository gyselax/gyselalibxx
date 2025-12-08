#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

int main(int argc, char** argv)
{
    ::testing::InitGoogleMock(&argc, argv);
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);
    return RUN_ALL_TESTS();
}
