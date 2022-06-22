// SPDX-License-Identifier: MIT

#include <type_traits>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace {

struct X;
struct Y;

using DDimX = UniformPointSampling<X>;
using DDimY = UniformPointSampling<Y>;
using NUDDimX = NonUniformPointSampling<X>;
using NUDDimY = NonUniformPointSampling<Y>;

using DDomNull = DiscreteDomain<std::nullptr_t>;
using DDomX = DiscreteDomain<DDimX>;
using NUDDomX = DiscreteDomain<NUDDimX>;
using DDomXY = DiscreteDomain<DDimX, DDimY>;
using NUDDomXY = DiscreteDomain<NUDDimX, NUDDimY>;
using DDomXNUDDomY = DiscreteDomain<DDimX, NUDDimY>;

} // namespace

TEST(RectilinearDomainTest, NonDiscreteDomainSpecializationValue)
{
    EXPECT_FALSE(is_rectilinear_domain_v<std::nullptr_t>);
    EXPECT_FALSE(is_rectilinear_domain_v<X>);
}

TEST(RectilinearDomainTest, DiscreteDomainSpecializationValue)
{
    EXPECT_FALSE(is_rectilinear_domain_v<DDomNull>);

    EXPECT_TRUE(is_rectilinear_domain_v<DDomX>);
    EXPECT_TRUE(is_rectilinear_domain_v<DDomXY>);
    EXPECT_TRUE(is_rectilinear_domain_v<NUDDomX>);
    EXPECT_TRUE(is_rectilinear_domain_v<NUDDomXY>);
    EXPECT_TRUE(is_rectilinear_domain_v<DDomXNUDDomY>);
}
