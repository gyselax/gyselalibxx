#include <gtest/gtest.h>

#include "blockview.h"

MDomainX dom(0., 10., 0, 2);

TEST(DBlockX, Constructor)
{
    DBlockX block(dom);
}

TEST(DBlockX, domain)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, block.domain());
}

TEST(DBlockX, domainX)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, block.domain<Dim::X>());
}

TEST(DBlockX, get_domainX)
{
    MDomainX dom(0., 10., 0, 2);
    DBlockX block(dom);
    ASSERT_EQ(dom, get_domain<Dim::X>(block));
}

TEST(DBlockX, deepcopy)
{
    constexpr auto NB_ITER = 10;
    DBlockX block(MDomainX(0., 10., 0, NB_ITER));
    for (auto&& ii : block.domain()) {
        block(ii) = 1.001 * ii;
    }
    DBlockX block2(block.domain());
    deepcopy(block2, block);
    for (auto&& ii : block.domain()) {
        // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
        ASSERT_EQ(block2(ii), block(ii));
    }
}
