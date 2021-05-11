#include <gtest/gtest.h>

#include "blockview.h"

class DBlockXTest : public ::testing::Test
{
protected:
    MDomainX dom = MDomainX(0., 10., 0, 2);
};



TEST_F(DBlockXTest, constructor)
{
    DBlockX block(dom);
}

TEST_F(DBlockXTest, domain)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, block.domain());
}

TEST_F(DBlockXTest, domainX)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, block.domain<Dim::X>());
}

TEST_F(DBlockXTest, get_domainX)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, get_domain<Dim::X>(block));
}

TEST_F(DBlockXTest, access)
{
    DBlockX block(dom);
    for (auto&& ii : block.domain()) {
        ASSERT_EQ(block(ii), block(std::array {ii}));
    }
}

TEST_F(DBlockXTest, deepcopy)
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
