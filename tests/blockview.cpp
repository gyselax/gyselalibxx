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

