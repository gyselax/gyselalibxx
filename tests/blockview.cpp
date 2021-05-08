#include <gtest/gtest.h>

#include "blockview.h"

TEST(DBlockX, Constructor)
{
    MDomainX dom(0., 10., 0, 2);
    DBlockX block(dom);
}

TEST(DBlockX, get_domain)
{
    MDomainX dom(0., 10., 0, 2);
    DBlockX block(dom);
    const MDomainX& x_dom = block.domain<Dim::X>();
    ASSERT_EQ(dom, x_dom);
    const MDomainX& xb_dom = get_domain<Dim::X>(block);
    ASSERT_EQ(dom, xb_dom);
}

