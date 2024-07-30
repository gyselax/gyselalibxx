// SPDX-License-Identifier: MIT

#include "9patches_2d_non_periodic_uniform.hpp"
#include <gtest/gtest.h>

TEST(MultipatchConnectivityTest, InterfaceConnections)
{
    EXPECT_TRUE(NorthInterface1::connected_to_patch<Patch1>());
    EXPECT_FALSE(NorthInterface1::connected_to_patch<Patch2>());
    EXPECT_FALSE(NorthInterface1::connected_to_patch<Patch3>());
    EXPECT_FALSE(NorthInterface1::connected_to_patch<Patch4>());
    EXPECT_FALSE(NorthInterface1::connected_to_patch<Patch5>());
    EXPECT_FALSE(NorthInterface1::connected_to_patch<Patch6>());
    EXPECT_FALSE(NorthInterface1::connected_to_patch<Patch7>());
    EXPECT_FALSE(NorthInterface1::connected_to_patch<Patch8>());
    EXPECT_FALSE(NorthInterface1::connected_to_patch<Patch9>());

    EXPECT_TRUE(Interface_1_2::connected_to_patch<Patch1>());
    EXPECT_TRUE(Interface_1_2::connected_to_patch<Patch2>());
    EXPECT_FALSE(Interface_1_2::connected_to_patch<Patch3>());
    EXPECT_FALSE(Interface_1_2::connected_to_patch<Patch4>());
    EXPECT_FALSE(Interface_1_2::connected_to_patch<Patch5>());
    EXPECT_FALSE(Interface_1_2::connected_to_patch<Patch6>());
    EXPECT_FALSE(Interface_1_2::connected_to_patch<Patch7>());
    EXPECT_FALSE(Interface_1_2::connected_to_patch<Patch8>());
    EXPECT_FALSE(Interface_1_2::connected_to_patch<Patch9>());

    EXPECT_TRUE(EastInterface3::connected_to_patch<Patch3>());
    EXPECT_FALSE(EastInterface3::connected_to_patch<Patch1>());
    EXPECT_FALSE(EastInterface3::connected_to_patch<Patch2>());
    EXPECT_FALSE(EastInterface3::connected_to_patch<Patch4>());
    EXPECT_FALSE(EastInterface3::connected_to_patch<Patch5>());
    EXPECT_FALSE(EastInterface3::connected_to_patch<Patch6>());
    EXPECT_FALSE(EastInterface3::connected_to_patch<Patch7>());
    EXPECT_FALSE(EastInterface3::connected_to_patch<Patch8>());
    EXPECT_FALSE(EastInterface3::connected_to_patch<Patch9>());
}

TEST(MultipatchConnectivityTest, GetConnections)
{
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch1>>, 4);
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch2>>, 4);
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch3>>, 4);
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch4>>, 4);
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch5>>, 4);
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch6>>, 4);
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch7>>, 4);
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch8>>, 4);
    EXPECT_EQ(std::tuple_size_v<Connectivity::template get_connections_t<Patch9>>, 4);

    using Patch1Interfaces = typename Connectivity::template get_type_seq_connections_t<Patch1>;
    bool north_interface_present = ddc::in_tags_v<NorthInterface1, Patch1Interfaces>;
    bool east_interface_present = ddc::in_tags_v<Interface_1_2, Patch1Interfaces>;
    bool west_interface_present = ddc::in_tags_v<Interface_1_4, Patch1Interfaces>;
    bool south_interface_present = ddc::in_tags_v<WestInterface1, Patch1Interfaces>;
    EXPECT_TRUE(north_interface_present);
    EXPECT_TRUE(east_interface_present);
    EXPECT_TRUE(west_interface_present);
    EXPECT_TRUE(south_interface_present);

    using Patch5Interfaces = typename Connectivity::template get_type_seq_connections_t<Patch5>;
    north_interface_present = ddc::in_tags_v<Interface_2_5, Patch5Interfaces>;
    east_interface_present = ddc::in_tags_v<Interface_5_6, Patch5Interfaces>;
    west_interface_present = ddc::in_tags_v<Interface_5_8, Patch5Interfaces>;
    south_interface_present = ddc::in_tags_v<Interface_4_5, Patch5Interfaces>;
    EXPECT_TRUE(north_interface_present);
    EXPECT_TRUE(east_interface_present);
    EXPECT_TRUE(west_interface_present);
    EXPECT_TRUE(south_interface_present);
}
