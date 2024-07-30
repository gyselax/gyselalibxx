// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>
#include "9patches_2d_non_periodic_uniform.hpp"
#include "ddc_helper.hpp"

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

TEST(MultipatchConnectivityTest, PatchCollection)
{
    EXPECT_EQ(ddcHelper::type_seq_length_v<Connectivity::all_patches>, 9);
    bool patch_1_found = ddc::in_tags_v<Patch1, Connectivity::all_patches>;
    bool patch_2_found = ddc::in_tags_v<Patch2, Connectivity::all_patches>;
    bool patch_3_found = ddc::in_tags_v<Patch3, Connectivity::all_patches>;
    bool patch_4_found = ddc::in_tags_v<Patch4, Connectivity::all_patches>;
    bool patch_5_found = ddc::in_tags_v<Patch5, Connectivity::all_patches>;
    bool patch_6_found = ddc::in_tags_v<Patch6, Connectivity::all_patches>;
    bool patch_7_found = ddc::in_tags_v<Patch7, Connectivity::all_patches>;
    bool patch_8_found = ddc::in_tags_v<Patch8, Connectivity::all_patches>;
    bool patch_9_found = ddc::in_tags_v<Patch9, Connectivity::all_patches>;
    EXPECT_TRUE(patch_1_found);
    EXPECT_TRUE(patch_2_found);
    EXPECT_TRUE(patch_3_found);
    EXPECT_TRUE(patch_4_found);
    EXPECT_TRUE(patch_5_found);
    EXPECT_TRUE(patch_6_found);
    EXPECT_TRUE(patch_7_found);
    EXPECT_TRUE(patch_8_found);
    EXPECT_TRUE(patch_9_found);
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
