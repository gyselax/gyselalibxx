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

template <class Patch>
typename Patch::IdxRange12 get_idx_range(int start1, int size1, int start2, int size2)
{
    using Grid1 = typename Patch::Grid1;
    using Grid2 = typename Patch::Grid2;
    Idx<Grid1> front1(start1);
    IdxStep<Grid1> length1(size1);
    Idx<Grid2> front2(start2);
    IdxStep<Grid2> length2(size2);
    IdxRange<Grid1> idx_range_1(front1, length1);
    IdxRange<Grid2> idx_range_2(front2, length2);
    return IdxRange<Grid1, Grid2>(idx_range_1, idx_range_2);
}

TEST(MultipatchConnectivityTest, GetAllIndexRangesAlongDim)
{
    Patch1::IdxRange12 idx_range_1 = get_idx_range<Patch1>(0, 5, 9, 10);
    Patch2::IdxRange12 idx_range_2 = get_idx_range<Patch2>(1, 5, 8, 10);
    Patch3::IdxRange12 idx_range_3 = get_idx_range<Patch3>(2, 5, 7, 10);
    Patch4::IdxRange12 idx_range_4 = get_idx_range<Patch4>(3, 5, 6, 10);
    Patch5::IdxRange12 idx_range_5 = get_idx_range<Patch5>(4, 5, 5, 10);
    Patch6::IdxRange12 idx_range_6 = get_idx_range<Patch6>(5, 5, 4, 10);
    Patch7::IdxRange12 idx_range_7 = get_idx_range<Patch7>(6, 5, 3, 10);
    Patch8::IdxRange12 idx_range_8 = get_idx_range<Patch8>(7, 5, 2, 10);
    Patch9::IdxRange12 idx_range_9 = get_idx_range<Patch9>(8, 5, 1, 10);
    std::tuple all_domains
            = {idx_range_1,
               idx_range_2,
               idx_range_3,
               idx_range_4,
               idx_range_5,
               idx_range_6,
               idx_range_7,
               idx_range_8,
               idx_range_9};

    using StartGrid = typename Patch2::Grid2;

    auto grids = Connectivity::template get_all_idx_ranges_along_direction<StartGrid>(all_domains);
    (void)(grids);
}
