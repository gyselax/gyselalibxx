// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>

#include "2patches_2d_non_periodic_uniform.hpp"
#include "5patches_figure_of_eight.hpp"
#include "9patches_2d_non_periodic_uniform.hpp"
#include "ddc_helper.hpp"

TEST(MultipatchConnectivityTest, InterfaceConnections)
{
    using namespace non_periodic_uniform_2d_9patches;
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

    EXPECT_TRUE(SouthInterface7::connected_to_patch<Patch7>());
    EXPECT_FALSE(SouthInterface7::connected_to_patch<Patch1>());
    EXPECT_FALSE(SouthInterface7::connected_to_patch<Patch2>());
    EXPECT_FALSE(SouthInterface7::connected_to_patch<Patch3>());
    EXPECT_FALSE(SouthInterface7::connected_to_patch<Patch4>());
    EXPECT_FALSE(SouthInterface7::connected_to_patch<Patch5>());
    EXPECT_FALSE(SouthInterface7::connected_to_patch<Patch6>());
    EXPECT_FALSE(SouthInterface7::connected_to_patch<Patch8>());
    EXPECT_FALSE(SouthInterface7::connected_to_patch<Patch9>());
}

TEST(MultipatchConnectivityTest, PatchCollection)
{
    using namespace non_periodic_uniform_2d_9patches;
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
    using namespace non_periodic_uniform_2d_9patches;
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
    bool south_interface_present = ddc::in_tags_v<Interface_3_1, Patch1Interfaces>;
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
    using namespace non_periodic_uniform_2d_9patches;
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

    using StartGridNonPeriodic = typename Patch2::Grid2;

    auto grids_y = Connectivity::template get_all_idx_ranges_along_direction<StartGridNonPeriodic>(
            all_domains);
    static_assert(std::tuple_size_v<decltype(grids_y)> == 3);
    EXPECT_EQ(std::get<0>(grids_y), ddc::select<GridY<8>>(idx_range_8));
    EXPECT_EQ(std::get<1>(grids_y), ddc::select<GridY<5>>(idx_range_5));
    EXPECT_EQ(std::get<2>(grids_y), ddc::select<GridY<2>>(idx_range_2));

    using StartGridPeriodic = typename Patch2::Grid1;

    auto grids_x = Connectivity::template get_all_idx_ranges_along_direction<StartGridPeriodic>(
            all_domains);
    static_assert(std::tuple_size_v<decltype(grids_x)> == 3);
    EXPECT_EQ(std::get<IdxRange<GridX<1>>>(grids_x), ddc::select<GridX<1>>(idx_range_1));
    EXPECT_EQ(std::get<IdxRange<GridX<2>>>(grids_x), ddc::select<GridX<2>>(idx_range_2));
    EXPECT_EQ(std::get<IdxRange<GridX<3>>>(grids_x), ddc::select<GridX<3>>(idx_range_3));
}

TEST(MultipatchConnectivityTest, GetAllIndexRangesAlongDimSimple)
{
    using namespace non_periodic_uniform_2d_2patches;
    Patch1::IdxRange12 idx_range_1 = get_idx_range<Patch1>(0, 5, 9, 10);
    Patch2::IdxRange12 idx_range_2 = get_idx_range<Patch2>(1, 5, 8, 10);
    std::tuple all_domains = {idx_range_1, idx_range_2};

    using StartGridY = typename Patch2::Grid2;

    auto grids_y
            = Connectivity::template get_all_idx_ranges_along_direction<StartGridY>(all_domains);
    static_assert(std::tuple_size_v<decltype(grids_y)> == 1);
    EXPECT_EQ(std::get<0>(grids_y), ddc::select<GridY2>(idx_range_2));

    using StartGridX = typename Patch2::Grid1;

    auto grids_x
            = Connectivity::template get_all_idx_ranges_along_direction<StartGridX>(all_domains);
    static_assert(std::tuple_size_v<decltype(grids_x)> == 2);
    static_assert(
            std::is_same_v<typename Patch1::IdxRange1, std::tuple_element_t<0, decltype(grids_x)>>);
    static_assert(
            std::is_same_v<typename Patch2::IdxRange1, std::tuple_element_t<1, decltype(grids_x)>>);
    EXPECT_EQ(std::get<0>(grids_x), ddc::select<GridX1>(idx_range_1));
    EXPECT_EQ(std::get<1>(grids_x), ddc::select<GridX2>(idx_range_2));
}

TEST(MultipatchConnectivityTest, GetAllIndexRangesAlongFigureOfEightDim)
{
    using namespace figure_of_eight_5patches;
    Patch1::IdxRange12 idx_range_1 = get_idx_range<Patch1>(0, 5, 9, 10);
    Patch2::IdxRange12 idx_range_2 = get_idx_range<Patch2>(1, 5, 8, 10);
    Patch3::IdxRange12 idx_range_3 = get_idx_range<Patch3>(2, 5, 7, 10);
    Patch4::IdxRange12 idx_range_4 = get_idx_range<Patch4>(3, 5, 6, 10);
    Patch5::IdxRange12 idx_range_5 = get_idx_range<Patch5>(4, 5, 5, 10);
    std::tuple all_domains = {idx_range_1, idx_range_2, idx_range_3, idx_range_4, idx_range_5};

    using MiddleGridX = typename Patch3::Grid1;

    auto figure_eight_grid
            = Connectivity::template get_all_idx_ranges_along_direction<MiddleGridX>(all_domains);
    static_assert(std::tuple_size_v<decltype(figure_eight_grid)> == 6);
    EXPECT_EQ(std::get<IdxRange<GridX<3>>>(figure_eight_grid), ddc::select<GridX<3>>(idx_range_3));
    EXPECT_EQ(std::get<IdxRange<GridX<4>>>(figure_eight_grid), ddc::select<GridX<4>>(idx_range_4));
    EXPECT_EQ(std::get<IdxRange<GridY<5>>>(figure_eight_grid), ddc::select<GridY<5>>(idx_range_5));
    EXPECT_EQ(std::get<IdxRange<GridY<3>>>(figure_eight_grid), ddc::select<GridY<3>>(idx_range_3));
    EXPECT_EQ(std::get<IdxRange<GridY<1>>>(figure_eight_grid), ddc::select<GridY<1>>(idx_range_1));
    EXPECT_EQ(std::get<IdxRange<GridX<2>>>(figure_eight_grid), ddc::select<GridX<2>>(idx_range_2));
}

TEST(MultipatchConnectivityDetailsTest, FilterEdges)
{
    using namespace non_periodic_uniform_2d_2patches;
    using FilteredEdges = Connectivity::inner_edges;
    static_assert(ddc::in_tags_v<NorthEdge1, FilteredEdges>);
    static_assert(ddc::in_tags_v<SouthEdge1, FilteredEdges>);
    static_assert(ddc::in_tags_v<EastEdge1, FilteredEdges>);
    static_assert(ddc::in_tags_v<WestEdge1, FilteredEdges>);
    static_assert(ddc::in_tags_v<NorthEdge2, FilteredEdges>);
    static_assert(ddc::in_tags_v<SouthEdge2, FilteredEdges>);
    static_assert(ddc::in_tags_v<EastEdge2, FilteredEdges>);
    static_assert(ddc::in_tags_v<WestEdge2, FilteredEdges>);
    static_assert(ddcHelper::type_seq_length_v<FilteredEdges> == 8);
}

TEST(MultipatchConnectivityDetailsTest, ExtractPatches)
{
    using namespace non_periodic_uniform_2d_2patches;
    using ExtractedPatches = Connectivity::all_patches;
    static_assert(ddcHelper::type_seq_length_v<ExtractedPatches> == 2);
    static_assert(ddc::in_tags_v<Patch1, ExtractedPatches>);
    static_assert(ddc::in_tags_v<Patch2, ExtractedPatches>);
}

TEST(MultipatchConnectivityDetailsTest, FindPatch)
{
    using namespace non_periodic_uniform_2d_9patches;
    using GridX1Patch = find_patch_t<GridX<1>, Connectivity::all_patches>;
    static_assert(std::is_same_v<GridX1Patch, Patch1>);
    using GridY7Patch = find_patch_t<GridY<7>, Connectivity::all_patches>;
    static_assert(std::is_same_v<GridY7Patch, Patch7>);
}

TEST(MultipatchConnectivityDetailsTest, FindAssociatedInterface)
{
    using namespace non_periodic_uniform_2d_2patches;
    using NorthEdge1Interface
            = find_associated_interface_t<NorthEdge1, Connectivity::interface_collection>;
    static_assert(std::is_same_v<NorthEdge1Interface, NorthInterface1>);
    using SouthEdge1Interface
            = find_associated_interface_t<SouthEdge1, Connectivity::interface_collection>;
    static_assert(std::is_same_v<SouthEdge1Interface, SouthInterface1>);
    using WestEdge2Interface
            = find_associated_interface_t<WestEdge2, Connectivity::interface_collection>;
    static_assert(std::is_same_v<WestEdge2Interface, Interface_1_2>);
}
