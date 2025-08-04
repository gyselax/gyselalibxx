// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_non_periodic_uniform.hpp"
#include "3patches_2d_non_periodic_non_uniform.hpp"
#include "5patches_figure_of_eight.hpp"
#include "9patches_2d_periodic_strips_uniform.hpp"
#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"
#include "types.hpp"

namespace {

template <class Patch>
typename Patch::IdxRange12 build_idx_range(int start1, int size1, int start2, int size2)
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

} // namespace

TEST(MultipatchConnectivityTest, InterfaceConnections)
{
    using namespace periodic_strips_uniform_2d_9patches;
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
    using namespace periodic_strips_uniform_2d_9patches;
    EXPECT_EQ(ddc::type_seq_size_v<Connectivity::all_patches>, 9);
    EXPECT_TRUE((ddc::in_tags_v<Patch1, Connectivity::all_patches>));
    EXPECT_TRUE((ddc::in_tags_v<Patch2, Connectivity::all_patches>));
    EXPECT_TRUE((ddc::in_tags_v<Patch3, Connectivity::all_patches>));
    EXPECT_TRUE((ddc::in_tags_v<Patch4, Connectivity::all_patches>));
    EXPECT_TRUE((ddc::in_tags_v<Patch5, Connectivity::all_patches>));
    EXPECT_TRUE((ddc::in_tags_v<Patch6, Connectivity::all_patches>));
    EXPECT_TRUE((ddc::in_tags_v<Patch7, Connectivity::all_patches>));
    EXPECT_TRUE((ddc::in_tags_v<Patch8, Connectivity::all_patches>));
    EXPECT_TRUE((ddc::in_tags_v<Patch9, Connectivity::all_patches>));
}

TEST(MultipatchConnectivityTest, GetConnections)
{
    using namespace periodic_strips_uniform_2d_9patches;
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
    EXPECT_TRUE((ddc::in_tags_v<NorthInterface1, Patch1Interfaces>));
    EXPECT_TRUE((ddc::in_tags_v<Interface_1_2, Patch1Interfaces>));
    EXPECT_TRUE((ddc::in_tags_v<Interface_1_4, Patch1Interfaces>));
    EXPECT_TRUE((ddc::in_tags_v<Interface_3_1, Patch1Interfaces>));

    using Patch5Interfaces = typename Connectivity::template get_type_seq_connections_t<Patch5>;
    EXPECT_TRUE((ddc::in_tags_v<Interface_2_5, Patch5Interfaces>));
    EXPECT_TRUE((ddc::in_tags_v<Interface_5_6, Patch5Interfaces>));
    EXPECT_TRUE((ddc::in_tags_v<Interface_5_8, Patch5Interfaces>));
    EXPECT_TRUE((ddc::in_tags_v<Interface_4_5, Patch5Interfaces>));
}

TEST(MultipatchConnectivityTest, FindConnections)
{
    using namespace periodic_strips_uniform_2d_9patches;
    EXPECT_TRUE((std::is_same_v<
                 typename Connectivity::template find_connections_t<Patch1, Patch2>,
                 ddc::detail::TypeSeq<Interface_1_2>>));
    EXPECT_TRUE((std::is_same_v<
                 typename Connectivity::template find_connections_t<Patch1, Patch5>,
                 ddc::detail::TypeSeq<>>));
    EXPECT_TRUE((std::is_same_v<
                 typename Connectivity::template find_connections_t<Patch7, Patch9>,
                 ddc::detail::TypeSeq<Interface_9_7>>));
}

TEST(MultipatchConnectivityTest, FindConnections2)
{
    using namespace figure_of_eight_5patches;
    using Interfaces_1_2 = typename Connectivity::template find_connections_t<Patch1, Patch2>;
    EXPECT_TRUE((ddc::in_tags_v<LoopInterface_2_1, Interfaces_1_2>));
    EXPECT_TRUE((ddc::in_tags_v<EightInterface_2_1, Interfaces_1_2>));
}

TEST(MultipatchConnectivityTest, GetAllInterfacesAlongDim)
{
    using namespace periodic_strips_uniform_2d_9patches;
    using StartGridNonPeriodic = typename Patch2::Grid2;

    using InterfaceTypeSeqY
            = Connectivity::template get_all_interfaces_along_direction_t<StartGridNonPeriodic>;
    static_assert(ddc::type_seq_size_v<InterfaceTypeSeqY> == 4);
    static_assert((std::is_same_v<
                   ddc::type_seq_element_t<0, InterfaceTypeSeqY>,
                   Interface<OutsideEdge, SouthEdge<8>, true>>));
    static_assert((std::is_same_v<
                   ddc::type_seq_element_t<1, InterfaceTypeSeqY>,
                   Interface<NorthEdge<8>, SouthEdge<5>, true>>));
    static_assert((std::is_same_v<
                   ddc::type_seq_element_t<2, InterfaceTypeSeqY>,
                   Interface<NorthEdge<5>, SouthEdge<2>, true>>));
    static_assert((std::is_same_v<
                   ddc::type_seq_element_t<3, InterfaceTypeSeqY>,
                   Interface<NorthEdge<2>, OutsideEdge, true>>));

    using StartGridPeriodic = typename Patch2::Grid1;

    using InterfaceTypeSeqX
            = Connectivity::template get_all_interfaces_along_direction_t<StartGridPeriodic>;
    static_assert(ddc::type_seq_size_v<InterfaceTypeSeqX> == 3);
    static_assert((ddc::in_tags_v<Interface_1_2, InterfaceTypeSeqX>));
    static_assert((ddc::in_tags_v<Interface_2_3, InterfaceTypeSeqX>));
    static_assert((ddc::in_tags_v<Interface<EastEdge<3>, WestEdge<1>, true>, InterfaceTypeSeqX>));
}

TEST(MultipatchConnectivityTest, GetAllIndexRangesAlongDim)
{
    using namespace periodic_strips_uniform_2d_9patches;
    Patch1::IdxRange12 idx_range_1 = build_idx_range<Patch1>(0, 5, 9, 10);
    Patch2::IdxRange12 idx_range_2 = build_idx_range<Patch2>(1, 5, 8, 10);
    Patch3::IdxRange12 idx_range_3 = build_idx_range<Patch3>(2, 5, 7, 10);
    Patch4::IdxRange12 idx_range_4 = build_idx_range<Patch4>(3, 5, 6, 10);
    Patch5::IdxRange12 idx_range_5 = build_idx_range<Patch5>(4, 5, 5, 10);
    Patch6::IdxRange12 idx_range_6 = build_idx_range<Patch6>(5, 5, 4, 10);
    Patch7::IdxRange12 idx_range_7 = build_idx_range<Patch7>(6, 5, 3, 10);
    Patch8::IdxRange12 idx_range_8 = build_idx_range<Patch8>(7, 5, 2, 10);
    Patch9::IdxRange12 idx_range_9 = build_idx_range<Patch9>(8, 5, 1, 10);
    std::tuple all_idx_ranges
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

    auto idx_ranges_y
            = Connectivity::template get_all_idx_ranges_along_direction<StartGridNonPeriodic>(
                    all_idx_ranges);
    static_assert(std::tuple_size_v<decltype(idx_ranges_y)> == 3);
    EXPECT_EQ(std::get<0>(idx_ranges_y), ddc::select<GridY<8>>(idx_range_8));
    EXPECT_EQ(std::get<1>(idx_ranges_y), ddc::select<GridY<5>>(idx_range_5));
    EXPECT_EQ(std::get<2>(idx_ranges_y), ddc::select<GridY<2>>(idx_range_2));

    using StartGridPeriodic = typename Patch2::Grid1;

    auto idx_ranges_x
            = Connectivity::template get_all_idx_ranges_along_direction<StartGridPeriodic>(
                    all_idx_ranges);
    static_assert(std::tuple_size_v<decltype(idx_ranges_x)> == 3);
    EXPECT_EQ(std::get<IdxRange<GridX<1>>>(idx_ranges_x), ddc::select<GridX<1>>(idx_range_1));
    EXPECT_EQ(std::get<IdxRange<GridX<2>>>(idx_ranges_x), ddc::select<GridX<2>>(idx_range_2));
    EXPECT_EQ(std::get<IdxRange<GridX<3>>>(idx_ranges_x), ddc::select<GridX<3>>(idx_range_3));
}

TEST(MultipatchConnectivityTest, GetAllIndexRangesAlongDimSimple)
{
    using namespace non_periodic_uniform_2d_2patches;
    Patch1::IdxRange12 idx_range_1 = build_idx_range<Patch1>(0, 5, 9, 10);
    Patch2::IdxRange12 idx_range_2 = build_idx_range<Patch2>(1, 5, 8, 10);
    std::tuple all_idx_ranges = {idx_range_1, idx_range_2};

    using StartGridY = typename Patch2::Grid2;

    auto idx_ranges_y
            = Connectivity::template get_all_idx_ranges_along_direction<StartGridY>(all_idx_ranges);
    static_assert(std::tuple_size_v<decltype(idx_ranges_y)> == 1);
    EXPECT_EQ(std::get<0>(idx_ranges_y), ddc::select<GridY<2>>(idx_range_2));

    using StartGridX = typename Patch2::Grid1;

    auto idx_ranges_x
            = Connectivity::template get_all_idx_ranges_along_direction<StartGridX>(all_idx_ranges);
    static_assert(std::tuple_size_v<decltype(idx_ranges_x)> == 2);
    static_assert(std::is_same_v<
                  typename Patch1::IdxRange1,
                  std::tuple_element_t<0, decltype(idx_ranges_x)>>);
    static_assert(std::is_same_v<
                  typename Patch2::IdxRange1,
                  std::tuple_element_t<1, decltype(idx_ranges_x)>>);
    EXPECT_EQ(std::get<0>(idx_ranges_x), ddc::select<GridX<1>>(idx_range_1));
    EXPECT_EQ(std::get<1>(idx_ranges_x), ddc::select<GridX<2>>(idx_range_2));
}

TEST(MultipatchConnectivityTest, GetAllIndexRangesAlongFigureOfEightDim)
{
    using namespace figure_of_eight_5patches;
    Patch1::IdxRange12 idx_range_1 = build_idx_range<Patch1>(0, 5, 9, 10);
    Patch2::IdxRange12 idx_range_2 = build_idx_range<Patch2>(1, 5, 8, 10);
    Patch3::IdxRange12 idx_range_3 = build_idx_range<Patch3>(2, 5, 7, 10);
    Patch4::IdxRange12 idx_range_4 = build_idx_range<Patch4>(3, 5, 6, 10);
    Patch5::IdxRange12 idx_range_5 = build_idx_range<Patch5>(4, 5, 5, 10);
    std::tuple all_idx_ranges = {idx_range_1, idx_range_2, idx_range_3, idx_range_4, idx_range_5};

    using MiddleGridX = typename Patch3::Grid1;

    auto figure_eight_grid = Connectivity::template get_all_idx_ranges_along_direction<MiddleGridX>(
            all_idx_ranges);
    static_assert(std::tuple_size_v<decltype(figure_eight_grid)> == 6);
    EXPECT_EQ(std::get<IdxRange<GridX<3>>>(figure_eight_grid), ddc::select<GridX<3>>(idx_range_3));
    EXPECT_EQ(std::get<IdxRange<GridX<4>>>(figure_eight_grid), ddc::select<GridX<4>>(idx_range_4));
    EXPECT_EQ(std::get<IdxRange<GridY<5>>>(figure_eight_grid), ddc::select<GridY<5>>(idx_range_5));
    EXPECT_EQ(std::get<IdxRange<GridY<3>>>(figure_eight_grid), ddc::select<GridY<3>>(idx_range_3));
    EXPECT_EQ(std::get<IdxRange<GridY<1>>>(figure_eight_grid), ddc::select<GridY<1>>(idx_range_1));
    EXPECT_EQ(std::get<IdxRange<GridX<2>>>(figure_eight_grid), ddc::select<GridX<2>>(idx_range_2));
}

TEST(MultipatchConnectivityTest, GetAllIndexRangesAlongFigureOfEightDimMultipatchType)
{
    using namespace figure_of_eight_5patches;
    Patch1::IdxRange12 idx_range_1 = build_idx_range<Patch1>(0, 5, 9, 10);
    Patch2::IdxRange12 idx_range_2 = build_idx_range<Patch2>(1, 5, 8, 10);
    Patch3::IdxRange12 idx_range_3 = build_idx_range<Patch3>(2, 5, 7, 10);
    Patch4::IdxRange12 idx_range_4 = build_idx_range<Patch4>(3, 5, 6, 10);
    Patch5::IdxRange12 idx_range_5 = build_idx_range<Patch5>(4, 5, 5, 10);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2, Patch3, Patch4, Patch5>
            all_domains(idx_range_1, idx_range_2, idx_range_3, idx_range_4, idx_range_5);

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

TEST(MultipatchConnectivityTest, GetAllInterfacesAlongFigureOfEightDim)
{
    using namespace figure_of_eight_5patches;

    using MiddleGridX = typename Patch3::Grid1;

    using InterfaceTypeSeq
            = Connectivity::template get_all_interfaces_along_direction_t<MiddleGridX>;
    static_assert(ddc::type_seq_size_v<InterfaceTypeSeq> == 6);
    static_assert((ddc::in_tags_v<Interface<EastEdge<3>, WestEdge<4>, true>, InterfaceTypeSeq>));
    static_assert((ddc::in_tags_v<Interface<EastEdge<4>, SouthEdge<5>, false>, InterfaceTypeSeq>));
    static_assert((ddc::in_tags_v<Interface<NorthEdge<5>, SouthEdge<3>, true>, InterfaceTypeSeq>));
    static_assert((ddc::in_tags_v<Interface<NorthEdge<3>, SouthEdge<1>, true>, InterfaceTypeSeq>));
    static_assert((ddc::in_tags_v<Interface<NorthEdge<1>, WestEdge<2>, false>, InterfaceTypeSeq>));
    static_assert((ddc::in_tags_v<Interface<EastEdge<2>, WestEdge<3>, true>, InterfaceTypeSeq>));
}

TEST(MultipatchConnectivityDetailsTest, FilterEdges)
{
    using namespace non_periodic_uniform_2d_2patches;
    using FilteredEdges = Connectivity::inner_edges;
    static_assert(ddc::in_tags_v<NorthEdge<1>, FilteredEdges>);
    static_assert(ddc::in_tags_v<SouthEdge<1>, FilteredEdges>);
    static_assert(ddc::in_tags_v<EastEdge<1>, FilteredEdges>);
    static_assert(ddc::in_tags_v<WestEdge<1>, FilteredEdges>);
    static_assert(ddc::in_tags_v<NorthEdge<2>, FilteredEdges>);
    static_assert(ddc::in_tags_v<SouthEdge<2>, FilteredEdges>);
    static_assert(ddc::in_tags_v<EastEdge<2>, FilteredEdges>);
    static_assert(ddc::in_tags_v<WestEdge<2>, FilteredEdges>);
    static_assert(ddc::type_seq_size_v<FilteredEdges> == 8);
}

TEST(MultipatchConnectivityDetailsTest, ExtractPatches)
{
    using namespace non_periodic_uniform_2d_2patches;
    using ExtractedPatches = Connectivity::all_patches;
    static_assert(ddc::type_seq_size_v<ExtractedPatches> == 2);
    static_assert(ddc::in_tags_v<Patch1, ExtractedPatches>);
    static_assert(ddc::in_tags_v<Patch2, ExtractedPatches>);
}

TEST(MultipatchConnectivityDetailsTest, FindPatch)
{
    using namespace periodic_strips_uniform_2d_9patches;
    using GridX1Patch = find_patch_t<GridX<1>, Connectivity::all_patches>;
    static_assert(std::is_same_v<GridX1Patch, Patch1>);
    using GridY7Patch = find_patch_t<GridY<7>, Connectivity::all_patches>;
    static_assert(std::is_same_v<GridY7Patch, Patch7>);
}

TEST(MultipatchConnectivityDetailsTest, FindAssociatedInterface)
{
    using namespace non_periodic_uniform_2d_2patches;
    using NorthEdge1Interface
            = find_associated_interface_t<NorthEdge<1>, Connectivity::interface_collection>;
    static_assert(std::is_same_v<NorthEdge1Interface, NorthInterface1>);
    using SouthEdge1Interface
            = find_associated_interface_t<SouthEdge<1>, Connectivity::interface_collection>;
    static_assert(std::is_same_v<SouthEdge1Interface, SouthInterface1>);
    using WestEdge2Interface
            = find_associated_interface_t<WestEdge<2>, Connectivity::interface_collection>;
    static_assert(std::is_same_v<WestEdge2Interface, Interface_1_2>);
}

TEST(MultipatchConnectivityDetailsTest, CollectGridsOnDim)
{
    using namespace non_periodic_non_uniform_2d_3patches;

    // INTERFACES ------------------------------------------------------------------------------------
    /* It is supposed to represent:
     * | 1 | 2 | 3 |
     */
    using NorthInterface1 = Interface<NorthEdge<1>, OutsideEdge, true>;
    using NorthInterface2 = Interface<NorthEdge<2>, OutsideEdge, true>;
    using NorthInterface3 = Interface<NorthEdge<3>, OutsideEdge, true>;

    using SouthInterface1 = Interface<OutsideEdge, SouthEdge<1>, true>;
    using SouthInterface2 = Interface<OutsideEdge, SouthEdge<2>, true>;
    using SouthInterface3 = Interface<OutsideEdge, SouthEdge<3>, true>;

    using WestInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
    using EastInterface3 = Interface<EastEdge<3>, OutsideEdge, true>;

    using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;
    using Interface_2_3 = Interface<EastEdge<2>, WestEdge<3>, true>;


    // CONNECTIVITY ----------------------------------------------------------------------------------
    using Connectivity = MultipatchConnectivity<
            NorthInterface1,
            NorthInterface2,
            NorthInterface3,
            SouthInterface1,
            SouthInterface2,
            SouthInterface3,
            WestInterface1,
            EastInterface3,
            Interface_1_2,
            Interface_2_3>;

    using interface_collection =
            typename Connectivity::get_all_interfaces_along_direction_t<GridX<1>>;
    using all_patches = typename Connectivity::all_patches;

    using expected_interfaces_along_gridx
            = ddc::detail::TypeSeq<WestInterface1, Interface_1_2, Interface_2_3, EastInterface3>;

    static_assert(std::is_same_v<interface_collection, expected_interfaces_along_gridx>);

    // Try from patch 1
    using Grid1DSeq1 = collect_grids_on_dim_t<Patch1, GridX<1>, interface_collection>;

    static_assert(std::is_same_v<Grid1DSeq1, ddc::detail::TypeSeq<GridX<1>, GridX<2>, GridX<3>>>);
    static_assert(ddc::in_tags_v<GridX<1>, Grid1DSeq1>);
    static_assert(ddc::in_tags_v<GridX<2>, Grid1DSeq1>);
    static_assert(ddc::in_tags_v<GridX<3>, Grid1DSeq1>);

    // Try from patch 2
    using Grid1DSeq2 = collect_grids_on_dim_t<Patch2, GridX<2>, interface_collection>;

    static_assert(std::is_same_v<Grid1DSeq2, ddc::detail::TypeSeq<GridX<1>, GridX<2>, GridX<3>>>);
    static_assert(ddc::in_tags_v<GridX<1>, Grid1DSeq2>);
    static_assert(ddc::in_tags_v<GridX<2>, Grid1DSeq2>);
    static_assert(ddc::in_tags_v<GridX<3>, Grid1DSeq2>);

    // Try from patch 3
    using Grid1DSeq3 = collect_grids_on_dim_t<Patch3, GridX<3>, interface_collection>;

    static_assert(std::is_same_v<Grid1DSeq3, ddc::detail::TypeSeq<GridX<1>, GridX<2>, GridX<3>>>);
    static_assert(ddc::in_tags_v<GridX<1>, Grid1DSeq3>);
    static_assert(ddc::in_tags_v<GridX<2>, Grid1DSeq3>);
    static_assert(ddc::in_tags_v<GridX<3>, Grid1DSeq3>);
}
