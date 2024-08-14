// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_non_periodic_uniform.hpp"
#include "edge.hpp"
#include "idx_range_slice.hpp"
#include "interface.hpp"
#include "matching_idx_slice.hpp"
#include "patch.hpp"


using namespace non_periodic_uniform_2d_2patches;

namespace {
class MatchingIdxSliceUniformGridTest : public ::testing::Test
{
protected:
    static constexpr Patch1::IdxStep1 x1_size = Patch1::IdxStep1(16 + 1);
    static constexpr Patch1::IdxStep2 y1_size = Patch1::IdxStep2(10 + 1 + 1);

    static constexpr Patch2::IdxStep1 x2_size = Patch2::IdxStep1(8 + 1 + 2);
    static constexpr Patch2::IdxStep2 y2_size = Patch2::IdxStep2(12 + 1 + 3);


    Patch1::IdxRange1 const idx_range_x1;
    Patch1::IdxRange2 const idx_range_y1;

    Patch2::IdxRange1 const idx_range_x2;
    Patch2::IdxRange2 const idx_range_y2;

public:
    MatchingIdxSliceUniformGridTest()
        : idx_range_x1(Patch1::Idx1(0), x1_size)
        , idx_range_y1(Patch1::Idx2(1), y1_size - 1)
        , idx_range_x2(Patch2::Idx1(2), x2_size - 2)
        , idx_range_y2(Patch2::Idx2(3), y2_size - 3) {};

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        Patch1::Coord1 const x1_min(0.0);
        Patch1::Coord1 const x1_max(2.0);

        Patch1::Coord2 const y1_min(2.5);
        Patch1::Coord2 const y1_max(7.0);

        ddc::init_discrete_space<GridX<1>>(GridX<1>::init(x1_min, x1_max, x1_size));
        ddc::init_discrete_space<GridY<1>>(GridY<1>::init(y1_min, y1_max, y1_size));


        // Patch 2
        Patch2::Coord1 const x2_min(1.0);
        Patch2::Coord1 const x2_max(3.0);

        Patch2::Coord2 const y2_min(-4.0);
        Patch2::Coord2 const y2_max(-3.5);

        ddc::init_discrete_space<GridX<2>>(GridX<2>::init(x2_min, x2_max, x2_size));
        ddc::init_discrete_space<GridY<2>>(GridY<2>::init(y2_min, y2_max, y2_size));
    }
};

} // namespace


/*
    Sticking between edges defined along 
        - GridX<1> : uniform of 16 cells; 
        - GridX<2> : uniform of 8 cells; 
    Conforming lines:
        0   2   4   6   8   10  12  14  16   
        0   1   2   3   4   5   6   7   8
*/
TEST_F(MatchingIdxSliceUniformGridTest, OneGridIncludeInTheOther)
{
    Patch1::IdxRange12 idx_range_1(idx_range_x1, idx_range_y1);
    Patch2::IdxRange12 idx_range_2(idx_range_x2, idx_range_y2);

    using SouthEdge1 = Edge<Patch1, GridY<1>, FRONT>;
    using NorthEdge2 = Edge<Patch2, GridY<2>, BACK>;

    using Interface_12 = Interface<SouthEdge1, NorthEdge2, true>;

    MatchingIdxSlice<Interface_12> idx_matching(idx_range_1, idx_range_2);

    // Get the conforming index ranges.
    IdxRangeSlice<GridX<1>> conforming_idx_range_x1 = idx_matching.get<GridX<1>>();
    IdxRangeSlice<GridX<2>> conforming_idx_range_x2 = idx_matching.get<GridX<2>>();

    // Compare
    // --- parameters for patch 1
    EXPECT_EQ(conforming_idx_range_x1.front(), idx_range_x1.front());
    EXPECT_EQ(conforming_idx_range_x1.extents(), Patch1::IdxStep1(idx_range_x2.size()));
    EXPECT_EQ(conforming_idx_range_x1.strides(), Patch1::IdxStep1(2));

    // --- parameters for patch 2
    EXPECT_EQ(conforming_idx_range_x2.front(), idx_range_x2.front());
    EXPECT_EQ(conforming_idx_range_x2.extents(), idx_range_x2.extents());
    EXPECT_EQ(conforming_idx_range_x2.strides(), Patch2::IdxStep1(1));
};



/*
    Sticking between edges defined along 
        - GridY<1> : uniform of 10 cells; 
        - GridY<2> : uniform of 12 cells; 
    Conforming lines:
        0   5   10
        0   6   12
*/
TEST_F(MatchingIdxSliceUniformGridTest, NoGridIncludeInTheOther)
{
    Patch1::IdxRange12 idx_range_1(idx_range_x1, idx_range_y1);
    Patch2::IdxRange12 idx_range_2(idx_range_x2, idx_range_y2);

    using WestEdge1 = Edge<Patch1, GridX<1>, BACK>;
    using EastEdge2 = Edge<Patch2, GridX<2>, FRONT>;

    using Interface_12 = Interface<WestEdge1, EastEdge2, true>;

    MatchingIdxSlice<Interface_12> idx_matching(idx_range_1, idx_range_2);

    // Get the conforming index ranges.
    IdxRangeSlice<GridY<1>> conforming_idx_range_y1 = idx_matching.get<GridY<1>>();
    IdxRangeSlice<GridY<2>> conforming_idx_range_y2 = idx_matching.get<GridY<2>>();

    // Compare
    // --- parameters for patch 1
    EXPECT_EQ(conforming_idx_range_y1.front(), idx_range_y1.front());
    EXPECT_EQ(conforming_idx_range_y1.extents(), Patch1::IdxStep2(3));
    EXPECT_EQ(conforming_idx_range_y1.strides(), Patch1::IdxStep2(5));

    // --- parameters for patch 2
    EXPECT_EQ(conforming_idx_range_y2.front(), idx_range_y2.front());
    EXPECT_EQ(conforming_idx_range_y2.extents(), Patch2::IdxStep2(3));
    EXPECT_EQ(conforming_idx_range_y2.strides(), Patch2::IdxStep2(6));
};



/// Check that the order of the edge in the interface does not matter.
TEST_F(MatchingIdxSliceUniformGridTest, OrderOfEdgeDoesntMatter)
{
    Patch1::IdxRange12 idx_range_1(idx_range_x1, idx_range_y1);
    Patch2::IdxRange12 idx_range_2(idx_range_x2, idx_range_y2);

    using SouthEdge1 = Edge<Patch1, GridY<1>, FRONT>;
    using NorthEdge2 = Edge<Patch2, GridY<2>, BACK>;

    using Interface_12 = Interface<NorthEdge2, SouthEdge1, true>;

    MatchingIdxSlice<Interface_12> idx_matching(idx_range_2, idx_range_1);

    // Get the conforming index ranges.
    IdxRangeSlice<GridX<1>> conforming_idx_range_x1 = idx_matching.get<GridX<1>>();
    IdxRangeSlice<GridX<2>> conforming_idx_range_x2 = idx_matching.get<GridX<2>>();

    // Compare
    // --- parameters for patch 1
    EXPECT_EQ(conforming_idx_range_x1.front(), idx_range_x1.front());
    EXPECT_EQ(conforming_idx_range_x1.extents(), Patch1::IdxStep1(idx_range_x2.size()));
    EXPECT_EQ(conforming_idx_range_x1.strides(), Patch1::IdxStep1(2));

    // --- parameters for patch 2
    EXPECT_EQ(conforming_idx_range_x2.front(), idx_range_x2.front());
    EXPECT_EQ(conforming_idx_range_x2.extents(), idx_range_x2.extents());
    EXPECT_EQ(conforming_idx_range_x2.strides(), Patch2::IdxStep1(1));
};


/// Check the get_from_perp function.
TEST_F(MatchingIdxSliceUniformGridTest, GetFromPerpTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_x1, idx_range_y1);
    Patch2::IdxRange12 idx_range_2(idx_range_x2, idx_range_y2);

    using SouthEdge1 = Edge<Patch1, GridY<1>, FRONT>;
    using NorthEdge2 = Edge<Patch2, GridY<2>, BACK>;

    using Interface_12 = Interface<SouthEdge1, NorthEdge2, true>;

    MatchingIdxSlice<Interface_12> idx_matching(idx_range_1, idx_range_2);

    // Get the conforming index ranges.
    IdxRangeSlice<GridX<1>> conforming_idx_range_x1 = idx_matching.get_from_perp<GridY<1>>();
    IdxRangeSlice<GridX<2>> conforming_idx_range_x2 = idx_matching.get_from_perp<GridY<2>>();

    // Compare
    // --- parameters for patch 1
    EXPECT_EQ(conforming_idx_range_x1.front(), idx_range_x1.front());
    EXPECT_EQ(conforming_idx_range_x1.extents(), Patch1::IdxStep1(idx_range_x2.size()));
    EXPECT_EQ(conforming_idx_range_x1.strides(), Patch1::IdxStep1(2));

    // --- parameters for patch 2
    EXPECT_EQ(conforming_idx_range_x2.front(), idx_range_x2.front());
    EXPECT_EQ(conforming_idx_range_x2.extents(), idx_range_x2.extents());
    EXPECT_EQ(conforming_idx_range_x2.strides(), Patch2::IdxStep1(1));
};


// Check we can access the indexes.
TEST_F(MatchingIdxSliceUniformGridTest, AccessingIndexes)
{
    Patch1::IdxRange12 idx_range_1(idx_range_x1, idx_range_y1);
    Patch2::IdxRange12 idx_range_2(idx_range_x2, idx_range_y2);

    using SouthEdge1 = Edge<Patch1, GridY<1>, FRONT>;
    using NorthEdge2 = Edge<Patch2, GridY<2>, BACK>;

    using Interface_12 = Interface<SouthEdge1, NorthEdge2, true>;

    MatchingIdxSlice<Interface_12> idx_matching(idx_range_1, idx_range_2);

    // Get the conforming index ranges.
    IdxRangeSlice<GridX<1>> conforming_idx_range_x1 = idx_matching.get<GridX<1>>();
    IdxRangeSlice<GridX<2>> conforming_idx_range_x2 = idx_matching.get<GridX<2>>();

    // Set the expected conforming index ranges.
    IdxRangeSlice<GridX<1>> expected_conforming_idx_range_x1(
            idx_range_x1.front(),
            Patch1::IdxStep1(idx_range_x2.size()),
            Patch1::IdxStep1(2));
    IdxRangeSlice<GridX<2>> expected_conforming_idx_range_x2(
            idx_range_x2.front(),
            idx_range_x2.extents(),
            Patch2::IdxStep1(1));

    // Check if all the elements in expected conforming idx range are in the computed conforming idx range
    for (Patch1::Idx1 const& idx : expected_conforming_idx_range_x1) {
        EXPECT_TRUE(conforming_idx_range_x1.contains(idx));
    };

    for (Patch2::Idx1 const& idx : expected_conforming_idx_range_x2) {
        EXPECT_TRUE(conforming_idx_range_x2.contains(idx));
    };

    // Check if all the elements in computed conforming idx range are in the expected conforming idx range
    // and are available.
    for (Patch1::Idx1 const& idx : conforming_idx_range_x1) {
        EXPECT_TRUE(expected_conforming_idx_range_x1.contains(idx));
    };

    for (Patch2::Idx1 const& idx : expected_conforming_idx_range_x2) {
        EXPECT_TRUE(expected_conforming_idx_range_x2.contains(idx));
    };
};
