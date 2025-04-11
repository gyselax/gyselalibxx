// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "3patches_2d_non_periodic_non_uniform.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "mesh_builder.hpp"
#include "patch.hpp"

using namespace non_periodic_non_uniform_2d_3patches;


namespace {
class IndexTransformationNonUniformTest : public ::testing::Test
{
private:
    static constexpr Patch1::IdxStep1 x1_size = Patch1::IdxStep1(16 + 1);
    static constexpr Patch1::IdxStep2 y1_size = Patch1::IdxStep2(10 + 1);

    static constexpr Patch2::IdxStep1 x2_size = Patch2::IdxStep1(8 + 1);
    static constexpr Patch2::IdxStep2 y2_size = Patch2::IdxStep2(12 + 1);

    static constexpr Patch3::IdxStep1 x3_size = Patch3::IdxStep1(8 + 1 + 2);
    static constexpr Patch3::IdxStep2 y3_size = Patch3::IdxStep2(12 + 1 + 3);

protected:
    Patch1::IdxRange1 const idx_range_x1;
    Patch1::IdxRange2 const idx_range_y1;

    Patch2::IdxRange1 const idx_range_x2;
    Patch2::IdxRange2 const idx_range_y2;

    Patch3::IdxRange1 const idx_range_x3;
    Patch3::IdxRange2 const idx_range_y3;

public:
    IndexTransformationNonUniformTest()
        : idx_range_x1(Patch1::Idx1(0), x1_size)
        , idx_range_y1(Patch1::Idx2(0), y1_size)
        , idx_range_x2(Patch2::Idx1(0), x2_size)
        , idx_range_y2(Patch2::Idx2(0), y2_size)
        , idx_range_x3(Patch3::Idx1(2), x3_size - 2)
        , idx_range_y3(Patch3::Idx2(3), y3_size - 3)
    {
    }


    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        Patch1::Coord1 const x1_min(0.0);
        Patch1::Coord1 const x1_max(2.0);

        Patch1::Coord2 const y1_min(2.5);
        Patch1::Coord2 const y1_max(7.0);

        std::vector<Patch1::Coord1> break_points_x1
                = build_uniform_break_points(x1_min, x1_max, x1_size - 1);
        std::vector<Patch1::Coord2> break_points_y1
                = build_uniform_break_points(y1_min, y1_max, y1_size - 1);

        ddc::init_discrete_space<GridX<1>>(break_points_x1);
        ddc::init_discrete_space<GridY<1>>(break_points_y1);

        // Patch 2
        Patch2::Coord1 const x2_min(1.0);
        Patch2::Coord1 const x2_max(3.0);

        Patch2::Coord2 const y2_min(-4.0);
        Patch2::Coord2 const y2_max(-3.5);

        std::vector<Patch2::Coord1> break_points_x2
                = build_uniform_break_points(x2_min, x2_max, x2_size - 1);
        std::vector<Patch2::Coord2> break_points_y2
                = build_uniform_break_points(y2_min, y2_max, y2_size - 1);

        ddc::init_discrete_space<GridX<2>>(break_points_x2);
        ddc::init_discrete_space<GridY<2>>(break_points_y2);

        // Patch 3
        Patch3::Coord1 const x3_min(1.0);
        Patch3::Coord1 const x3_max(3.0);

        Patch3::Coord2 const y3_min(-4.0);
        Patch3::Coord2 const y3_max(-3.5);

        std::vector<Patch3::Coord1> break_points_x3
                = build_uniform_break_points(x3_min, x3_max, x3_size - 1);
        std::vector<Patch3::Coord2> break_points_y3
                = build_uniform_break_points(y3_min, y3_max, y3_size - 1);

        ddc::init_discrete_space<GridX<3>>(break_points_x3);
        ddc::init_discrete_space<GridY<3>>(break_points_y3);
    }
};

} // namespace


TEST_F(IndexTransformationNonUniformTest, IndexAvailability)
{
    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeY2F = Edge<Patch2, GridY<2>, FRONT>;
    using Interface12 = Interface<EdgeY1B, EdgeY2F, false>;

    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeY3F = Edge<Patch3, GridY<3>, FRONT>;
    using Interface13 = Interface<EdgeY1B, EdgeY3F, false>;

    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeY3F = Edge<Patch3, GridY<3>, FRONT>;
    using Interface31 = Interface<EdgeY3F, EdgeY1B, false>;

    // Coordinate transformation .................................................................
    // Index starting at 0
    EdgeTransformation<Interface12> index_transformation_12(idx_range_x1, idx_range_x2);

    Patch1::Idx1 test_idx_x1(9);
    bool is_available = index_transformation_12.is_match_available(test_idx_x1);
    EXPECT_FALSE(is_available);

    test_idx_x1 = Patch1::Idx1(8);
    is_available = index_transformation_12.is_match_available(test_idx_x1);
    EXPECT_TRUE(is_available);

    // Target index not starting at 0
    EdgeTransformation<Interface13> index_transformation_13(idx_range_x1, idx_range_x3);

    test_idx_x1 = Patch1::Idx1(9);
    is_available = index_transformation_13.is_match_available(test_idx_x1);
    EXPECT_FALSE(is_available);

    test_idx_x1 = Patch1::Idx1(8);
    is_available = index_transformation_13.is_match_available(test_idx_x1);
    EXPECT_TRUE(is_available);

    // Current index not starting at 0
    EdgeTransformation<Interface31> index_transformation_31(idx_range_x3, idx_range_x1);

    test_idx_x1 = Patch1::Idx1(8);
    is_available = index_transformation_31.is_match_available(test_idx_x1);
    EXPECT_TRUE(is_available);
}


TEST_F(IndexTransformationNonUniformTest, InvertedOrientation)
{
    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeY2F = Edge<Patch2, GridY<2>, FRONT>;
    using Interface12 = Interface<EdgeY1B, EdgeY2F, false>;

    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeY3F = Edge<Patch3, GridY<3>, FRONT>;
    using Interface13 = Interface<EdgeY1B, EdgeY3F, false>;

    // Coordinate transformation .................................................................
    // Index starting at 0
    EdgeTransformation<Interface12> index_transformation_12(idx_range_x1, idx_range_x2);

    Patch1::Idx1 test_idx_x1(12);
    Patch2::Idx1 test_idx_x2(index_transformation_12(test_idx_x1));

    EXPECT_EQ(test_idx_x2, Patch2::Idx1(2));

    // Target index not starting at 0
    EdgeTransformation<Interface13> index_transformation_13(idx_range_x1, idx_range_x3);

    test_idx_x1 = Patch1::Idx1(12);
    Patch3::Idx1 test_idx_x3(index_transformation_13(test_idx_x1));

    EXPECT_EQ(test_idx_x3, Patch3::Idx1(2 + 2));
}



TEST_F(IndexTransformationNonUniformTest, StickingDifferentDimensions)
{
    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeX2F = Edge<Patch2, GridX<2>, FRONT>;
    using Interface12 = Interface<EdgeY1B, EdgeX2F, true>;

    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeX3F = Edge<Patch3, GridX<3>, FRONT>;
    using Interface13 = Interface<EdgeY1B, EdgeX3F, true>;

    // Coordinate transformation .................................................................
    // Index starting at 0
    EdgeTransformation<Interface12> index_transformation_12(idx_range_x1, idx_range_y2);

    Patch1::Idx1 test_idx_x1(4);
    Patch2::Idx2 test_idx_y2(index_transformation_12(test_idx_x1));

    EXPECT_EQ(test_idx_y2, Patch2::Idx2(3));

    // Target index not starting at 0
    EdgeTransformation<Interface13> index_transformation_13(idx_range_x1, idx_range_y3);

    test_idx_x1 = Patch1::Idx1(4);
    Patch3::Idx2 test_idx_y3(index_transformation_13(test_idx_x1));

    // The first index of the index range is 3. 
    EXPECT_EQ(test_idx_y3, Patch3::Idx2(3) +  Patch3::IdxStep2(3));
}



TEST_F(IndexTransformationNonUniformTest, ReverseTransformation)
{
    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeY2F = Edge<Patch2, GridY<2>, FRONT>;
    using Interface12 = Interface<EdgeY1B, EdgeY2F, false>;

    using EdgeY1B = Edge<Patch1, GridY<1>, BACK>;
    using EdgeY3F = Edge<Patch3, GridY<3>, FRONT>;
    using Interface13 = Interface<EdgeY1B, EdgeY3F, false>;

    // Coordinate transformation .................................................................
    EdgeTransformation<Interface12> index_transformation_12(idx_range_x1, idx_range_x2);

    Patch2::Idx1 test_idx_x2(2);
    Patch1::Idx1 test_idx_x1(index_transformation_12(test_idx_x2));

    EXPECT_EQ(test_idx_x1, Patch1::Idx1(12));

    EdgeTransformation<Interface13> index_transformation_13(idx_range_x1, idx_range_x3);

    Patch3::Idx1 test_idx_x3(2);
    // The first index of the index range is 2. 
    test_idx_x3 += Patch3::IdxStep1(2); 
    test_idx_x1 = index_transformation_13(test_idx_x3);

    EXPECT_EQ(test_idx_x1, Patch1::Idx1(12));
}
