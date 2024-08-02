// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_non_periodic_uniform.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "patch.hpp"


namespace {
class CoordinateTransformationTest : public ::testing::Test
{
private:
    static constexpr Patch1::IdxStep1 x1_size = Patch1::IdxStep1(10);
    static constexpr Patch1::IdxStep2 y1_size = Patch1::IdxStep2(10);

    static constexpr Patch2::IdxStep1 x2_size = Patch2::IdxStep1(10);
    static constexpr Patch2::IdxStep2 y2_size = Patch2::IdxStep2(10);

protected:
    Patch1::IdxRange1 const idx_range_x1;
    Patch1::IdxRange2 const idx_range_y1;

    Patch2::IdxRange1 const idx_range_x2;
    Patch2::IdxRange2 const idx_range_y2;

public:
    CoordinateTransformationTest()
        : idx_range_x1(Patch1::Idx1(0), x1_size)
        , idx_range_y1(Patch1::Idx2(0), y1_size)
        , idx_range_x2(Patch2::Idx1(0), x2_size)
        , idx_range_y2(Patch2::Idx2(0), y2_size) {};

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        Patch1::Coord1 const x1_min(0.0);
        Patch1::Coord1 const x1_max(2.0);

        Patch1::Coord2 const y1_min(2.5);
        Patch1::Coord2 const y1_max(7.0);

        ddc::init_discrete_space<GridX1>(GridX1::init(x1_min, x1_max, x1_size));
        ddc::init_discrete_space<GridY1>(GridY1::init(y1_min, y1_max, y1_size));

        // Patch 2
        Patch2::Coord1 const x2_min(1.0);
        Patch2::Coord1 const x2_max(3.0);

        Patch2::Coord2 const y2_min(-4.0);
        Patch2::Coord2 const y2_max(-3.5);

        ddc::init_discrete_space<GridX2>(GridX2::init(x2_min, x2_max, x2_size));
        ddc::init_discrete_space<GridY2>(GridY2::init(y2_min, y2_max, y2_size));
    }
};

} // namespace



TEST_F(CoordinateTransformationTest, InvertedOrientation)
{
    using EgdeX1B = Edge<Patch1, GridX1, BACK>;
    using EgdeX2F = Edge<Patch2, GridX2, FRONT>;
    using Interface12 = Interface<EgdeX1B, EgdeX2F, false>;

    // Coordinate transformation .................................................................
    EdgeTransformation<Interface12> coord_transformation(idx_range_x1, idx_range_x2);

    Patch1::Coord1 test_coord_x1(1.5);
    Patch2::Coord1 test_coord_x2(coord_transformation(test_coord_x1));

    EXPECT_NEAR(double(test_coord_x2), 1.5, 1e-14);
}



TEST_F(CoordinateTransformationTest, StickingDifferentDimensions)
{
    using EgdeX1B = Edge<Patch1, GridX1, BACK>;
    using EgdeY2F = Edge<Patch2, GridY2, FRONT>;
    using Interface12 = Interface<EgdeX1B, EgdeY2F, true>;

    // Coordinate transformation .................................................................
    EdgeTransformation<Interface12> coord_transformation(idx_range_x1, idx_range_y2);

    Patch1::Coord1 test_coord_x1(0.7);
    Patch2::Coord2 test_coord_y2(coord_transformation(test_coord_x1));

    EXPECT_NEAR(double(test_coord_y2), -3.825, 1e-14);
}



TEST_F(CoordinateTransformationTest, ReverseTransformation)
{
    using EgdeX1B = Edge<Patch1, GridX1, BACK>;
    using EgdeX2F = Edge<Patch2, GridX2, FRONT>;
    using Interface12 = Interface<EgdeX1B, EgdeX2F, false>;

    // Coordinate transformation .................................................................
    EdgeTransformation<Interface12> coord_transformation(idx_range_x1, idx_range_x2);

    Patch2::Coord1 test_coord_x2(1.5);
    Patch1::Coord1 test_coord_x1 = coord_transformation(test_coord_x2);

    EXPECT_NEAR(double(test_coord_x1), 1.5, 1e-14);
}