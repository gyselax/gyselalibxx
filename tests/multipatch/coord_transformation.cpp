// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "coord_transformation.hpp"
#include "edge.hpp"
#include "interface.hpp"

namespace {
// Continuous dimension of patch 1
struct RDimX1
{
    static bool constexpr PERIODIC = false;
};

struct RDimY1
{
    static bool constexpr PERIODIC = false;
};

// Continuous dimension of patch 2
struct RDimX2
{
    static bool constexpr PERIODIC = false;
};

struct RDimY2
{
    static bool constexpr PERIODIC = false;
};


using CoordX1 = ddc::Coordinate<RDimX1>;
using CoordY1 = ddc::Coordinate<RDimY1>;
using CoordX2 = ddc::Coordinate<RDimX2>;
using CoordY2 = ddc::Coordinate<RDimY2>;

// Discrete dimensions
struct IDimX1 : ddc::UniformPointSampling<RDimX1>
{
};
struct IDimY1 : ddc::UniformPointSampling<RDimY1>
{
};
struct IDimX2 : ddc::UniformPointSampling<RDimX2>
{
};
struct IDimY2 : ddc::UniformPointSampling<RDimY2>
{
};

using IndexX1 = ddc::DiscreteElement<IDimX1>;
using IndexY1 = ddc::DiscreteElement<IDimY1>;
using IndexX2 = ddc::DiscreteElement<IDimX2>;
using IndexY2 = ddc::DiscreteElement<IDimY2>;


using IVectX1 = ddc::DiscreteVector<IDimX1>;
using IVectY1 = ddc::DiscreteVector<IDimY1>;
using IVectX2 = ddc::DiscreteVector<IDimX2>;
using IVectY2 = ddc::DiscreteVector<IDimY2>;

using IDomainX1 = ddc::DiscreteDomain<IDimX1>;
using IDomainY1 = ddc::DiscreteDomain<IDimY1>;
using IDomainX2 = ddc::DiscreteDomain<IDimX2>;
using IDomainY2 = ddc::DiscreteDomain<IDimY2>;


class CoordinateTransformationTest : public ::testing::Test
{
private:
    static constexpr IVectX1 x1_size = IVectX1(10);
    static constexpr IVectY1 y1_size = IVectY1(10);

    static constexpr IVectX2 x2_size = IVectX2(10);
    static constexpr IVectY2 y2_size = IVectY2(10);

protected:
    IDomainX1 const domain_x1;
    IDomainY1 const domain_y1;

    IDomainX2 const domain_x2;
    IDomainY2 const domain_y2;

public:
    CoordinateTransformationTest()
        : domain_x1(IndexX1(0), x1_size)
        , domain_y1(IndexY1(0), y1_size)
        , domain_x2(IndexX2(0), x2_size)
        , domain_y2(IndexY2(0), y2_size) {};

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        CoordX1 const x1_min(0.0);
        CoordX1 const x1_max(2.0);

        CoordY1 const y1_min(2.5);
        CoordY1 const y1_max(7.0);

        ddc::init_discrete_space<IDimX1>(IDimX1::init(x1_min, x1_max, x1_size));
        ddc::init_discrete_space<IDimY1>(IDimY1::init(y1_min, y1_max, y1_size));

        // Patch 2
        CoordX2 const x2_min(1.0);
        CoordX2 const x2_max(3.0);

        CoordY2 const y2_min(-4.0);
        CoordY2 const y2_max(-3.5);

        ddc::init_discrete_space<IDimX2>(IDimX2::init(x2_min, x2_max, x2_size));
        ddc::init_discrete_space<IDimY2>(IDimY2::init(y2_min, y2_max, y2_size));
    }
};


} // namespace



TEST_F(CoordinateTransformationTest, InvertedOrientation)
{
    // --- TEST 1 ---
    // Creating interfaces .......................................................................
    Edge<IDimX1> edge_1 = {1, domain_x1, BACK};
    Edge<IDimX2> edge_2 = {2, domain_x2, FRONT};

    Interface<IDimX1, IDimX2> interface = {edge_1, edge_2, false};

    // Coordinate transformation .................................................................
    EdgeCoordinatesTransformation coord_transformation(interface);

    CoordX1 test_coord_x1(1.5);
    CoordX2 test_coord_x2(coord_transformation(test_coord_x1));

    EXPECT_NEAR(double(test_coord_x2), 1.5, 1e-14);
}


TEST_F(CoordinateTransformationTest, StickingDifferentDimensions)
{
    // --- TEST 2 ---
    // Creating interfaces .......................................................................
    Edge<IDimX1> edge_1 = {1, domain_x1, BACK};
    Edge<IDimY2> edge_2 = {2, domain_y2, FRONT};

    Interface<IDimX1, IDimY2> interface = {edge_1, edge_2, true};

    // Coordinate transformation .................................................................
    EdgeCoordinatesTransformation coord_transformation(interface);

    CoordX1 test_coord_x1(0.7);
    CoordY2 test_coord_y2(coord_transformation(test_coord_x1));

    EXPECT_NEAR(double(test_coord_y2), -3.825, 1e-14);
}


TEST_F(CoordinateTransformationTest, ReverseTransformation)
{
    // --- TEST 3 ---
    // Creating interfaces .......................................................................
    Edge<IDimX1> edge_1 = {1, domain_x1, BACK};
    Edge<IDimX2> edge_2 = {2, domain_x2, FRONT};

    Interface<IDimX1, IDimX2> interface = {edge_1, edge_2, false};

    // Coordinate transformation .................................................................
    EdgeCoordinatesTransformation coord_transformation(interface);

    CoordX2 test_coord_x2(1.5);
    CoordX1 test_coord_x1 = coord_transformation(test_coord_x2);

    EXPECT_NEAR(double(test_coord_x1), 1.5, 1e-14);
}