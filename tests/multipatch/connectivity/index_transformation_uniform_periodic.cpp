// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_uniform.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "patch.hpp"

using namespace onion_shape_uniform_2d_2patches;

namespace {
class IndexTransformationUniformPeriodicTest : public ::testing::Test
{
private:
    using SplineInterpPointsTheta1 = ddc::KnotsAsInterpolationPoints<
            BSplinesTheta<1>,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC>;
    using SplineInterpPointsTheta2 = ddc::KnotsAsInterpolationPoints<
            BSplinesTheta<2>,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC>;

    static constexpr Patch1::IdxStep1 r1_npoints = Patch1::IdxStep1(16 + 1);
    static constexpr Patch1::IdxStep2 theta1_npoints = Patch1::IdxStep2(8);

    static constexpr Patch2::IdxStep1 r2_npoints = Patch2::IdxStep1(8 + 1);
    static constexpr Patch2::IdxStep2 theta2_npoints = Patch2::IdxStep2(12);

protected:
    Patch1::IdxRange1 const idx_range_r1;
    Patch1::IdxRange2 const idx_range_theta1;

    Patch2::IdxRange1 const idx_range_r2;
    Patch2::IdxRange2 const idx_range_theta2;

public:
    IndexTransformationUniformPeriodicTest()
        : idx_range_r1(Patch1::Idx1(0), r1_npoints)
        , idx_range_theta1(Patch1::Idx2(0), theta1_npoints)
        , idx_range_r2(Patch2::Idx1(0), r2_npoints)
        , idx_range_theta2(Patch2::Idx2(1), theta2_npoints)
    {
    }

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        Patch1::Coord1 const r1_min(0.0);
        Patch1::Coord1 const r1_max(2.0);

        Patch1::Coord2 const theta1_min(2.5);
        Patch1::Coord2 const theta1_max(7.0);

        ddc::init_discrete_space<BSplinesTheta<1>>(theta1_min, theta1_max, theta1_npoints);

        ddc::init_discrete_space<GridR<1>>(GridR<1>::init(r1_min, r1_max, r1_npoints));
        ddc::init_discrete_space<GridTheta<1>>(
                SplineInterpPointsTheta1::get_sampling<GridTheta<1>>());

        // Patch 2
        Patch2::Coord1 const r2_min(1.0);
        Patch2::Coord1 const r2_max(3.0);

        Patch2::Coord2 const theta2_min(-4.0);
        Patch2::Coord2 const theta2_max(-3.5);

        ddc::init_discrete_space<BSplinesTheta<2>>(theta2_min, theta2_max, theta2_npoints);

        ddc::init_discrete_space<GridR<2>>(GridR<2>::init(r2_min, r2_max, r2_npoints));
        ddc::init_discrete_space<GridTheta<2>>(
                SplineInterpPointsTheta2::get_sampling<GridTheta<2>>());
    }
};

} // namespace


TEST_F(IndexTransformationUniformPeriodicTest, IndexAvailability)
{
    using EdgeR1B = Edge<Patch1, GridR<1>, BACK>;
    using EdgeR2F = Edge<Patch2, GridR<2>, FRONT>;
    using Interface12 = Interface<EdgeR1B, EdgeR2F, false>;

    // Coordinate transformation .................................................................
    // Current index starting at 0
    EdgeTransformation<Interface12> index_transformation(idx_range_theta1, idx_range_theta2);

    Patch1::Idx2 test_idx_theta1(7);
    bool is_available = index_transformation.is_match_available(test_idx_theta1);
    EXPECT_FALSE(is_available);

    test_idx_theta1 = Patch1::Idx2(6);
    is_available = index_transformation.is_match_available(test_idx_theta1);
    EXPECT_TRUE(is_available);
}


TEST_F(IndexTransformationUniformPeriodicTest, InvertedOrientation)
{
    using EdgeR1B = Edge<Patch1, GridR<1>, BACK>;
    using EdgeR2F = Edge<Patch2, GridR<2>, FRONT>;
    using Interface12 = Interface<EdgeR1B, EdgeR2F, false>;

    // Coordinate transformation .................................................................
    EdgeTransformation<Interface12> index_transformation(idx_range_theta1, idx_range_theta2);

    Patch1::Idx2 test_idx_theta1(6);
    Patch2::Idx2 test_idx_theta2(index_transformation(test_idx_theta1));

    EXPECT_EQ((test_idx_theta2 - idx_range_theta2.front()).value(), 3);
}


TEST_F(IndexTransformationUniformPeriodicTest, ReverseTransformation)
{
    using EdgeR1B = Edge<Patch1, GridR<1>, BACK>;
    using EdgeR2F = Edge<Patch2, GridR<2>, FRONT>;
    using Interface12 = Interface<EdgeR1B, EdgeR2F, false>;

    // Coordinate transformation .................................................................
    EdgeTransformation<Interface12> index_transformation(idx_range_theta1, idx_range_theta2);

    Patch2::Idx2 test_idx_theta2(3 + 1);
    Patch1::Idx2 test_idx_theta1(index_transformation(test_idx_theta2));

    EXPECT_EQ((test_idx_theta1 - idx_range_theta1.front()).value(), 6);
}


TEST_F(IndexTransformationUniformPeriodicTest, PeriodicitySameOrientation)
{
    using EdgeR1B = Edge<Patch1, GridR<1>, BACK>;
    using EdgeR2F = Edge<Patch2, GridR<2>, FRONT>;
    using Interface12 = Interface<EdgeR1B, EdgeR2F, true>;

    // Coordinate transformation .................................................................
    EdgeTransformation<Interface12> index_transformation(idx_range_theta1, idx_range_theta2);

    Patch1::Idx2 test_idx_theta1(8);
    Patch2::Idx2 test_idx_theta2 = index_transformation(test_idx_theta1);

    EXPECT_EQ((test_idx_theta2 - idx_range_theta2.front()).value(), 0);
}


TEST_F(IndexTransformationUniformPeriodicTest, PeriodicityInvertedOrientation)
{
    using EdgeR1B = Edge<Patch1, GridR<1>, BACK>;
    using EdgeR2F = Edge<Patch2, GridR<2>, FRONT>;
    using Interface12 = Interface<EdgeR1B, EdgeR2F, false>;

    // Coordinate transformation .................................................................
    EdgeTransformation<Interface12> index_transformation(idx_range_theta1, idx_range_theta2);

    Patch1::Idx2 test_idx_theta1(0);
    Patch2::Idx2 test_idx_theta2 = index_transformation(test_idx_theta1);

    EXPECT_EQ((test_idx_theta2 - idx_range_theta2.front()).value(), 0);
}