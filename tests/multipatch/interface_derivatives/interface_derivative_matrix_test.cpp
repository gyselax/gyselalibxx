// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

// #include "2patches_2d_onion_shape_non_uniform.hpp"
// #include "../../test_utils.hpp"

#include "9patches_2d_periodic_strips_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
// #include "global_2d_onion_shape_non_uniform.hpp"


#include <typeinfo>

#include <boost/type_index.hpp>

#include "interface.hpp"
#include "interface_derivative_matrix.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "interface_derivatives_test_utils.hpp"


/*
    Test: 
    * 3 patches
    * compute derivatives x1, x2, x1x2? (or just function + derivatives?)
    * which interpolation methods to choose? Additional interpolation points? 
*/

namespace {
// Multi-patch ---
using namespace periodic_strips_non_uniform_2d_9patches;

// Equivalent global mesh ---
struct Xg
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Xg;
};

struct Yg
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Yg;
};

struct GridXg : NonUniformGridBase<Xg>
{
};
struct GridYg : NonUniformGridBase<Yg>
{
};

struct BSplinesXg : ddc::NonUniformBSplines<Xg, 3>
{
};
struct BSplinesYg : ddc::NonUniformBSplines<Yg, 3>
{
};


using HostExecSpace = Kokkos::DefaultHostExecutionSpace;


// Interpolation points type for the patches
template <std::size_t PatchIdx>
using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
        BSplinesX<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

template <std::size_t PatchIdx, ddc::BoundCond BoundCondMin, ddc::BoundCond BoundCondMax>
using SplineInterpPointsY
        = ddcHelper::NonUniformInterpolationPoints<BSplinesY<1>, BoundCondMin, BoundCondMax>;

// Interpolation points type for the equivalent global spline.
using SplineInterpPointsXg = ddcHelper::NonUniformInterpolationPoints<
        BSplinesXg,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;
using SplineInterpPointsYg = ddcHelper::NonUniformInterpolationPoints<
        BSplinesYg,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE>;


// Operators on the equivalent global spline
using SplineRThetagBuilder = ddc::SplineBuilder2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::SplineSolver::LAPACK>;

using SplineRThetagEvaluator = ddc::SplineEvaluator2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::PeriodicExtrapolationRule<Xg>,
        ddc::PeriodicExtrapolationRule<Xg>,
        ddc::ConstantExtrapolationRule<Yg, Xg>,
        ddc::ConstantExtrapolationRule<Yg, Xg>>;


// TODO: Add the functions defined here in another file. 

// /**
//  *  @brief Get interpolation points from the break points by placing 
//  * the interpolation points on the break points and adding one on the
//  * left boundary cell at 2/3 of the cell. 
//  */
// template <class CoordType>
// std::vector<CoordType> get_interpolation_points_add_one_on_left(
//         std::vector<CoordType> const& break_points)
// {
//     CoordType additional_point(break_points[0] * 2. / 3. + break_points[1] * 1. / 3.);
//     std::vector<CoordType> interpolation_points(break_points);
//     interpolation_points.insert(interpolation_points.begin() + 1, additional_point);
//     return interpolation_points;
// }

// /**
//  *  @brief Get interpolation points from the break points by placing 
//  * the interpolation points on the break points and adding one on the
//  * right boundary cell at 1/3 of the cell. 
//  */
// template <class CoordType>
// std::vector<CoordType> get_interpolation_points_add_one_on_right(
//         std::vector<CoordType> const& break_points)
// {
//     int n_bpoints = break_points.size();
//     CoordType additional_point(
//             break_points[n_bpoints - 1] * 2. / 3. + break_points[n_bpoints - 2] * 1. / 3.);
//     std::vector<CoordType> interpolation_points(break_points);
//     interpolation_points.insert(interpolation_points.end() - 1, additional_point);
//     return interpolation_points;
// }

// /**
//  * @brief Fill in a vector of points for the equivalent global mesh 
//  * by conserving the same order of the given points.
//  */
// template <class CoordTypeG, class CoordTypeP>
// void fill_in(std::vector<CoordTypeG>& points_global, std::vector<CoordTypeP> const& points_patch)
// {
//     for (CoordTypeP pt : points_patch) {
//         points_global.push_back(CoordTypeG {double(pt)});
//     }
// }

// /**
//  * @brief Fill in a vector of points for the equivalent global mesh
//  *  by reversing the order of the given points.
//  */
// template <class CoordTypeG, class CoordTypeP>
// void fill_in_reverse(
//         std::vector<CoordTypeG>& points_global,
//         std::vector<CoordTypeP> const& points_patch)
// {
//     std::size_t const n_pt = points_patch.size();
//     CoordTypeP const max = points_patch[n_pt - 1];
//     CoordTypeP const min = points_patch[0];
//     for (int i(0); i < n_pt; ++i) {
//         points_global.push_back(CoordTypeG {double(min + max - points_patch[n_pt - 1 - i])});
//     }
// }

/// @brief Initialise the function with f(r,theta) = r(3-r)sin(theta).
template <class Grid1, class Grid2>
void initialise_2D_function(host_t<DField<IdxRange<Grid1, Grid2>>> function)
{
    ddc::for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        double const rg = ddc::coordinate(Idx<Grid1>(idx));
        double const thetag = ddc::coordinate(Idx<Grid2>(idx));
        function(idx) = rg * (3. - rg) * Kokkos::sin(thetag);
    });
}



struct InterfaceDerivativeMatrixTest : public ::testing::Test
{
    // DEFINE BOUNDARIES OF THE DOMAINS ----------------------------------------------------------
    // patches 1 | 4 | 7  dim X ------------------
    static constexpr Coord<X> x1_min = Coord<X>(0.0);
    static constexpr Coord<X> x1_max = Coord<X>(1.0);
    static constexpr IdxStep<GridX<1>> x1_ncells = IdxStep<GridX<1>>(10);

    static constexpr Coord<X> x4_min = x1_min;
    static constexpr Coord<X> x4_max = x1_max;
    static constexpr IdxStep<GridX<4>> x4_ncells = IdxStep<GridX<4>>(10);

    static constexpr Coord<X> x7_min = x1_min;
    static constexpr Coord<X> x7_max = x1_max;
    static constexpr IdxStep<GridX<7>> x7_ncells = IdxStep<GridX<7>>(10);

    // patches 2 | 5 | 8  dim X ------------------
    static constexpr Coord<X> x2_min = Coord<X>(1.0);
    static constexpr Coord<X> x2_max = Coord<X>(2.0);
    static constexpr IdxStep<GridX<2>> x2_ncells = IdxStep<GridX<2>>(10);

    static constexpr Coord<X> x5_min = x2_min;
    static constexpr Coord<X> x5_max = x2_max;
    static constexpr IdxStep<GridX<5>> x5_ncells = IdxStep<GridX<5>>(10);

    static constexpr Coord<X> x8_min = x2_min;
    static constexpr Coord<X> x8_max = x2_max;
    static constexpr IdxStep<GridX<8>> x8_ncells = IdxStep<GridX<8>>(10);

    // patches 3 | 6 | 9  dim X ------------------
    static constexpr Coord<X> x3_min = Coord<X>(2.0);
    static constexpr Coord<X> x3_max = Coord<X>(3.0);
    static constexpr IdxStep<GridX<3>> x3_ncells = IdxStep<GridX<3>>(10);

    static constexpr Coord<X> x6_min = x3_min;
    static constexpr Coord<X> x6_max = x3_max;
    static constexpr IdxStep<GridX<6>> x6_ncells = IdxStep<GridX<6>>(10);

    static constexpr Coord<X> x9_min = x3_min;
    static constexpr Coord<X> x9_max = x3_max;
    static constexpr IdxStep<GridX<9>> x9_ncells = IdxStep<GridX<9>>(10);

    // patches 1 | 2 | 3  dim Y --------------------
    static constexpr Coord<Y> y1_min = Coord<Y>(2.0);
    static constexpr Coord<Y> y1_max = Coord<Y>(3.0);
    static constexpr IdxStep<GridY<1>> y1_ncells = IdxStep<GridY<1>>(10);

    static constexpr Coord<Y> y2_min = y1_min;
    static constexpr Coord<Y> y2_max = y1_max;
    static constexpr IdxStep<GridY<2>> y2_ncells = IdxStep<GridY<2>>(10);

    static constexpr Coord<Y> y3_min = y1_min;
    static constexpr Coord<Y> y3_max = y1_max;
    static constexpr IdxStep<GridY<3>> y3_ncells = IdxStep<GridY<3>>(10);

    // patches 4 | 5 | 6  dim Y --------------------
    static constexpr Coord<Y> y4_min = Coord<Y>(1.0);
    static constexpr Coord<Y> y4_max = Coord<Y>(2.0);
    static constexpr IdxStep<GridY<4>> y4_ncells = IdxStep<GridY<4>>(10);

    static constexpr Coord<Y> y5_min = y4_min;
    static constexpr Coord<Y> y5_max = y4_max;
    static constexpr IdxStep<GridY<5>> y5_ncells = IdxStep<GridY<5>>(10);

    static constexpr Coord<Y> y6_min = y4_min;
    static constexpr Coord<Y> y6_max = y4_max;
    static constexpr IdxStep<GridY<6>> y6_ncells = IdxStep<GridY<6>>(10);

    // patches 7 | 8 | 9  dim Y --------------------
    static constexpr Coord<Y> y7_min = Coord<Y>(0.0);
    static constexpr Coord<Y> y7_max = Coord<Y>(1.0);
    static constexpr IdxStep<GridY<7>> y7_ncells = IdxStep<GridY<7>>(10);

    static constexpr Coord<Y> y8_min = y7_min;
    static constexpr Coord<Y> y8_max = y7_max;
    static constexpr IdxStep<GridY<8>> y8_ncells = IdxStep<GridY<8>>(10);

    static constexpr Coord<Y> y9_min = y7_min;
    static constexpr Coord<Y> y9_max = y7_max;
    static constexpr IdxStep<GridY<9>> y9_ncells = IdxStep<GridY<9>>(10);


    // global ------------------------------------
    static constexpr Coord<Xg> xg_min = Coord<Xg> {double(x1_min)};
    static constexpr Coord<Xg> xg_max = Coord<Xg> {double(x3_max)};
    static constexpr IdxStep<GridXg> xg_ncells
            = IdxStep<GridXg>(x1_ncells.value() + x2_ncells.value() + x3_ncells.value());

    static constexpr Coord<Yg> yg_min = Coord<Yg> {double(x1_min)};
    static constexpr Coord<Yg> yg_max = Coord<Yg> {double(x3_max)};
    static constexpr IdxStep<GridYg> yg_ncells
            = IdxStep<GridYg>(y1_ncells.value() + y4_ncells.value() + y7_ncells.value());

    static constexpr ddc::BoundCond BcH = ddc::BoundCond::HERMITE;
    static constexpr ddc::BoundCond BcG = ddc::BoundCond::GREVILLE;

protected:
    const IdxRange<GridX<1>> idx_range_x1;
    const IdxRange<GridX<2>> idx_range_x2;
    const IdxRange<GridX<3>> idx_range_x3;
    const IdxRange<GridX<4>> idx_range_x4;
    const IdxRange<GridX<5>> idx_range_x5;
    const IdxRange<GridX<6>> idx_range_x6;
    const IdxRange<GridX<7>> idx_range_x7;
    const IdxRange<GridX<8>> idx_range_x8;
    const IdxRange<GridX<9>> idx_range_x9;


    const IdxRange<GridY<1>> idx_range_y1;
    const IdxRange<GridY<2>> idx_range_y2;
    const IdxRange<GridY<3>> idx_range_y3;
    const IdxRange<GridY<4>> idx_range_y4;
    const IdxRange<GridY<5>> idx_range_y5;
    const IdxRange<GridY<6>> idx_range_y6;
    const IdxRange<GridY<7>> idx_range_y7;
    const IdxRange<GridY<8>> idx_range_y8;
    const IdxRange<GridY<9>> idx_range_y9;


    const IdxRange<GridX<1>, GridY<1>> idx_range_xy1;
    const IdxRange<GridX<2>, GridY<2>> idx_range_xy2;
    const IdxRange<GridX<3>, GridY<3>> idx_range_xy3;
    const IdxRange<GridX<4>, GridY<4>> idx_range_xy4;
    const IdxRange<GridX<5>, GridY<5>> idx_range_xy5;
    const IdxRange<GridX<6>, GridY<6>> idx_range_xy6;
    const IdxRange<GridX<7>, GridY<7>> idx_range_xy7;
    const IdxRange<GridX<8>, GridY<8>> idx_range_xy8;
    const IdxRange<GridX<9>, GridY<9>> idx_range_xy9;

public:
    InterfaceDerivativeMatrixTest()
        : idx_range_x1(SplineInterpPointsX<1>::template get_domain<GridX<1>>())
        , idx_range_x2(SplineInterpPointsX<2>::template get_domain<GridX<2>>())
        , idx_range_x3(SplineInterpPointsX<3>::template get_domain<GridX<3>>())
        , idx_range_x4(SplineInterpPointsX<4>::template get_domain<GridX<4>>())
        , idx_range_x5(SplineInterpPointsX<5>::template get_domain<GridX<5>>())
        , idx_range_x6(SplineInterpPointsX<6>::template get_domain<GridX<6>>())
        , idx_range_x7(SplineInterpPointsX<7>::template get_domain<GridX<7>>())
        , idx_range_x8(SplineInterpPointsX<8>::template get_domain<GridX<8>>())
        , idx_range_x9(SplineInterpPointsX<9>::template get_domain<GridX<9>>())
        , idx_range_y1(SplineInterpPointsY<1, BcH, BcG>::template get_domain<GridY<1>>())
        , idx_range_y2(SplineInterpPointsY<2, BcH, BcG>::template get_domain<GridY<2>>())
        , idx_range_y3(SplineInterpPointsY<3, BcH, BcG>::template get_domain<GridY<3>>())
        , idx_range_y4(SplineInterpPointsY<4, BcH, BcH>::template get_domain<GridY<4>>())
        , idx_range_y5(SplineInterpPointsY<5, BcH, BcH>::template get_domain<GridY<5>>())
        , idx_range_y6(SplineInterpPointsY<6, BcH, BcH>::template get_domain<GridY<6>>())
        , idx_range_y7(SplineInterpPointsY<7, BcG, BcH>::template get_domain<GridY<7>>())
        , idx_range_y8(SplineInterpPointsY<8, BcG, BcH>::template get_domain<GridY<8>>())
        , idx_range_y9(SplineInterpPointsY<9, BcG, BcH>::template get_domain<GridY<9>>())
        , idx_range_xy1(idx_range_x1, idx_range_y1)
        , idx_range_xy2(idx_range_x2, idx_range_y2)
        , idx_range_xy3(idx_range_x3, idx_range_y3)
        , idx_range_xy4(idx_range_x4, idx_range_y4)
        , idx_range_xy5(idx_range_x5, idx_range_y5)
        , idx_range_xy6(idx_range_x6, idx_range_y6)
        , idx_range_xy7(idx_range_x7, idx_range_y7)
        , idx_range_xy8(idx_range_x8, idx_range_y8)
        , idx_range_xy9(idx_range_x9, idx_range_y9)
    {
    }


    // INITIALISE DOMAINS ------------------------------------------------------------------------
    static void SetUpTestSuite()
    {
        // Creating of meshes and supports .......................................................
        std::vector<Coord<X>> break_points_x147
                = build_random_non_uniform_break_points(x1_min, x1_max, x1_ncells);
        std::vector<Coord<X>> break_points_x258
                = build_random_non_uniform_break_points(x2_min, x2_max, x2_ncells);
        std::vector<Coord<X>> break_points_x369
                = build_random_non_uniform_break_points(x3_min, x3_max, x3_ncells);

        std::vector<Coord<Y>> break_points_y123
                = build_random_non_uniform_break_points(y1_min, y1_max, y1_ncells);
        std::vector<Coord<Y>> break_points_y456
                = build_random_non_uniform_break_points(y4_min, y4_max, y4_ncells);
        std::vector<Coord<Y>> break_points_y789
                = build_random_non_uniform_break_points(y7_min, y7_max, y7_ncells);

        std::vector<Coord<X>> interpolation_points_x147 = break_points_x147;
        std::vector<Coord<X>> interpolation_points_x258 = break_points_x258;
        std::vector<Coord<X>> interpolation_points_x369 = break_points_x369;

        std::vector<Coord<Y>> interpolation_points_y123
                = get_interpolation_points_add_one_on_right(break_points_y123);
        std::vector<Coord<Y>> interpolation_points_y456 = break_points_y456;
        std::vector<Coord<Y>> interpolation_points_y789
                = get_interpolation_points_add_one_on_left(break_points_y789);



        // Patch 1 ...............................................................................
        ddc::init_discrete_space<BSplinesX<1>>(break_points_x147);
        ddc::init_discrete_space<BSplinesY<1>>(break_points_y123);

        ddc::init_discrete_space<GridX<1>>(interpolation_points_x147);
        ddc::init_discrete_space<GridY<1>>(interpolation_points_y123);

        // Patch 2 ...............................................................................
        ddc::init_discrete_space<BSplinesX<2>>(break_points_x258);
        ddc::init_discrete_space<BSplinesY<2>>(break_points_y123);

        ddc::init_discrete_space<GridX<2>>(interpolation_points_x258);
        ddc::init_discrete_space<GridY<2>>(interpolation_points_y123);

        // Patch 3 ...............................................................................
        ddc::init_discrete_space<BSplinesX<3>>(break_points_x369);
        ddc::init_discrete_space<BSplinesY<3>>(break_points_y123);

        ddc::init_discrete_space<GridX<3>>(interpolation_points_x369);
        ddc::init_discrete_space<GridY<3>>(interpolation_points_y123);

        // Patch 4 ...............................................................................
        ddc::init_discrete_space<BSplinesX<4>>(break_points_x147);
        ddc::init_discrete_space<BSplinesY<4>>(break_points_y456);

        ddc::init_discrete_space<GridX<4>>(interpolation_points_x147);
        ddc::init_discrete_space<GridY<4>>(interpolation_points_y456);

        // Patch 5 ...............................................................................
        ddc::init_discrete_space<BSplinesX<5>>(break_points_x258);
        ddc::init_discrete_space<BSplinesY<5>>(break_points_y456);

        ddc::init_discrete_space<GridX<5>>(interpolation_points_x258);
        ddc::init_discrete_space<GridY<5>>(interpolation_points_y456);

        // Patch 6 ...............................................................................
        ddc::init_discrete_space<BSplinesX<6>>(break_points_x369);
        ddc::init_discrete_space<BSplinesY<6>>(break_points_y456);

        ddc::init_discrete_space<GridX<6>>(interpolation_points_x369);
        ddc::init_discrete_space<GridY<6>>(interpolation_points_y456);

        // Patch 7 ...............................................................................
        ddc::init_discrete_space<BSplinesX<7>>(break_points_x147);
        ddc::init_discrete_space<BSplinesY<7>>(break_points_y789);

        ddc::init_discrete_space<GridX<7>>(interpolation_points_x147);
        ddc::init_discrete_space<GridY<7>>(interpolation_points_y789);

        // Patch 8 ...............................................................................
        ddc::init_discrete_space<BSplinesX<8>>(break_points_x258);
        ddc::init_discrete_space<BSplinesY<8>>(break_points_y789);

        ddc::init_discrete_space<GridX<8>>(interpolation_points_x258);
        ddc::init_discrete_space<GridY<8>>(interpolation_points_y789);

        // Patch 9 ...............................................................................
        ddc::init_discrete_space<BSplinesX<9>>(break_points_x369);
        ddc::init_discrete_space<BSplinesY<9>>(break_points_y789);

        ddc::init_discrete_space<GridX<9>>(interpolation_points_x369);
        ddc::init_discrete_space<GridY<9>>(interpolation_points_y789);


        // Equivalent global domain ..............................................................
        std::vector<Coord<Xg>> break_points_xg;
        std::vector<Coord<Xg>> interpolation_points_xg;

        break_points_x147.pop_back();
        interpolation_points_x147.pop_back();
        fill_in(break_points_xg, break_points_x147);
        fill_in(interpolation_points_xg, interpolation_points_x147);

        break_points_x258.pop_back();
        interpolation_points_x258.pop_back();
        fill_in(break_points_xg, break_points_x258);
        fill_in(interpolation_points_xg, interpolation_points_x258);

        fill_in(break_points_xg, break_points_x369);
        fill_in(interpolation_points_xg, interpolation_points_x369);


        std::vector<Coord<Yg>> break_points_yg;
        std::vector<Coord<Yg>> interpolation_points_yg;


        break_points_y789.pop_back();
        interpolation_points_y789.pop_back();
        fill_in(break_points_yg, break_points_y789);
        fill_in(interpolation_points_yg, interpolation_points_y789);

        break_points_y456.pop_back();
        interpolation_points_y456.pop_back();
        fill_in(break_points_yg, break_points_y456);
        fill_in(interpolation_points_yg, interpolation_points_y456);

        fill_in(break_points_yg, break_points_y123);
        fill_in(interpolation_points_yg, interpolation_points_y123);


        ddc::init_discrete_space<BSplinesXg>(break_points_xg);
        ddc::init_discrete_space<BSplinesYg>(break_points_yg);

        ddc::init_discrete_space<GridXg>(interpolation_points_xg);
        ddc::init_discrete_space<GridYg>(interpolation_points_yg);
    }


    // TEST OPERATORS ----------------------------------------------------------------------------
};

} // end namespace



// Check that the local grids and the equivalent global grid match together.
TEST_F(InterfaceDerivativeMatrixTest, InterpolationPointsCheck)
{
    // Check for Patch 1 ---
    ddc::for_each(idx_range_xy1, [&](Patch1::Idx12 const& idx) {
        Patch1::IdxStep1 idx_x(Patch1::Idx1(idx) - idx_range_x1.front());
        Patch1::IdxStep2 idx_y(Patch1::Idx2(idx) - idx_range_y1.front());
        Idx<GridXg, GridYg>
                idx_g(idx_x.value(), idx_y.value() + y7_ncells.value() + y4_ncells.value() + 1);
        EXPECT_NEAR(ddc::coordinate(Patch1::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch1::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });

    // Check for Patch 2 ---
    ddc::for_each(idx_range_xy2, [&](Patch2::Idx12 const& idx) {
        Patch2::IdxStep1 idx_x(Patch2::Idx1(idx) - idx_range_x2.front());
        Patch2::IdxStep2 idx_y(Patch2::Idx2(idx) - idx_range_y2.front());
        Idx<GridXg, GridYg>
                idx_g(idx_x.value() + x1_ncells.value(),
                      idx_y.value() + y7_ncells.value() + y4_ncells.value() + 1);
        EXPECT_NEAR(ddc::coordinate(Patch2::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch2::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });

    // Check for Patch 3 ---
    ddc::for_each(idx_range_xy3, [&](Patch3::Idx12 const& idx) {
        Patch3::IdxStep1 idx_x(Patch3::Idx1(idx) - idx_range_x3.front());
        Patch3::IdxStep2 idx_y(Patch3::Idx2(idx) - idx_range_y3.front());
        Idx<GridXg, GridYg>
                idx_g(idx_x.value() + x1_ncells.value() + x2_ncells.value(),
                      idx_y.value() + y7_ncells.value() + y4_ncells.value() + 1);
        EXPECT_NEAR(ddc::coordinate(Patch3::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch3::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });

    // Check for Patch 4 ---
    ddc::for_each(idx_range_xy4, [&](Patch4::Idx12 const& idx) {
        Patch4::IdxStep1 idx_x(Patch4::Idx1(idx) - idx_range_x4.front());
        Patch4::IdxStep2 idx_y(Patch4::Idx2(idx) - idx_range_y4.front());
        Idx<GridXg, GridYg> idx_g(idx_x.value(), idx_y.value() + y4_ncells.value() + 1);
        EXPECT_NEAR(ddc::coordinate(Patch4::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch4::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });

    // Check for Patch 5 ---
    ddc::for_each(idx_range_xy5, [&](Patch5::Idx12 const& idx) {
        Patch5::IdxStep1 idx_x(Patch5::Idx1(idx) - idx_range_x5.front());
        Patch5::IdxStep2 idx_y(Patch5::Idx2(idx) - idx_range_y5.front());
        Idx<GridXg, GridYg>
                idx_g(idx_x.value() + x1_ncells.value(), idx_y.value() + y4_ncells.value() + 1);
        EXPECT_NEAR(ddc::coordinate(Patch5::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch5::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });

    // Check for Patch 6 ---
    ddc::for_each(idx_range_xy6, [&](Patch6::Idx12 const& idx) {
        Patch6::IdxStep1 idx_x(Patch6::Idx1(idx) - idx_range_x6.front());
        Patch6::IdxStep2 idx_y(Patch6::Idx2(idx) - idx_range_y6.front());
        Idx<GridXg, GridYg>
                idx_g(idx_x.value() + x1_ncells.value() + x2_ncells.value(),
                      idx_y.value() + y4_ncells.value() + 1);
        EXPECT_NEAR(ddc::coordinate(Patch6::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch6::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });

    // Check for Patch 7 ---
    ddc::for_each(idx_range_xy7, [&](Patch7::Idx12 const& idx) {
        Patch7::IdxStep1 idx_x(Patch7::Idx1(idx) - idx_range_x7.front());
        Patch7::IdxStep2 idx_y(Patch7::Idx2(idx) - idx_range_y7.front());
        Idx<GridXg, GridYg> idx_g(idx_x.value(), idx_y.value());
        EXPECT_NEAR(ddc::coordinate(Patch7::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch7::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });

    // Check for Patch 8 ---
    ddc::for_each(idx_range_xy8, [&](Patch8::Idx12 const& idx) {
        Patch8::IdxStep1 idx_x(Patch8::Idx1(idx) - idx_range_x8.front());
        Patch8::IdxStep2 idx_y(Patch8::Idx2(idx) - idx_range_y8.front());
        Idx<GridXg, GridYg> idx_g(idx_x.value() + x1_ncells.value(), idx_y.value());
        EXPECT_NEAR(ddc::coordinate(Patch8::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch8::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });

    // Check for Patch 9 ---
    ddc::for_each(idx_range_xy9, [&](Patch9::Idx12 const& idx) {
        Patch9::IdxStep1 idx_x(Patch9::Idx1(idx) - idx_range_x9.front());
        Patch9::IdxStep2 idx_y(Patch9::Idx2(idx) - idx_range_y9.front());
        Idx<GridXg, GridYg>
                idx_g(idx_x.value() + x1_ncells.value() + x2_ncells.value(), idx_y.value());
        EXPECT_NEAR(ddc::coordinate(Patch9::Idx1(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Patch9::Idx2(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });
}


TEST_F(InterfaceDerivativeMatrixTest, InterfaceDerivativeMatrixCheck)
{
    // SingleInterfaceDerivativesCalculators along y.
    SingleInterfaceDerivativesCalculator<Interface_1_2> const
            derivatives_calculator_1_2(idx_range_xy1, idx_range_xy2);

    SingleInterfaceDerivativesCalculator<Interface_2_3> const
            derivatives_calculator_2_3(idx_range_xy2, idx_range_xy3);

    SingleInterfaceDerivativesCalculator<Interface_3_1> const
            derivatives_calculator_3_1(idx_range_xy3, idx_range_xy1);

    SingleInterfaceDerivativesCalculator<Interface_4_5> const
            derivatives_calculator_4_5(idx_range_xy4, idx_range_xy5);

    SingleInterfaceDerivativesCalculator<Interface_5_6> const
            derivatives_calculator_5_6(idx_range_xy5, idx_range_xy6);

    SingleInterfaceDerivativesCalculator<Interface_6_4> const
            derivatives_calculator_6_4(idx_range_xy6, idx_range_xy4);

    SingleInterfaceDerivativesCalculator<Interface_7_8> const
            derivatives_calculator_7_8(idx_range_xy7, idx_range_xy8);

    SingleInterfaceDerivativesCalculator<Interface_8_9> const
            derivatives_calculator_8_9(idx_range_xy8, idx_range_xy9);

    SingleInterfaceDerivativesCalculator<Interface_9_7> const
            derivatives_calculator_9_7(idx_range_xy9, idx_range_xy7);

    // SingleInterfaceDerivativesCalculators along x.
    SingleInterfaceDerivativesCalculator<
            Interface_1_4,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE> const derivatives_calculator_1_4(idx_range_xy1, idx_range_xy4);

    SingleInterfaceDerivativesCalculator<
            Interface_4_7,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE> const
            derivatives_calculator_4_7(idx_range_xy4, idx_range_xy7);

    SingleInterfaceDerivativesCalculator<
            Interface_2_5,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE> const derivatives_calculator_2_5(idx_range_xy2, idx_range_xy5);

    SingleInterfaceDerivativesCalculator<
            Interface_5_8,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE> const
            derivatives_calculator_5_8(idx_range_xy5, idx_range_xy8);

    SingleInterfaceDerivativesCalculator<
            Interface_3_6,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE> const derivatives_calculator_3_6(idx_range_xy3, idx_range_xy6);

    SingleInterfaceDerivativesCalculator<
            Interface_6_9,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE> const
            derivatives_calculator_6_9(idx_range_xy6, idx_range_xy9);



    // Order in sequences.
    const std::tuple derivative_calculators_123 = std::
            tie(derivatives_calculator_1_2, derivatives_calculator_2_3, derivatives_calculator_3_1);
    const std::tuple derivative_calculators_456 = std::
            tie(derivatives_calculator_4_5, derivatives_calculator_5_6, derivatives_calculator_6_4);
    const std::tuple derivative_calculators_789 = std::
            tie(derivatives_calculator_7_8, derivatives_calculator_8_9, derivatives_calculator_9_7);

    const std::tuple derivative_calculators_147
            = std::tie(derivatives_calculator_4_7, derivatives_calculator_1_4);
    const std::tuple derivative_calculators_258
            = std::tie(derivatives_calculator_5_8, derivatives_calculator_2_5);
    const std::tuple derivative_calculators_369
            = std::tie(derivatives_calculator_6_9, derivatives_calculator_3_6);


    // Order index ranges.
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2, Patch3>
            idx_ranges_123(idx_range_xy1, idx_range_xy2, idx_range_xy3);
    MultipatchType<IdxRangeOnPatch, Patch4, Patch5, Patch6>
            idx_ranges_456(idx_range_xy4, idx_range_xy5, idx_range_xy6);
    MultipatchType<IdxRangeOnPatch, Patch7, Patch8, Patch9>
            idx_ranges_789(idx_range_xy7, idx_range_xy8, idx_range_xy9);

    MultipatchType<IdxRangeOnPatch, Patch1, Patch4, Patch7>
            idx_ranges_147(idx_range_xy1, idx_range_xy4, idx_range_xy7);
    MultipatchType<IdxRangeOnPatch, Patch2, Patch5, Patch8>
            idx_ranges_258(idx_range_xy2, idx_range_xy5, idx_range_xy8);
    MultipatchType<IdxRangeOnPatch, Patch3, Patch6, Patch9>
            idx_ranges_369(idx_range_xy3, idx_range_xy6, idx_range_xy9);


    // Define the matrix calculators.
    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<1>,
            DFieldOnPatch_host,
            ConstDeriv1_OnPatch_2D,
            true,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            Kokkos::DefaultHostExecutionSpace,
            Patch1,
            Patch2,
            Patch3>
            matrix_123(idx_ranges_123, derivative_calculators_123);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<4>,
            DFieldOnPatch_host,
            ConstDeriv1_OnPatch_2D,
            true,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            Kokkos::DefaultHostExecutionSpace,
            Patch4,
            Patch5,
            Patch6>
            matrix_456(idx_ranges_456, derivative_calculators_456);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<7>,
            DFieldOnPatch_host,
            ConstDeriv1_OnPatch_2D,
            true,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            Kokkos::DefaultHostExecutionSpace,
            Patch7,
            Patch8,
            Patch9>
            matrix_789(idx_ranges_789, derivative_calculators_789);



    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<1>,
            DFieldOnPatch_host,
            ConstDeriv1_OnPatch_2D,
            false,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            Kokkos::DefaultHostExecutionSpace,
            Patch1,
            Patch4,
            Patch7>
            matrix_147(idx_ranges_147, derivative_calculators_147);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<2>,
            DFieldOnPatch_host,
            ConstDeriv1_OnPatch_2D,
            false,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            Kokkos::DefaultHostExecutionSpace,
            Patch2,
            Patch5,
            Patch8>
            matrix_258(idx_ranges_258, derivative_calculators_258);


    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<3>,
            DFieldOnPatch_host,
            ConstDeriv1_OnPatch_2D,
            false,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            Kokkos::DefaultHostExecutionSpace,
            Patch3,
            Patch6,
            Patch9>
            matrix_369(idx_ranges_369, derivative_calculators_369);



    // using interface_collection =
    //         typename Connectivity::get_all_interfaces_along_direction_t<GridX<1>>;

    // std::cout << boost::typeindex::type_id_with_cvr<ddc::type_seq_element_t<0,interface_collection>>().pretty_name() << std::endl;
    // std::cout << boost::typeindex::type_id_with_cvr<ddc::type_seq_element_t<1,interface_collection>>().pretty_name() << std::endl;
    // std::cout << boost::typeindex::type_id_with_cvr<ddc::type_seq_element_t<2,interface_collection>>().pretty_name() << std::endl;
    // std::cout << boost::typeindex::type_id_with_cvr<ddc::type_seq_element_t<3,interface_collection>>().pretty_name() << std::endl;
}

// // Check the values of the computed interface derivatives.
// TYPED_TEST(
//         SingleInterfaceDerivativesCalculatorFixture,
//         InterfaceDerivativesExactAndApproximationFormulae)
// {
//     // Get parameters of the test.
//     constexpr ddc::BoundCond Interpolation_v = TestFixture::Interpolation_v;
//     using Edge1 = typename TestFixture::Edge1;
//     using Edge2 = typename TestFixture::Edge2;
//     using Interface_1_2 = typename TestFixture::Interface_1_2;

//     // Define EdgeTransformation operator.
//     using IdxRangePar2 = std::conditional_t<
//             (std::is_same_v<Edge2, SouthEdge2>),
//             typename Patch2::IdxRange1,
//             typename Patch2::IdxRange2>;
//     IdxRangePar2 idx_range_par2;
//     if constexpr (std::is_same_v<Edge2, SouthEdge2>) {
//         idx_range_par2 = TestFixture::idx_range_eta2;
//     } else {
//         idx_range_par2 = TestFixture::idx_range_xi2;
//     }

//     EdgeTransformation<Interface_1_2>
//             idx_convertor_12(TestFixture::idx_range_theta1, idx_range_par2);

//     // Initialise functions values ===============================================================
//     // --- patch 1
//     host_t<DFieldMem<Patch1::IdxRange12>> function_1_alloc(TestFixture::idx_range_rtheta1);
//     host_t<DField<Patch1::IdxRange12>> function_1 = get_field(function_1_alloc);

//     // --- patch 2
//     host_t<DFieldMem<Patch2::IdxRange12>> function_2_alloc(TestFixture::idx_range_etaxi2);
//     host_t<DField<Patch2::IdxRange12>> function_2 = get_field(function_2_alloc);

//     // --- global
//     host_t<DFieldMem<IdxRange<GridRg, GridThetag>>> function_g_alloc(
//             TestFixture::idx_range_rtheta_g);
//     host_t<DField<IdxRange<GridRg, GridThetag>>> function_g = get_field(function_g_alloc);

//     // Fill in with the correct value.
//     // --- patch 1
//     if constexpr (std::is_same_v<Edge1, EastEdge1>) {
//         // Same orientation than the equivalent global domain.
//         initialise_2D_function<Patch1::Grid1, Patch1::Grid2>(function_1);
//     } else {
//         // Different orientation: global: ↑→ | local: ←↓
//         ddc::for_each(get_idx_range(function_1), [&](Idx<Patch1::Grid1, Patch1::Grid2> idx) {
//             // Get the coordinate on the equivalent global domain.
//             double const rg = TestFixture::r1_max - ddc::coordinate(Idx<Patch1::Grid1>(idx));
//             double const thetag
//                     = TestFixture::theta1_max - ddc::coordinate(Idx<Patch1::Grid2>(idx));
//             function_1(idx) = rg * (3. - rg) * Kokkos::sin(thetag);
//         });
//     }

//     // --- patch 2
//     if constexpr (std::is_same_v<Edge2, WestEdge2>) {
//         // Same orientation than the equivalent global domain.
//         initialise_2D_function<Patch2::Grid1, Patch2::Grid2>(function_2);
//     } else if (std::is_same_v<Edge2, EastEdge2>) {
//         // Different orientation: global: ↑→ | local: ←↓
//         ddc::for_each(get_idx_range(function_2), [&](Idx<Patch2::Grid1, Patch2::Grid2> idx) {
//             // Get the coordinate on the equivalent global domain.
//             double const rg = TestFixture::eta2_max - ddc::coordinate(Idx<Patch2::Grid1>(idx))
//                               + TestFixture::r1_max;
//             double const thetag = TestFixture::xi2_max - ddc::coordinate(Idx<Patch2::Grid2>(idx));
//             function_2(idx) = rg * (3. - rg) * Kokkos::sin(thetag);
//         });
//     } else if (std::is_same_v<Edge2, SouthEdge2>) {
//         // Different orientation: global: ↑→ | local: ↓→
//         ddc::for_each(get_idx_range(function_2), [&](Idx<Patch2::Grid1, Patch2::Grid2> idx) {
//             // Get the coordinate on the equivalent global domain.
//             double const rg
//                     = double(ddc::coordinate(Idx<Patch2::Grid2>(idx))) + TestFixture::r1_max;
//             double const thetag
//                     = double((TestFixture::eta2_max - ddc::coordinate(Idx<Patch2::Grid1>(idx))));
//             function_2(idx) = rg * (3. - rg) * Kokkos::sin(thetag);
//         });
//     }

//     // --- global
//     initialise_2D_function<GridRg, GridThetag>(function_g);


//     // Build an equivalent global spline =========================================================
//     SplineRThetagBuilder<Interpolation_v> builder_g(TestFixture::idx_range_rtheta_g);

//     host_t<DFieldMem<IdxRange<BSplinesRg, BSplinesThetag>>> function_g_coef_alloc(
//             builder_g.batched_spline_domain(TestFixture::idx_range_rtheta_g));
//     host_t<DField<IdxRange<BSplinesRg, BSplinesThetag>>> function_g_coef
//             = get_field(function_g_coef_alloc);

//     if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
//         // --- Spline builder
//         builder_g(function_g_coef, get_const_field(function_g));
//     } else {
//         // --- Set derivatives
//         Idx<ddc::Deriv<Rg>> first_deriv_rg(1);
//         IdxStep<ddc::Deriv<Rg>> n_deriv_rg(1);
//         IdxRange<ddc::Deriv<Rg>> deriv_rg_idx_range(first_deriv_rg, n_deriv_rg);

//         IdxRange<ddc::Deriv<Rg>, GridThetag>
//                 derivs_rg_idx_range(deriv_rg_idx_range, TestFixture::idx_range_theta_g);

//         host_t<DFieldMem<IdxRange<ddc::Deriv<Rg>, GridThetag>>> derivs_rgmin_alloc(
//                 derivs_rg_idx_range);
//         host_t<DField<IdxRange<ddc::Deriv<Rg>, GridThetag>>> derivs_rgmin
//                 = get_field(derivs_rgmin_alloc);

//         host_t<DFieldMem<IdxRange<ddc::Deriv<Rg>, GridThetag>>> derivs_rgmax_alloc(
//                 derivs_rg_idx_range);
//         host_t<DField<IdxRange<ddc::Deriv<Rg>, GridThetag>>> derivs_rgmax
//                 = get_field(derivs_rgmax_alloc);

//         ddc::for_each(TestFixture::idx_range_theta_g, [&](Idx<GridThetag> const& idx_thetag) {
//             derivs_rgmin(first_deriv_rg, idx_thetag) = 3 * Kokkos::sin(ddc::coordinate(idx_thetag));
//             derivs_rgmax(first_deriv_rg, idx_thetag)
//                     = -3 * Kokkos::sin(ddc::coordinate(idx_thetag));
//         });

//         // --- Spline builder
//         builder_g(
//                 function_g_coef,
//                 get_const_field(function_g),
//                 std::optional(get_const_field(derivs_rgmin)),
//                 std::optional(get_const_field(derivs_rgmax)));
//     }

//     ddc::ConstantExtrapolationRule<Rg, Thetag> bc_rmin_g(TestFixture::rg_min);
//     ddc::ConstantExtrapolationRule<Rg, Thetag> bc_rmax_g(TestFixture::rg_max);
//     ddc::PeriodicExtrapolationRule<Thetag> bc_theta_g;
//     SplineRThetagEvaluator evaluator_g(bc_rmin_g, bc_rmax_g, bc_theta_g, bc_theta_g);


//     // Check the local derivatives with the global ones at the interfaces ========================
//     if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
//         // We test if the boundaries are well treated => only work with 5 cells to better identify an error.
//         // 5 cells -------------------------------------------------------------------------------
//         SingleInterfaceDerivativesCalculator<
//                 Interface_1_2,
//                 ddc::BoundCond::GREVILLE,
//                 ddc::BoundCond::GREVILLE> const
//                 derivatives_calculator(
//                         TestFixture::idx_range_rtheta1,
//                         TestFixture::idx_range_etaxi2);


//         ddc::for_each(TestFixture::idx_range_theta1, [&](Patch1::Idx2 const& idx2_1) {
//             Patch2::Idx2 idx2_2 = idx_convertor_12(idx2_1);

//             // Coordinate at the interface.
//             Coord<Rg, Thetag> interface_coord;
//             if constexpr (std::is_same_v<Edge1, EastEdge1>) {
//                 interface_coord = Coord<
//                         Rg,
//                         Thetag>(double(TestFixture::r1_max), double(ddc::coordinate(idx2_1)));
//             } else {
//                 interface_coord = Coord<Rg, Thetag>(
//                         double(TestFixture::r1_max),
//                         double(TestFixture::theta1_max - ddc::coordinate(idx2_1)));
//             }

//             // Coefficient c.
//             double const sum_values = derivatives_calculator.get_function_coefficients(
//                     get_const_field(function_2[idx2_2]),
//                     get_const_field(function_1[idx2_1]));

//             double const global_deriv
//                     = evaluator_g.deriv_dim_1(interface_coord, get_const_field(function_g_coef));

//             // Exact formula ---------------------------------------------------------------------
//             double const local_deriv = sum_values;
//             EXPECT_NEAR(local_deriv, global_deriv, 5e-13);
//         });
//     } else {
//         // 30 cells ------------------------------------------------------------------------------
//         TestFixture::check_exact_and_approximation(
//                 30,
//                 5e-13,
//                 function_1,
//                 function_2,
//                 evaluator_g,
//                 function_g_coef,
//                 idx_convertor_12);

//         // 20 cells ------------------------------------------------------------------------------
//         TestFixture::check_exact_and_approximation(
//                 20,
//                 1e-10,
//                 function_1,
//                 function_2,
//                 evaluator_g,
//                 function_g_coef,
//                 idx_convertor_12);

//         // 15 cells ------------------------------------------------------------------------------
//         TestFixture::check_exact_and_approximation(
//                 15,
//                 1e-7,
//                 function_1,
//                 function_2,
//                 evaluator_g,
//                 function_g_coef,
//                 idx_convertor_12);

//         // 10 cells ------------------------------------------------------------------------------
//         TestFixture::check_exact_and_approximation(
//                 10,
//                 5e-5,
//                 function_1,
//                 function_2,
//                 evaluator_g,
//                 function_g_coef,
//                 idx_convertor_12);

//         // 5 cells -------------------------------------------------------------------------------
//         TestFixture::check_exact_and_approximation(
//                 5,
//                 2e-2,
//                 function_1,
//                 function_2,
//                 evaluator_g,
//                 function_g_coef,
//                 idx_convertor_12);
//     }
// }