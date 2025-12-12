// SPDX-License-Identifier: MIT

#include <typeinfo>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <boost/type_index.hpp>
#include <gtest/gtest.h>

#include "3patches_2d_non_periodic_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "interface.hpp"
#include "interface_derivative_matrix.hpp"
#include "interface_derivatives_matrix_test_utils.hpp"
#include "interface_derivatives_test_utils.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "types.hpp"
#include "view.hpp"


/*
    Test InterfaceExactDerivativeMatrix on the following geometry:

    #if defined(REVERSE_PATCH1)
        |  1  |  2  |  3  | 
           ←↓   ↑→     ↑→

    #elif defined(REVERSE_PATCH2)
        |  1  |  2  |  3  | 
          ↑→     ←↓   ↑→

    #elif defined(REVERSE_PATCH3)
        |  1  |  2  |  3  | 
           ↑→    ↑→    ←↓  
    
    #elif defined(CHANGE_BOUND1)
        |  1  |  2  |  3  | 
           ←↑    ↑→    ↑→  
           
    #elif defined(CHANGE_BOUND3)
        |  1  |  2  |  3  | 
           ↑→    ↑→    ↓→  

        with the global X dimension with Hermite boundary conditions 
        and the global Y spline with additional points as closure condition 
        (ddc::BoundCond::GREVILLE).

    > test ddc::BoundCond::HERMITE boundary conditions. 
    > test application on the X direction. 
    > test application to compute first derivatives and cross-derivatives. 
    > test on a non-uniform patches. 
    > test with a middle ill-oriented interface. 
    > test with not matching directions of the patches.
    > test with exact formulation in SingleInterfaceDerivativeCalculator. 
    > test agreement between computed and global spline derivatives. 
    > test agreement between local and global splines.
*/

namespace {
// Multi-patch tags ---
using namespace non_periodic_non_uniform_2d_3patches;

// INTERFACES ------------------------------------------------------------------------------------
using NorthInterface1 = Interface<NorthEdge<1>, OutsideEdge, true>;
using NorthInterface2 = Interface<NorthEdge<2>, OutsideEdge, true>;
using NorthInterface3 = Interface<NorthEdge<3>, OutsideEdge, true>;

using SouthInterface1 = Interface<OutsideEdge, SouthEdge<1>, true>;
using SouthInterface2 = Interface<OutsideEdge, SouthEdge<2>, true>;
using SouthInterface3 = Interface<OutsideEdge, SouthEdge<3>, true>;

using EastInterface1 = Interface<OutsideEdge, EastEdge<1>, true>;
using EastInterface3 = Interface<OutsideEdge, EastEdge<3>, true>;

using WestInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
using WestInterface3 = Interface<OutsideEdge, WestEdge<3>, true>;

#if defined(REVERSE_PATCH1)
using Interface_1_2 = Interface<WestEdge<1>, WestEdge<2>, false>;
using Interface_2_3 = Interface<EastEdge<2>, WestEdge<3>, true>;

using OutsideInterface1 = Interface<OutsideEdge, EastEdge<1>, true>;
using OutsideInterface3 = Interface<EastEdge<3>, OutsideEdge, true>;

#elif defined(REVERSE_PATCH2)
using Interface_1_2 = Interface<EastEdge<1>, EastEdge<2>, false>;
using Interface_2_3 = Interface<WestEdge<2>, WestEdge<3>, false>;

using OutsideInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
using OutsideInterface3 = Interface<EastEdge<3>, OutsideEdge, true>;

#elif defined(REVERSE_PATCH3)
using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;
using Interface_2_3 = Interface<EastEdge<3>, EastEdge<2>, false>;

using OutsideInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
using OutsideInterface3 = Interface<WestEdge<3>, OutsideEdge, true>;

#elif defined(CHANGE_BOUND1)
using Interface_1_2 = Interface<SouthEdge<1>, WestEdge<2>, true>;
using Interface_2_3 = Interface<EastEdge<2>, WestEdge<3>, true>;

using OutsideInterface1 = Interface<OutsideEdge, NorthEdge<1>, true>;
using OutsideInterface3 = Interface<EastEdge<3>, OutsideEdge, true>;

#else
using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;
using Interface_2_3 = Interface<EastEdge<2>, SouthEdge<3>, false>;

using OutsideInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
using OutsideInterface3 = Interface<NorthEdge<3>, OutsideEdge, true>;

#endif


// CONNECTIVITY ----------------------------------------------------------------------------------
using Connectivity = MultipatchConnectivity<
#if defined(CHANGE_BOUND1)
        EastInterface1,
        WestInterface1,
#else
        NorthInterface1,
        SouthInterface1,
#endif
        NorthInterface2,
        SouthInterface2,
#if defined(CHANGE_BOUND3)
        EastInterface3,
        WestInterface3,
#else
        NorthInterface3,
        SouthInterface3,
#endif
        OutsideInterface1,
        OutsideInterface3,
        Interface_1_2,
        Interface_2_3>;

// Equivalent global mesh tags ---
struct Xg
{
    static bool constexpr PERIODIC = false;
};
struct Yg
{
    static bool constexpr PERIODIC = false;
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

using DerivXg = ddc::Deriv<Xg>;
using DerivYg = ddc::Deriv<Yg>;

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

// Interpolation points type for the patches.
template <std::size_t PatchIdx>
using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
        BSplinesX<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

template <std::size_t PatchIdx>
using SplineInterpPointsY = ddcHelper::NonUniformInterpolationPoints<
        BSplinesY<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

// Interpolation points type for the equivalent global spline.
using SplineInterpPointsXg = ddcHelper::
        NonUniformInterpolationPoints<BSplinesXg, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;
using SplineInterpPointsYg = ddcHelper::
        NonUniformInterpolationPoints<BSplinesYg, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;

// Operators on the equivalent global spline.
using SplineRThetagBuilder = ddc::SplineBuilder2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::SplineSolver::LAPACK>;

using SplineRThetagBuilderDerivField = SplineBuliderDerivField2D<
        HostExecSpace,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

using SplineRThetagEvaluator = ddc::SplineEvaluator2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::ConstantExtrapolationRule<Xg, Yg>,
        ddc::ConstantExtrapolationRule<Xg, Yg>,
        ddc::ConstantExtrapolationRule<Yg, Xg>,
        ddc::ConstantExtrapolationRule<Yg, Xg>>;

// Tools to get the equivalent global coordinate  of a local coordinate.
using CoordTransform1 = CoordTransform<Xg, Yg, X<1>, Y<1>>;
using CoordTransform2 = CoordTransform<Xg, Yg, X<2>, Y<2>>;
using CoordTransform3 = CoordTransform<Xg, Yg, X<3>, Y<3>>;

struct InterfaceExactDerivativeMatrixHermiteTest : public ::testing::Test
{
    // DEFINE BOUNDARIES OF THE DOMAINS ----------------------------------------------------------
    // patches 1 ---------------------------------
    static constexpr Coord<X<1>> x1_min = Coord<X<1>>(0.0);
    static constexpr Coord<X<1>> x1_max = Coord<X<1>>(1.0);
    static constexpr IdxStep<GridX<1>> x1_ncells = IdxStep<GridX<1>>(6);

    static constexpr Coord<Y<1>> y1_min = Coord<Y<1>>(0.0);
    static constexpr Coord<Y<1>> y1_max = Coord<Y<1>>(1.0);
    static constexpr IdxStep<GridY<1>> y1_ncells = IdxStep<GridY<1>>(6);

    // patches 2 ---------------------------------
    static constexpr Coord<X<2>> x2_min = Coord<X<2>>(1.0);
    static constexpr Coord<X<2>> x2_max = Coord<X<2>>(2.0);
    static constexpr IdxStep<GridX<2>> x2_ncells = IdxStep<GridX<2>>(6);

    static constexpr Coord<Y<2>> y2_min = Coord<Y<2>>(0.0);
    static constexpr Coord<Y<2>> y2_max = Coord<Y<2>>(1.0);
    static constexpr IdxStep<GridY<2>> y2_ncells = IdxStep<GridY<2>>(6);

// patches 3 ---------------------------------
#if defined(CHANGE_BOUND3)
    static constexpr Coord<X<3>> x3_min = Coord<X<3>>(0.0);
    static constexpr Coord<X<3>> x3_max = Coord<X<3>>(1.0);
    static constexpr IdxStep<GridX<3>> x3_ncells = IdxStep<GridX<3>>(6);

    static constexpr Coord<Y<3>> y3_min = Coord<Y<3>>(2.0);
    static constexpr Coord<Y<3>> y3_max = Coord<Y<3>>(3.0);
    static constexpr IdxStep<GridY<3>> y3_ncells = IdxStep<GridY<3>>(6);
#else
    static constexpr Coord<X<3>> x3_min = Coord<X<3>>(2.0);
    static constexpr Coord<X<3>> x3_max = Coord<X<3>>(3.0);
    static constexpr IdxStep<GridX<3>> x3_ncells = IdxStep<GridX<3>>(6);

    static constexpr Coord<Y<3>> y3_min = Coord<Y<3>>(0.0);
    static constexpr Coord<Y<3>> y3_max = Coord<Y<3>>(1.0);
    static constexpr IdxStep<GridY<3>> y3_ncells = IdxStep<GridY<3>>(6);
#endif

    // global ------------------------------------
    static constexpr Coord<Xg> xg_min = Coord<Xg> {double(x1_min)};
    static constexpr Coord<Xg> xg_max = Coord<Xg> {double(x3_max)};
    static constexpr IdxStep<GridXg> xg_ncells
            = IdxStep<GridXg>(x1_ncells.value() + x2_ncells.value() + x3_ncells.value());

    static constexpr Coord<Yg> yg_min = Coord<Yg> {double(y1_min)};
    static constexpr Coord<Yg> yg_max = Coord<Yg> {double(y1_max)};
    static constexpr IdxStep<GridYg> yg_ncells = IdxStep<GridYg>(y1_ncells.value());

    // coordinate transformation -----------------
    // --- for Patch1
#if defined(REVERSE_PATCH1)
    static constexpr bool is_reverse_x_1 = false;
    static constexpr bool is_reverse_y_1 = false;
    static constexpr bool are_exchange_x_y_1 = true;
#elif defined(CHANGE_BOUND1)
    static constexpr bool is_reverse_x_1 = true;
    static constexpr bool is_reverse_y_1 = false;
    static constexpr bool are_exchange_x_y_1 = false;
#else
    static constexpr bool is_reverse_x_1 = true;
    static constexpr bool is_reverse_y_1 = true;
    static constexpr bool are_exchange_x_y_1 = true;
#endif

    // --- for Patch2
#if defined(REVERSE_PATCH2)
    static constexpr bool is_reverse_x_2 = false;
    static constexpr bool is_reverse_y_2 = false;
    static constexpr bool are_exchange_x_y_2 = true;
#else
    static constexpr bool is_reverse_x_2 = true;
    static constexpr bool is_reverse_y_2 = true;
    static constexpr bool are_exchange_x_y_2 = true;
#endif

    // --- for Patch3
#if defined(REVERSE_PATCH3)
    static constexpr bool is_reverse_x_3 = false;
    static constexpr bool is_reverse_y_3 = false;
    static constexpr bool are_exchange_x_y_3 = true;
#elif defined(CHANGE_BOUND3)
    static constexpr bool is_reverse_x_3 = false;
    static constexpr bool is_reverse_y_3 = true;
    static constexpr bool are_exchange_x_y_3 = false;
#else
    static constexpr bool is_reverse_x_3 = true;
    static constexpr bool is_reverse_y_3 = true;
    static constexpr bool are_exchange_x_y_3 = true;
#endif

protected:
    const IdxRange<GridX<1>> idx_range_x1;
    const IdxRange<GridX<2>> idx_range_x2;
    const IdxRange<GridX<3>> idx_range_x3;

    const IdxRange<GridY<1>> idx_range_y1;
    const IdxRange<GridY<2>> idx_range_y2;
    const IdxRange<GridY<3>> idx_range_y3;

    const IdxRange<GridX<1>, GridY<1>> idx_range_xy1;
    const IdxRange<GridX<2>, GridY<2>> idx_range_xy2;
    const IdxRange<GridX<3>, GridY<3>> idx_range_xy3;

    const IdxRange<GridXg> idx_range_xg;
    const IdxRange<GridYg> idx_range_yg;
    const IdxRange<GridXg, GridYg> idx_range_xy_g;

    const CoordTransform1 coord_transform_1;
    const CoordTransform2 coord_transform_2;
    const CoordTransform3 coord_transform_3;

public:
    InterfaceExactDerivativeMatrixHermiteTest()
        : idx_range_x1(SplineInterpPointsX<1>::template get_domain<GridX<1>>())
        , idx_range_x2(SplineInterpPointsX<2>::template get_domain<GridX<2>>())
        , idx_range_x3(SplineInterpPointsX<3>::template get_domain<GridX<3>>())
        , idx_range_y1(SplineInterpPointsY<1>::template get_domain<GridY<1>>())
        , idx_range_y2(SplineInterpPointsY<2>::template get_domain<GridY<2>>())
        , idx_range_y3(SplineInterpPointsY<3>::template get_domain<GridY<3>>())
        , idx_range_xy1(idx_range_x1, idx_range_y1)
        , idx_range_xy2(idx_range_x2, idx_range_y2)
        , idx_range_xy3(idx_range_x3, idx_range_y3)
        , idx_range_xg(SplineInterpPointsXg::get_domain<GridXg>())
        , idx_range_yg(SplineInterpPointsYg::get_domain<GridYg>())
        , idx_range_xy_g(idx_range_xg, idx_range_yg)
        , coord_transform_1(
                  is_reverse_x_1,
                  is_reverse_y_1,
                  are_exchange_x_y_1,
                  x1_min,
                  x1_max,
                  y1_min,
                  y1_max)
        , coord_transform_2(
                  is_reverse_x_2,
                  is_reverse_y_2,
                  are_exchange_x_y_2,
                  x2_min,
                  x2_max,
                  y2_min,
                  y2_max)
        , coord_transform_3(
                  is_reverse_x_3,
                  is_reverse_y_3,
                  are_exchange_x_y_3,
                  x3_min,
                  x3_max,
                  y3_min,
                  y3_max)
    {
    }

    // INITIALISE DOMAINS ------------------------------------------------------------------------
    static void SetUpTestSuite()
    {
        // Creating of meshes and supports .......................................................
        // The patches are conforming between each other.
        std::vector<Coord<X<1>>> break_points_x1
                = build_random_non_uniform_break_points(x1_min, x1_max, x1_ncells);
        std::vector<Coord<X<2>>> break_points_x2
                = build_random_non_uniform_break_points(x2_min, x2_max, x2_ncells);
#if defined(CHANGE_BOUND3)
        std::vector<Coord<Y<3>>> break_points_y3
                = build_random_non_uniform_break_points(y3_min, y3_max, y3_ncells);
#else
        std::vector<Coord<X<3>>> break_points_x3
                = build_random_non_uniform_break_points(x3_min, x3_max, x3_ncells);
#endif

        std::vector<Coord<Y<1>>> break_points_y
                = build_random_non_uniform_break_points(y1_min, y1_max, y1_ncells);

#if defined(REVERSE_PATCH1)
        std::vector<Coord<Y<1>>> break_points_y1;
        fill_in_reverse(break_points_y1, break_points_y);
        std::vector<Coord<Y<2>>> break_points_y2 = convert_dim<Y<2>, Y<1>>(break_points_y);
        std::vector<Coord<Y<3>>> break_points_y3 = convert_dim<Y<3>, Y<1>>(break_points_y);
#elif defined(REVERSE_PATCH2)
        std::vector<Coord<Y<1>>> break_points_y1 = break_points_y;
        std::vector<Coord<Y<2>>> break_points_y2;
        fill_in_reverse(break_points_y2, break_points_y);
        std::vector<Coord<Y<3>>> break_points_y3 = convert_dim<Y<3>, Y<1>>(break_points_y);
#elif defined(REVERSE_PATCH3)
        std::vector<Coord<Y<1>>> break_points_y1 = break_points_y;
        std::vector<Coord<Y<2>>> break_points_y2 = convert_dim<Y<2>, Y<1>>(break_points_y);
        std::vector<Coord<Y<3>>> break_points_y3;
        fill_in_reverse(break_points_y3, break_points_y);
#elif (CHANGE_BOUND1)
        std::vector<Coord<Y<1>>> break_points_y1 = convert_dim<Y<1>, X<1>>(break_points_x1);
        break_points_x1 = convert_dim<X<1>, Y<1>>(break_points_y);
        std::vector<Coord<Y<2>>> break_points_y2 = convert_dim<Y<2>, Y<1>>(break_points_y);
        std::vector<Coord<Y<3>>> break_points_y3 = convert_dim<Y<3>, Y<1>>(break_points_y);
#else
        std::vector<Coord<Y<1>>> break_points_y1 = break_points_y;
        std::vector<Coord<Y<2>>> break_points_y2 = convert_dim<Y<2>, Y<1>>(break_points_y);
        std::vector<Coord<X<3>>> break_points_x3;
        fill_in_reverse(break_points_x3, break_points_y);
#endif

        std::vector<Coord<X<1>>> interpolation_points_x1 = break_points_x1;
        std::vector<Coord<X<2>>> interpolation_points_x2 = break_points_x2;
        std::vector<Coord<X<3>>> interpolation_points_x3 = break_points_x3;

        std::vector<Coord<Y<1>>> interpolation_points_y1 = break_points_y1;
        std::vector<Coord<Y<2>>> interpolation_points_y2 = break_points_y2;
        std::vector<Coord<Y<3>>> interpolation_points_y3 = break_points_y3;

        // Patch 1 ...............................................................................
        ddc::init_discrete_space<BSplinesX<1>>(break_points_x1);
        ddc::init_discrete_space<BSplinesY<1>>(break_points_y1);

        ddc::init_discrete_space<GridX<1>>(interpolation_points_x1);
        ddc::init_discrete_space<GridY<1>>(interpolation_points_y1);

        // Patch 2 ...............................................................................
        ddc::init_discrete_space<BSplinesX<2>>(break_points_x2);
        ddc::init_discrete_space<BSplinesY<2>>(break_points_y2);

        ddc::init_discrete_space<GridX<2>>(interpolation_points_x2);
        ddc::init_discrete_space<GridY<2>>(interpolation_points_y2);

        // Patch 3 ...............................................................................
        ddc::init_discrete_space<BSplinesX<3>>(break_points_x3);
        ddc::init_discrete_space<BSplinesY<3>>(break_points_y3);

        ddc::init_discrete_space<GridX<3>>(interpolation_points_x3);
        ddc::init_discrete_space<GridY<3>>(interpolation_points_y3);

        // Equivalent global domain ..............................................................
        std::vector<Coord<Xg>> break_points_xg;
        std::vector<Coord<Xg>> interpolation_points_xg;

        // --- break points from Patch1
#if defined(REVERSE_PATCH1)
        std::vector<Coord<X<1>>> break_points_x1_reverse;
        std::vector<Coord<X<1>>> interpolation_points_x1_reverse;
        fill_in_reverse(break_points_x1_reverse, break_points_x1);
        fill_in_reverse(interpolation_points_x1_reverse, interpolation_points_x1);
        break_points_x1_reverse.pop_back();
        interpolation_points_x1_reverse.pop_back();
        fill_in(break_points_xg, break_points_x1_reverse);
        fill_in(interpolation_points_xg, interpolation_points_x1_reverse);
#elif (CHANGE_BOUND1)
        std::vector<Coord<X<1>>> break_points_y1_reverse;
        std::vector<Coord<X<1>>> interpolation_points_y1_reverse;
        fill_in_reverse(break_points_y1_reverse, break_points_y1);
        fill_in_reverse(interpolation_points_y1_reverse, interpolation_points_y1);
        break_points_y1_reverse.pop_back();
        interpolation_points_y1_reverse.pop_back();
        fill_in(break_points_xg, break_points_y1_reverse);
        fill_in(interpolation_points_xg, interpolation_points_y1_reverse);
#else
        break_points_x1.pop_back();
        interpolation_points_x1.pop_back();
        fill_in(break_points_xg, break_points_x1);
        fill_in(interpolation_points_xg, interpolation_points_x1);
#endif

        // --- break points from Patch2
#if defined(REVERSE_PATCH2)
        std::vector<Coord<X<2>>> break_points_x2_reverse;
        std::vector<Coord<X<2>>> interpolation_points_x2_reverse;
        fill_in_reverse(break_points_x2_reverse, break_points_x2);
        fill_in_reverse(interpolation_points_x2_reverse, interpolation_points_x2);
        break_points_x2_reverse.pop_back();
        interpolation_points_x2_reverse.pop_back();
        fill_in(break_points_xg, break_points_x2_reverse);
        fill_in(interpolation_points_xg, interpolation_points_x2_reverse);
#else
        break_points_x2.pop_back();
        interpolation_points_x2.pop_back();
        fill_in(break_points_xg, break_points_x2);
        fill_in(interpolation_points_xg, interpolation_points_x2);
#endif

        // --- break points from Patch3
#if defined(REVERSE_PATCH3)
        fill_in_reverse(break_points_xg, break_points_x3);
        fill_in_reverse(interpolation_points_xg, interpolation_points_x3);
#elif (CHANGE_BOUND3)
        fill_in(break_points_xg, break_points_y3);
        fill_in(interpolation_points_xg, interpolation_points_y3);
#else
        fill_in(break_points_xg, break_points_x3);
        fill_in(interpolation_points_xg, interpolation_points_x3);
#endif

        std::vector<Coord<Yg>> break_points_yg;
        std::vector<Coord<Yg>> interpolation_points_yg;

        fill_in(break_points_yg, break_points_y);
        fill_in(interpolation_points_yg, break_points_y);

        ddc::init_discrete_space<BSplinesXg>(break_points_xg);
        ddc::init_discrete_space<BSplinesYg>(break_points_yg);

        ddc::init_discrete_space<GridXg>(interpolation_points_xg);
        ddc::init_discrete_space<GridYg>(interpolation_points_yg);
    }
};
} // end namespace



// Check that the local grids and the equivalent global grid match together.
TEST_F(InterfaceExactDerivativeMatrixHermiteTest, InterpolationPointsCheck)
{
    int const x_shift1 = x1_ncells.value();
    int const x_shift2 = x1_ncells.value() + x2_ncells.value();

    check_interpolation_grids<Patch1, GridXg, GridYg>(idx_range_xy1, 0, 0, coord_transform_1);
    check_interpolation_grids<
            Patch2,
            GridXg,
            GridYg>(idx_range_xy2, x_shift1, 0, coord_transform_2);
    check_interpolation_grids<
            Patch3,
            GridXg,
            GridYg>(idx_range_xy3, x_shift2, 0, coord_transform_3);
}



TEST_F(InterfaceExactDerivativeMatrixHermiteTest, CheckForHermiteBc)
{
    std::tuple coord_transforms(coord_transform_1, coord_transform_2, coord_transform_3);

    // Instantiate the derivatives calculators ---------------------------------------------------
    // SingleInterfaceDerivativesCalculators for interfaces along y (periodic).
    SingleInterfaceDerivativesCalculator<Interface_1_2> const
            derivatives_calculator_1_2(idx_range_xy1, idx_range_xy2);
    SingleInterfaceDerivativesCalculator<Interface_2_3> const
            derivatives_calculator_2_3(idx_range_xy2, idx_range_xy3);

    // Collect the derivative calculators --------------------------------------------------------
    // We do not follow the physical order to test the operator.
    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect(derivatives_calculator_2_3, derivatives_calculator_1_2);

    // Collect the index ranges ------------------------------------------------------------------
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2, Patch3>
            idx_ranges(idx_range_xy1, idx_range_xy2, idx_range_xy3);

    // Instantiate the matrix calculators --------------------------------------------------------

    InterfaceExactDerivativeMatrix<
            Connectivity,
#if defined(CHANGE_BOUND1)
            GridY<1>,
#else
            GridX<1>,
#endif
            ddc::detail::TypeSeq<Patch1, Patch2, Patch3>,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            SingleInterfaceDerivativesCalculatorCollection<Interface_2_3, Interface_1_2>>
            matrix(idx_ranges, deriv_calculators_collect);

    // Instantiate DerivField ====================================================================
    // Instantiate index range slices ------------------------------------------------------------
    IdxRangeSlice<GridX<1>> idx_range_slice_dx1 = get_bound_idx_range_slice(idx_range_x1);
    IdxRangeSlice<GridX<2>> idx_range_slice_dx2 = get_bound_idx_range_slice(idx_range_x2);
    IdxRangeSlice<GridX<3>> idx_range_slice_dx3 = get_bound_idx_range_slice(idx_range_x3);

    IdxRangeSlice<GridY<1>> idx_range_slice_dy1 = get_bound_idx_range_slice(idx_range_y1);
    IdxRangeSlice<GridY<2>> idx_range_slice_dy2 = get_bound_idx_range_slice(idx_range_y2);
    IdxRangeSlice<GridY<3>> idx_range_slice_dy3 = get_bound_idx_range_slice(idx_range_y3);

    // Collect the index range slices.
    MultipatchType<IdxRange1SliceOnPatch, Patch1, Patch2, Patch3>
            idx_ranges_slice_dx(idx_range_slice_dx1, idx_range_slice_dx2, idx_range_slice_dx3);

    MultipatchType<IdxRange2SliceOnPatch, Patch1, Patch2, Patch3>
            idx_ranges_slice_dy(idx_range_slice_dy1, idx_range_slice_dy2, idx_range_slice_dy3);

    // Instantiate DerivField --------------------------------------------------------------------
    DerivFieldMemOnPatch_host<Patch1>
            function_and_derivs_1_alloc(idx_range_xy1, idx_range_slice_dx1, idx_range_slice_dy1);
    DerivFieldMemOnPatch_host<Patch2>
            function_and_derivs_2_alloc(idx_range_xy2, idx_range_slice_dx2, idx_range_slice_dy2);
    DerivFieldMemOnPatch_host<Patch3>
            function_and_derivs_3_alloc(idx_range_xy3, idx_range_slice_dx3, idx_range_slice_dy3);

    DerivFieldOnPatch_host<Patch1> function_and_derivs_1(function_and_derivs_1_alloc);
    DerivFieldOnPatch_host<Patch2> function_and_derivs_2(function_and_derivs_2_alloc);
    DerivFieldOnPatch_host<Patch3> function_and_derivs_3(function_and_derivs_3_alloc);

    // Collect the fields with derivatives.
    MultipatchField<DerivFieldOnPatch_host, Patch1, Patch2, Patch3> functions_and_derivs(
            function_and_derivs_1,
            function_and_derivs_2,
            function_and_derivs_3);

    // Instantiate the global function.
    IdxRangeSlice<GridXg> idx_range_slice_dxg = get_bound_idx_range_slice(idx_range_xg);
    IdxRangeSlice<GridYg> idx_range_slice_dyg = get_bound_idx_range_slice(idx_range_yg);

    DerivFieldMem<double, IdxRange<DerivXg, GridXg, DerivYg, GridYg>, 1>
            function_and_derivs_g_alloc(idx_range_xy_g, idx_range_slice_dxg, idx_range_slice_dyg);
    DerivField<double, IdxRange<DerivXg, GridXg, DerivYg, GridYg>> function_and_derivs_g(
            function_and_derivs_g_alloc);

    // Initialise the data =======================================================================
    // --- the function values.
    initialise_all_functions<Xg, Yg>(functions_and_derivs, coord_transforms);
    initialise_2D_function<GridXg, GridYg, CoordTransform<Xg, Yg, Xg, Yg>>(
            function_and_derivs_g.get_values_field());

    // --- the derivatives of the equivalent global spline.
    Idx<DerivXg> first_dxg(1);
    Idx<DerivYg> first_dyg(1);

    Idx<GridXg> idx_slice_xg_min(idx_range_slice_dxg.front());
    Idx<GridXg> idx_slice_xg_max(idx_range_slice_dxg.back());
    Idx<GridYg> idx_slice_yg_min(idx_range_slice_dyg.front());
    Idx<GridYg> idx_slice_yg_max(idx_range_slice_dyg.back());

    Idx<DerivXg, GridXg> idx_dxg_min(first_dxg, idx_slice_xg_min);
    Idx<DerivXg, GridXg> idx_dxg_max(first_dxg, idx_slice_xg_max);
    Idx<DerivYg, GridYg> idx_dyg_min(first_dyg, idx_slice_yg_min);
    Idx<DerivYg, GridYg> idx_dyg_max(first_dyg, idx_slice_yg_max);

    Idx<DerivXg, GridXg, DerivYg, GridYg> idx_dxgdyg_min_min(idx_dxg_min, idx_dyg_min);
    Idx<DerivXg, GridXg, DerivYg, GridYg> idx_dxgdyg_max_min(idx_dxg_max, idx_dyg_min);
    Idx<DerivXg, GridXg, DerivYg, GridYg> idx_dxgdyg_min_max(idx_dxg_min, idx_dyg_max);
    Idx<DerivXg, GridXg, DerivYg, GridYg> idx_dxgdyg_max_max(idx_dxg_max, idx_dyg_max);

    ddc::for_each(idx_range_yg, [&](Idx<GridYg> const idx) {
        double const xgmin = xg_min;
        double const xgmax = xg_max;
        double const yg = ddc::coordinate(idx);
        function_and_derivs_g[idx_dxg_min](idx)
                = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xgmin + 0.25) * std::sin(yg);
        function_and_derivs_g[idx_dxg_max](idx)
                = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xgmax + 0.25) * std::sin(yg);
    });
    ddc::for_each(idx_range_xg, [&](Idx<GridXg> const idx) {
        double const ygmin = yg_min;
        double const ygmax = yg_max;
        double const xg = ddc::coordinate(Idx<GridXg>(idx));
        function_and_derivs_g[idx_dyg_min](idx)
                = std::cos(2. / 3 * M_PI * xg + 0.25) * std ::cos(ygmin);
        function_and_derivs_g[idx_dyg_max](idx)
                = std::cos(2. / 3 * M_PI * xg + 0.25) * std ::cos(ygmax);
    });
    function_and_derivs_g(idx_dxgdyg_min_min)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xg_min + 0.25) * std::sin(yg_min);
    function_and_derivs_g(idx_dxgdyg_max_min)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xg_max + 0.25) * std::sin(yg_min);
    function_and_derivs_g(idx_dxgdyg_min_max)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xg_min + 0.25) * std::sin(yg_max);
    function_and_derivs_g(idx_dxgdyg_max_max)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xg_max + 0.25) * std::sin(yg_max);

    // --- the local derivatives from an equivalent global spline.
    // ------- build global spline representation
    SplineRThetagBuilder builder_g(idx_range_xy_g);
    SplineRThetagBuilderDerivField apply_builder_g(builder_g);

    host_t<DFieldMem<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef_alloc(
            builder_g.batched_spline_domain(idx_range_xy_g));
    host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef
            = get_field(function_g_coef_alloc);

    apply_builder_g(function_g_coef, function_and_derivs_g);

    host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const_function_g_coef
            = get_const_field(function_g_coef);

    // ------ global spline evaluator
    ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymin_g(yg_min, xg_min, xg_max);
    ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymax_g(yg_max, xg_min, xg_max);
    ddc::ConstantExtrapolationRule<Xg, Yg> bc_xmin_g(xg_min, yg_min, yg_max);
    ddc::ConstantExtrapolationRule<Xg, Yg> bc_xmax_g(xg_max, yg_min, yg_max);
    SplineRThetagEvaluator evaluator_g(bc_xmin_g, bc_xmax_g, bc_ymin_g, bc_ymax_g);

    // ------ initialise the boundary first derivatives from the global spline
    // X bound on Patch1 ---
#if defined(REVERSE_PATCH1)
    Idx<ddc::Deriv<X<1>>, GridX<1>>
            idx_slice_xmax_1(Idx<ddc::Deriv<X<1>>>(1), idx_range_slice_dx1.back());

    DField<IdxRange<typename Patch1::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_x1max_extracted = function_and_derivs_1[idx_slice_xmax_1];
    ddc::for_each(idx_range_y1, [&](Idx<GridY<1>> const& idx_par) {
        Idx<GridX<1>, GridY<1>> idx_max(idx_range_x1.back(), idx_par);
        Coord<Xg, Yg> interface_coord_max(
                coord_transform_1.get_global_coord(ddc::coordinate(idx_max)));
        derivs_x1max_extracted(idx_par)
                = -evaluator_g.deriv_dim_1(interface_coord_max, const_function_g_coef);
    });
#elif defined(CHANGE_BOUND1)
    Idx<ddc::Deriv<Y<1>>, GridY<1>>
            idx_slice_ymax_1(Idx<ddc::Deriv<Y<1>>>(1), idx_range_slice_dy1.back());

    DField<IdxRange<typename Patch1::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_y1max_extracted = function_and_derivs_1[idx_slice_ymax_1];
    ddc::for_each(idx_range_x1, [&](Idx<GridX<1>> const& idx_par) {
        double y = ddc::coordinate(idx_par);
        double x = y1_min;
        Coord<Xg, Yg> interface_coord_max(x, y);
        derivs_y1max_extracted(idx_par)
                = -evaluator_g.deriv_dim_1(interface_coord_max, const_function_g_coef);
    });
#else
    Idx<ddc::Deriv<X<1>>, GridX<1>>
            idx_slice_xmin_1(Idx<ddc::Deriv<X<1>>>(1), idx_range_slice_dx1.front());

    DField<IdxRange<typename Patch1::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_x1min_extracted = function_and_derivs_1[idx_slice_xmin_1];
    ddc::for_each(idx_range_y1, [&](Idx<GridY<1>> const& idx_par) {
        Idx<GridX<1>, GridY<1>> idx_min(idx_range_x1.front(), idx_par);
        Coord<Xg, Yg> interface_coord_min(
                coord_transform_1.get_global_coord(ddc::coordinate(idx_min)));
        derivs_x1min_extracted(idx_par)
                = evaluator_g.deriv_dim_1(interface_coord_min, const_function_g_coef);
    });
#endif

    // X bound on Patch3 ---
#if defined(REVERSE_PATCH3)
    Idx<ddc::Deriv<X<3>>, GridX<3>>
            idx_slice_xmin_3(Idx<ddc::Deriv<X<3>>>(1), idx_range_slice_dx3.front());

    DField<IdxRange<typename Patch3::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_x3min_extracted = function_and_derivs_3[idx_slice_xmin_3];
    ddc::for_each(idx_range_y3, [&](Idx<GridY<3>> const& idx_par) {
        Idx<GridX<3>, GridY<3>> idx_min(idx_range_x3.front(), idx_par);
        Coord<Xg, Yg> interface_coord_min(
                coord_transform_3.get_global_coord(ddc::coordinate(idx_min)));
        derivs_x3min_extracted(idx_par)
                = -evaluator_g.deriv_dim_1(interface_coord_min, const_function_g_coef);
    });
#elif (CHANGE_BOUND3)
    Idx<ddc::Deriv<Y<3>>, GridY<3>>
            idx_slice_ymax_3(Idx<ddc::Deriv<Y<3>>>(1), idx_range_slice_dy3.back());

    DField<IdxRange<typename Patch3::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_y3max_extracted = function_and_derivs_3[idx_slice_ymax_3];
    ddc::for_each(idx_range_x3, [&](Idx<GridX<3>> const& idx_par) {
        double y = x3_min + x3_max - ddc::coordinate(idx_par);
        double x = ddc::coordinate(idx_range_y3.back());
        Coord<Xg, Yg> interface_coord_max(x, y);
        derivs_y3max_extracted(idx_par)
                = evaluator_g.deriv_dim_1(interface_coord_max, const_function_g_coef);
    });
#else
    Idx<ddc::Deriv<X<3>>, GridX<3>>
            idx_slice_xmax_3(Idx<ddc::Deriv<X<3>>>(1), idx_range_slice_dx3.back());

    DField<IdxRange<typename Patch3::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_x3max_extracted = function_and_derivs_3[idx_slice_xmax_3];
    ddc::for_each(idx_range_y3, [&](Idx<GridY<3>> const& idx_par) {
        Idx<GridX<3>, GridY<3>> idx_max(idx_range_x3.back(), idx_par);
        Coord<Xg, Yg> interface_coord_max(
                coord_transform_3.get_global_coord(ddc::coordinate(idx_max)));
        derivs_x3max_extracted(idx_par)
                = evaluator_g.deriv_dim_1(interface_coord_max, const_function_g_coef);
    });
#endif

    // Y bound on Patch1 ---
#if (CHANGE_BOUND1)
    initialise_x_derivatives<Patch1>(
            function_and_derivs_1,
            idx_range_slice_dx1,
            evaluator_g,
            const_function_g_coef,
            coord_transform_1);
#else
    initialise_y_derivatives<Patch1>(
            function_and_derivs_1,
            idx_range_slice_dy1,
            evaluator_g,
            const_function_g_coef,
            coord_transform_1);
#endif

    // Y bound on Patch2 ---
    initialise_y_derivatives<Patch2>(
            function_and_derivs_2,
            idx_range_slice_dy2,
            evaluator_g,
            const_function_g_coef,
            coord_transform_2);

    // Y bound on Patch3 ---
#if (CHANGE_BOUND3)
    initialise_x_derivatives<Patch3>(
            function_and_derivs_3,
            idx_range_slice_dx3,
            evaluator_g,
            const_function_g_coef,
            coord_transform_3);
#else
    initialise_y_derivatives<Patch3>(
            function_and_derivs_3,
            idx_range_slice_dy3,
            evaluator_g,
            const_function_g_coef,
            coord_transform_3);
#endif

    // ------ initialise the cross-derivatives from the global spline
    initialise_all_cross_derivatives(
            functions_and_derivs,
            idx_ranges_slice_dx,
            idx_ranges_slice_dy,
            evaluator_g,
            const_function_g_coef,
            coord_transforms);

    // --- the first derivatives (on inner interfaces) from the function values.
    matrix.solve_deriv(functions_and_derivs);

    // --- the cross-derivatives from the first derivatives.
    /*
        Here, it is not needed to compute the cross-derivatives because
        they are given by the boundary conditions. However, we want to 
        check that the matrix computes correctly the values. 
    */
    matrix.solve_cross_deriv(functions_and_derivs);

    // Test the values of the derivatives ========================================================
    using EmptyPatchSeq = ddc::detail::TypeSeq<>;

    // Check each derivatives ---
    check_all_x_derivatives(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            idx_ranges,
            idx_ranges_slice_dx,
            coord_transforms);
    check_all_y_derivatives<EmptyPatchSeq, EmptyPatchSeq>(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            idx_ranges,
            idx_ranges_slice_dy,
            coord_transforms);
    check_all_xy_derivatives<EmptyPatchSeq, EmptyPatchSeq>(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            idx_ranges,
            idx_ranges_slice_dx,
            idx_ranges_slice_dy,
            coord_transforms);

    // Check the whole spline representations ---
    check_all_spline_representation_agreement<EmptyPatchSeq, EmptyPatchSeq>(
            idx_ranges,
            idx_ranges_slice_dx,
            idx_ranges_slice_dy,
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            coord_transforms);
}
