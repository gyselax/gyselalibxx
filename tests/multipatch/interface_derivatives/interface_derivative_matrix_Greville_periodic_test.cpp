// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "9patches_2d_periodic_strips_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "interface.hpp"
#include "interface_derivative_matrix.hpp"
#include "interface_derivatives_matrix_test_utils.hpp"
#include "interface_derivatives_test_utils.hpp"
#include "linear_coord_transform.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "orthogonal_coord_transforms.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "types.hpp"
#include "view.hpp"


/*
    Test InterfaceDerivativeMatrix on the following geometry:

        |  1  |  2  |  3  |  1 ...
        -------------------
        |  4  |  5  |  6  |  4 ...
        -------------------
        |  7  |  8  |  9  |  7 ... 

        with the global X dimension periodic and the global Y spline
        with additional points as closure condition (ddc::BoundCond::GREVILLE).

    > test ddc::BoundCond::PERIODIC boundary conditions. 
    > test ddc::BoundCond::GREVILLE boundary conditions. 
    > test application on the X and Y directions. 
    > test application to compute first derivatives and cross-derivatives. 
    > test on a non-uniform patches. 
    > test with exact formulation in SingleInterfaceDerivativeCalculator. 
    > test agreement between computed and global spline derivatives. 
    > test agreement between local and global splines.
*/

namespace {
// Multi-patch tags ---
using namespace periodic_strips_non_uniform_2d_9patches;

// Equivalent global mesh tags ---
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

// Interpolation points type for the patches.
template <std::size_t PatchIdx>
using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
        BSplinesX<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

template <std::size_t PatchIdx, ddc::BoundCond BoundCondMin, ddc::BoundCond BoundCondMax>
using SplineInterpPointsY
        = ddcHelper::NonUniformInterpolationPoints<BSplinesY<PatchIdx>, BoundCondMin, BoundCondMax>;

// Interpolation points type for the equivalent global spline.
using SplineInterpPointsXg = ddcHelper::NonUniformInterpolationPoints<
        BSplinesXg,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;
using SplineInterpPointsYg = ddcHelper::NonUniformInterpolationPoints<
        BSplinesYg,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE>;

// Operators on the equivalent global spline.
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

template<std::size_t I>
struct CoordTransformGroup {
    using XTransform = LinearCoordTransform<Xg, X<I>>;
    using YTransform = LinearCoordTransform<Yg, Y<I>>;
    XTransform x_transform;
    YTransform y_transform;
    OrthogonalCoordTransforms<
            Coord<Yg, Xg>,
            Coord<X<I>, Y<I>>,
            Coord<Xg, Yg>,
            XTransform,
            YTransform>
            coord_transform;
    CoordTransformGroup(Coord<Xg> xg_min, Coord<Yg> yg_min, Coord<X<I>> xl_min, Coord<Y<I>> yl_min)
        : x_transform(xg_min, xl_min, 1.0)
        , y_transform(yg_min, yl_min, 1.0)
        , coord_transform(x_transform, y_transform) {}
};

struct InterfaceDerivativeMatrixGrevillePeriodicTest : public ::testing::Test
{
    // DEFINE BOUNDARIES OF THE DOMAINS ----------------------------------------------------------
    // We define only the patch of the first row or column as the patches are sharing data.
    // patches 1 | 4 | 7  dim X ------------------
    static constexpr Coord<X<1>> x1_min = Coord<X<1>>(0.0);
    static constexpr Coord<X<1>> x1_max = Coord<X<1>>(1.0);
    static constexpr IdxStep<GridX<1>> x1_ncells = IdxStep<GridX<1>>(10);

    // patches 2 | 5 | 8  dim X ------------------
    static constexpr Coord<X<2>> x2_min = Coord<X<2>>(1.0);
    static constexpr Coord<X<2>> x2_max = Coord<X<2>>(2.0);
    static constexpr IdxStep<GridX<2>> x2_ncells = IdxStep<GridX<2>>(10);

    // patches 3 | 6 | 9  dim X ------------------
    static constexpr Coord<X<3>> x3_min = Coord<X<3>>(2.0);
    static constexpr Coord<X<3>> x3_max = Coord<X<3>>(3.0);
    static constexpr IdxStep<GridX<3>> x3_ncells = IdxStep<GridX<3>>(10);

    // patches 1 | 2 | 3  dim Y --------------------
    static constexpr Coord<Y<1>> y1_min = Coord<Y<1>>(2.0);
    static constexpr Coord<Y<1>> y1_max = Coord<Y<1>>(3.0);
    static constexpr IdxStep<GridY<1>> y1_ncells = IdxStep<GridY<1>>(10);

    // patches 4 | 5 | 6  dim Y --------------------
    static constexpr Coord<Y<4>> y4_min = Coord<Y<4>>(1.0);
    static constexpr Coord<Y<4>> y4_max = Coord<Y<4>>(2.0);
    static constexpr IdxStep<GridY<4>> y4_ncells = IdxStep<GridY<4>>(10);

    // patches 7 | 8 | 9  dim Y --------------------
    static constexpr Coord<Y<7>> y7_min = Coord<Y<7>>(0.0);
    static constexpr Coord<Y<7>> y7_max = Coord<Y<7>>(1.0);
    static constexpr IdxStep<GridY<7>> y7_ncells = IdxStep<GridY<7>>(10);

    // global ------------------------------------
    static constexpr Coord<Yg> yg_min = convert_dim<Yg, Y<7>>(y7_min);
    static constexpr Coord<Yg> yg_max = convert_dim<Yg, Y<1>>(y1_max);

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

    const IdxRange<GridXg> idx_range_xg;
    const IdxRange<GridYg> idx_range_yg;
    const IdxRange<GridXg, GridYg> idx_range_xy_g;

public:
    InterfaceDerivativeMatrixGrevillePeriodicTest()
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
        , idx_range_xg(SplineInterpPointsXg::get_domain<GridXg>())
        , idx_range_yg(SplineInterpPointsYg::get_domain<GridYg>())
        , idx_range_xy_g(idx_range_xg, idx_range_yg)
    {
    }

    // INITIALISE DOMAINS ------------------------------------------------------------------------
    static void SetUpTestSuite()
    {
        // Creating of meshes and supports .......................................................
        // The patches are conforming between each other.
        std::vector<Coord<X<1>>> break_points_x147
                = build_random_non_uniform_break_points(x1_min, x1_max, x1_ncells);
        std::vector<Coord<X<2>>> break_points_x258
                = build_random_non_uniform_break_points(x2_min, x2_max, x2_ncells);
        std::vector<Coord<X<3>>> break_points_x369
                = build_random_non_uniform_break_points(x3_min, x3_max, x3_ncells);

        std::vector<Coord<Y<1>>> break_points_y123
                = build_random_non_uniform_break_points(y1_min, y1_max, y1_ncells);
        std::vector<Coord<Y<4>>> break_points_y456
                = build_random_non_uniform_break_points(y4_min, y4_max, y4_ncells);
        std::vector<Coord<Y<7>>> break_points_y789
                = build_random_non_uniform_break_points(y7_min, y7_max, y7_ncells);

        std::vector<Coord<X<1>>> interpolation_points_x147 = break_points_x147;
        std::vector<Coord<X<2>>> interpolation_points_x258 = break_points_x258;
        std::vector<Coord<X<3>>> interpolation_points_x369 = break_points_x369;

        std::vector<Coord<Y<1>>> interpolation_points_y123
                = get_interpolation_points_add_one_on_right(break_points_y123);
        std::vector<Coord<Y<4>>> interpolation_points_y456 = break_points_y456;
        std::vector<Coord<Y<7>>> interpolation_points_y789
                = get_interpolation_points_add_one_on_left(break_points_y789);

        // Patch 1 ...............................................................................
        init_space<1>(break_points_x147, break_points_y123, interpolation_points_y123);
        init_space<2>(break_points_x258, break_points_y123, interpolation_points_y123);
        init_space<3>(break_points_x369, break_points_y123, interpolation_points_y123);
        init_space<4>(break_points_x147, break_points_y456, interpolation_points_y456);
        init_space<5>(break_points_x258, break_points_y456, interpolation_points_y456);
        init_space<6>(break_points_x369, break_points_y456, interpolation_points_y456);
        init_space<7>(break_points_x147, break_points_y789, interpolation_points_y789);
        init_space<8>(break_points_x258, break_points_y789, interpolation_points_y789);
        init_space<9>(break_points_x369, break_points_y789, interpolation_points_y789);

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

    template <int I, int XI, int YI>
    static void init_space(
            std::vector<Coord<X<XI>>> break_points_x,
            std::vector<Coord<Y<YI>>> break_points_y,
            std::vector<Coord<Y<YI>>> interpolation_points_y)
    {
        static_assert((I - 1) % 3 + 1 == XI);
        static_assert(I - (I - 1) % 3 == YI);
        ddc::init_discrete_space<BSplinesX<I>>(convert_dim<X<I>, X<XI>>(break_points_x));
        ddc::init_discrete_space<BSplinesY<I>>(convert_dim<Y<I>, Y<YI>>(break_points_y));

        ddc::init_discrete_space<GridX<I>>(convert_dim<X<I>, X<XI>>(break_points_x));
        ddc::init_discrete_space<GridY<I>>(convert_dim<Y<I>, Y<YI>>(interpolation_points_y));
    }
};

} // end namespace



TEST_F(InterfaceDerivativeMatrixGrevillePeriodicTest, CheckForPeriodicAndGrevilleBC)
{
    Coord<Xg> xg_min = ddc::coordinate(idx_range_xg.front());
    std::tuple coord_transforms {
        CoordTransformGroup<1>(xg_min, yg_min, x1_min, convert_dim<Y<1>>(y7_min)),
        CoordTransformGroup<2>(xg_min, yg_min, convert_dim<X<2>>(x1_min), convert_dim<Y<2>>(y7_min)),
        CoordTransformGroup<3>(xg_min, yg_min, convert_dim<X<3>>(x1_min), convert_dim<Y<3>>(y7_min)),
        CoordTransformGroup<4>(xg_min, yg_min, convert_dim<X<4>>(x1_min), convert_dim<Y<4>>(y7_min)),
        CoordTransformGroup<5>(xg_min, yg_min, convert_dim<X<5>>(x1_min), convert_dim<Y<5>>(y7_min)),
        CoordTransformGroup<6>(xg_min, yg_min, convert_dim<X<6>>(x1_min), convert_dim<Y<6>>(y7_min)),
        CoordTransformGroup<7>(xg_min, yg_min, convert_dim<X<7>>(x1_min), y7_min),
        CoordTransformGroup<8>(xg_min, yg_min, convert_dim<X<8>>(x1_min), convert_dim<Y<8>>(y7_min)),
        CoordTransformGroup<9>(xg_min, yg_min, convert_dim<X<9>>(x1_min), convert_dim<Y<9>>(y7_min))};

    // Instantiate the derivatives calculators ---------------------------------------------------
    // SingleInterfaceDerivativesCalculators for interfaces along y (periodic).
    SingleInterfaceDerivativesCalculator<Interface_1_2> const derivatives_calculator_1_2(
            idx_range_x1.take_last(IdxStep<GridX<1>>(10)),
            idx_range_x2.take_first(IdxStep<GridX<2>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_2_3> const derivatives_calculator_2_3(
            idx_range_x2.take_last(IdxStep<GridX<2>>(10)),
            idx_range_x3.take_first(IdxStep<GridX<3>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_3_1> const derivatives_calculator_3_1(
            idx_range_x3.take_last(IdxStep<GridX<3>>(10)),
            idx_range_x1.take_first(IdxStep<GridX<1>>(10)));

    SingleInterfaceDerivativesCalculator<Interface_4_5> const derivatives_calculator_4_5(
            idx_range_x4.take_last(IdxStep<GridX<4>>(10)),
            idx_range_x5.take_first(IdxStep<GridX<5>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_5_6> const derivatives_calculator_5_6(
            idx_range_x5.take_last(IdxStep<GridX<5>>(10)),
            idx_range_x6.take_first(IdxStep<GridX<6>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_6_4> const derivatives_calculator_6_4(
            idx_range_x6.take_last(IdxStep<GridX<6>>(10)),
            idx_range_x4.take_first(IdxStep<GridX<4>>(10)));

    SingleInterfaceDerivativesCalculator<Interface_7_8> const derivatives_calculator_7_8(
            idx_range_x7.take_last(IdxStep<GridX<7>>(10)),
            idx_range_x8.take_first(IdxStep<GridX<8>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_8_9> const derivatives_calculator_8_9(
            idx_range_x8.take_last(IdxStep<GridX<8>>(10)),
            idx_range_x9.take_first(IdxStep<GridX<9>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_9_7> const derivatives_calculator_9_7(
            idx_range_x9.take_last(IdxStep<GridX<9>>(10)),
            idx_range_x7.take_first(IdxStep<GridX<7>>(10)));

    // SingleInterfaceDerivativesCalculators for interfaces along x.
    SingleInterfaceDerivativesCalculator<Interface_1_4> const derivatives_calculator_1_4(
            idx_range_y1.take_first(IdxStep<GridY<1>>(10)),
            idx_range_y4.take_last(IdxStep<GridY<4>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_4_7> const derivatives_calculator_4_7(
            idx_range_y4.take_first(IdxStep<GridY<4>>(10)),
            idx_range_y7.take_last(IdxStep<GridY<7>>(10)));

    SingleInterfaceDerivativesCalculator<Interface_2_5> const derivatives_calculator_2_5(
            idx_range_y2.take_first(IdxStep<GridY<2>>(10)),
            idx_range_y5.take_last(IdxStep<GridY<5>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_5_8> const derivatives_calculator_5_8(
            idx_range_y5.take_first(IdxStep<GridY<5>>(10)),
            idx_range_y8.take_last(IdxStep<GridY<8>>(10)));

    SingleInterfaceDerivativesCalculator<Interface_3_6> const derivatives_calculator_3_6(
            idx_range_y3.take_first(IdxStep<GridY<3>>(10)),
            idx_range_y6.take_last(IdxStep<GridY<6>>(10)));
    SingleInterfaceDerivativesCalculator<Interface_6_9> const derivatives_calculator_6_9(
            idx_range_y6.take_first(IdxStep<GridY<6>>(10)),
            idx_range_y9.take_last(IdxStep<GridY<9>>(10)));

    // TODO: WAITING FOR THE MERGE OF PR 499.
    // constexpr std::size_t nb_chosen_cells = 9;

    // // SingleInterfaceDerivativesCalculators for interfaces along y (periodic).
    // SingleInterfaceDerivativesCalculator<Interface_1_2> const
    //         derivatives_calculator_1_2(idx_range_xy1, idx_range_xy2, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_2_3> const
    //         derivatives_calculator_2_3(idx_range_xy2, idx_range_xy3, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_3_1> const
    //         derivatives_calculator_3_1(idx_range_xy3, idx_range_xy1, nb_chosen_cells);

    // SingleInterfaceDerivativesCalculator<Interface_4_5> const
    //         derivatives_calculator_4_5(idx_range_xy4, idx_range_xy5, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_5_6> const
    //         derivatives_calculator_5_6(idx_range_xy5, idx_range_xy6, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_6_4> const
    //         derivatives_calculator_6_4(idx_range_xy6, idx_range_xy4, nb_chosen_cells);

    // SingleInterfaceDerivativesCalculator<Interface_7_8> const
    //         derivatives_calculator_7_8(idx_range_xy7, idx_range_xy8, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_8_9> const
    //         derivatives_calculator_8_9(idx_range_xy8, idx_range_xy9, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_9_7> const
    //         derivatives_calculator_9_7(idx_range_xy9, idx_range_xy7, nb_chosen_cells);

    // // SingleInterfaceDerivativesCalculators for interfaces along x.
    // SingleInterfaceDerivativesCalculator<Interface_1_4> const
    //         derivatives_calculator_1_4(idx_range_xy1, idx_range_xy4, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_4_7> const
    //         derivatives_calculator_4_7(idx_range_xy4, idx_range_xy7, nb_chosen_cells);

    // SingleInterfaceDerivativesCalculator<Interface_2_5> const
    //         derivatives_calculator_2_5(idx_range_xy2, idx_range_xy5, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_5_8> const
    //         derivatives_calculator_5_8(idx_range_xy5, idx_range_xy8, nb_chosen_cells);

    // SingleInterfaceDerivativesCalculator<Interface_3_6> const
    //         derivatives_calculator_3_6(idx_range_xy3, idx_range_xy6, nb_chosen_cells);
    // SingleInterfaceDerivativesCalculator<Interface_6_9> const
    //         derivatives_calculator_6_9(idx_range_xy6, idx_range_xy9, nb_chosen_cells);

    // Collect the derivative calculators --------------------------------------------------------
    SingleInterfaceDerivativesCalculatorCollection deriv_calculators_collect_123(
            derivatives_calculator_1_2,
            derivatives_calculator_2_3,
            derivatives_calculator_3_1);

    SingleInterfaceDerivativesCalculatorCollection deriv_calculators_collect_456(
            derivatives_calculator_4_5,
            derivatives_calculator_5_6,
            derivatives_calculator_6_4);

    SingleInterfaceDerivativesCalculatorCollection deriv_calculators_collect_789(
            derivatives_calculator_7_8,
            derivatives_calculator_8_9,
            derivatives_calculator_9_7);

    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect_147(derivatives_calculator_1_4, derivatives_calculator_4_7);

    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect_258(derivatives_calculator_2_5, derivatives_calculator_5_8);

    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect_369(derivatives_calculator_3_6, derivatives_calculator_6_9);

    // Collect the index ranges ------------------------------------------------------------------
    MultipatchType<
            IdxRangeOnPatch,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            idx_ranges(
                    idx_range_xy1,
                    idx_range_xy2,
                    idx_range_xy3,
                    idx_range_xy4,
                    idx_range_xy5,
                    idx_range_xy6,
                    idx_range_xy7,
                    idx_range_xy8,
                    idx_range_xy9);

    MultipatchType<IdxRangeOnPatch, Patch1, Patch2, Patch3> idx_ranges_123(idx_ranges);
    MultipatchType<IdxRangeOnPatch, Patch4, Patch5, Patch6> idx_ranges_456(idx_ranges);
    MultipatchType<IdxRangeOnPatch, Patch7, Patch8, Patch9> idx_ranges_789(idx_ranges);

    // For 1|4|7, we artificially add another patch to test if the method will modify the correct patches.
    MultipatchType<IdxRangeOnPatch, Patch1, Patch4, Patch7, Patch2> idx_ranges_147(idx_ranges);
    MultipatchType<IdxRangeOnPatch, Patch2, Patch5, Patch8> idx_ranges_258(idx_ranges);
    MultipatchType<IdxRangeOnPatch, Patch3, Patch6, Patch9> idx_ranges_369(idx_ranges);

    // Instantiate the matrix calculators --------------------------------------------------------
    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<1>,
            ddc::detail::TypeSeq<Patch1, Patch2, Patch3>,
            SingleInterfaceDerivativesCalculatorCollection<
                    Interface_1_2,
                    Interface_2_3,
                    Interface_3_1>>
            matrix_123(idx_ranges_123, deriv_calculators_collect_123);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<4>,
            ddc::detail::TypeSeq<Patch4, Patch5, Patch6>,
            SingleInterfaceDerivativesCalculatorCollection<
                    Interface_4_5,
                    Interface_5_6,
                    Interface_6_4>>
            matrix_456(idx_ranges_456, deriv_calculators_collect_456);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<7>,
            ddc::detail::TypeSeq<Patch7, Patch8, Patch9>,
            SingleInterfaceDerivativesCalculatorCollection<
                    Interface_7_8,
                    Interface_8_9,
                    Interface_9_7>>
            matrix_789(idx_ranges_789, deriv_calculators_collect_789);

    // Test with an extra patch (Patch2) to check it will only take the needed patches.
    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<1>,
            ddc::detail::TypeSeq<Patch1, Patch4, Patch7, Patch2>,
            SingleInterfaceDerivativesCalculatorCollection<Interface_1_4, Interface_4_7>>
            matrix_147(idx_ranges_147, deriv_calculators_collect_147);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<2>,
            ddc::detail::TypeSeq<Patch2, Patch5, Patch8>,
            SingleInterfaceDerivativesCalculatorCollection<Interface_2_5, Interface_5_8>>
            matrix_258(idx_ranges_258, deriv_calculators_collect_258);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<3>,
            ddc::detail::TypeSeq<Patch3, Patch6, Patch9>,
            SingleInterfaceDerivativesCalculatorCollection<Interface_3_6, Interface_6_9>>
            matrix_369(idx_ranges_369, deriv_calculators_collect_369);

    // Instantiate DerivField ====================================================================
    // Instantiate index range slices ------------------------------------------------------------
    IdxRangeSlice<GridX<1>> idx_range_slice_dx1 = get_bound_idx_range_slice(idx_range_x1);
    IdxRangeSlice<GridX<2>> idx_range_slice_dx2 = get_bound_idx_range_slice(idx_range_x2);
    IdxRangeSlice<GridX<3>> idx_range_slice_dx3 = get_bound_idx_range_slice(idx_range_x3);
    IdxRangeSlice<GridX<4>> idx_range_slice_dx4 = get_bound_idx_range_slice(idx_range_x4);
    IdxRangeSlice<GridX<5>> idx_range_slice_dx5 = get_bound_idx_range_slice(idx_range_x5);
    IdxRangeSlice<GridX<6>> idx_range_slice_dx6 = get_bound_idx_range_slice(idx_range_x6);
    IdxRangeSlice<GridX<7>> idx_range_slice_dx7 = get_bound_idx_range_slice(idx_range_x7);
    IdxRangeSlice<GridX<8>> idx_range_slice_dx8 = get_bound_idx_range_slice(idx_range_x8);
    IdxRangeSlice<GridX<9>> idx_range_slice_dx9 = get_bound_idx_range_slice(idx_range_x9);

    IdxRangeSlice<GridY<1>> idx_range_slice_dy1 = get_bound_idx_range_slice(idx_range_y1);
    IdxRangeSlice<GridY<2>> idx_range_slice_dy2 = get_bound_idx_range_slice(idx_range_y2);
    IdxRangeSlice<GridY<3>> idx_range_slice_dy3 = get_bound_idx_range_slice(idx_range_y3);
    IdxRangeSlice<GridY<4>> idx_range_slice_dy4 = get_bound_idx_range_slice(idx_range_y4);
    IdxRangeSlice<GridY<5>> idx_range_slice_dy5 = get_bound_idx_range_slice(idx_range_y5);
    IdxRangeSlice<GridY<6>> idx_range_slice_dy6 = get_bound_idx_range_slice(idx_range_y6);
    IdxRangeSlice<GridY<7>> idx_range_slice_dy7 = get_bound_idx_range_slice(idx_range_y7);
    IdxRangeSlice<GridY<8>> idx_range_slice_dy8 = get_bound_idx_range_slice(idx_range_y8);
    IdxRangeSlice<GridY<9>> idx_range_slice_dy9 = get_bound_idx_range_slice(idx_range_y9);

    // Instantiate DerivField --------------------------------------------------------------------
    DerivFieldMemOnPatch_host<Patch1>
            function_and_derivs_1_alloc(idx_range_xy1, idx_range_slice_dx1, idx_range_slice_dy1);
    DerivFieldMemOnPatch_host<Patch2>
            function_and_derivs_2_alloc(idx_range_xy2, idx_range_slice_dx2, idx_range_slice_dy2);
    DerivFieldMemOnPatch_host<Patch3>
            function_and_derivs_3_alloc(idx_range_xy3, idx_range_slice_dx3, idx_range_slice_dy3);
    DerivFieldMemOnPatch_host<Patch4>
            function_and_derivs_4_alloc(idx_range_xy4, idx_range_slice_dx4, idx_range_slice_dy4);
    DerivFieldMemOnPatch_host<Patch5>
            function_and_derivs_5_alloc(idx_range_xy5, idx_range_slice_dx5, idx_range_slice_dy5);
    DerivFieldMemOnPatch_host<Patch6>
            function_and_derivs_6_alloc(idx_range_xy6, idx_range_slice_dx6, idx_range_slice_dy6);
    DerivFieldMemOnPatch_host<Patch7>
            function_and_derivs_7_alloc(idx_range_xy7, idx_range_slice_dx7, idx_range_slice_dy7);
    DerivFieldMemOnPatch_host<Patch8>
            function_and_derivs_8_alloc(idx_range_xy8, idx_range_slice_dx8, idx_range_slice_dy8);
    DerivFieldMemOnPatch_host<Patch9>
            function_and_derivs_9_alloc(idx_range_xy9, idx_range_slice_dx9, idx_range_slice_dy9);

    DerivFieldOnPatch_host<Patch1> function_and_derivs_1(function_and_derivs_1_alloc);
    DerivFieldOnPatch_host<Patch2> function_and_derivs_2(function_and_derivs_2_alloc);
    DerivFieldOnPatch_host<Patch3> function_and_derivs_3(function_and_derivs_3_alloc);
    DerivFieldOnPatch_host<Patch4> function_and_derivs_4(function_and_derivs_4_alloc);
    DerivFieldOnPatch_host<Patch5> function_and_derivs_5(function_and_derivs_5_alloc);
    DerivFieldOnPatch_host<Patch6> function_and_derivs_6(function_and_derivs_6_alloc);
    DerivFieldOnPatch_host<Patch7> function_and_derivs_7(function_and_derivs_7_alloc);
    DerivFieldOnPatch_host<Patch8> function_and_derivs_8(function_and_derivs_8_alloc);
    DerivFieldOnPatch_host<Patch9> function_and_derivs_9(function_and_derivs_9_alloc);

    // Collect the fields with derivatives.
    MultipatchField<
            DerivFieldOnPatch_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            functions_and_derivs(
                    function_and_derivs_1,
                    function_and_derivs_2,
                    function_and_derivs_3,
                    function_and_derivs_4,
                    function_and_derivs_5,
                    function_and_derivs_6,
                    function_and_derivs_7,
                    function_and_derivs_8,
                    function_and_derivs_9);

    // Instantiate the field of global function values. No derivatives needed.
    host_t<DFieldMem<IdxRange<GridXg, GridYg>>> function_g_alloc(idx_range_xy_g);
    host_t<DField<IdxRange<GridXg, GridYg>>> function_g = get_field(function_g_alloc);

    // Initialise the data =======================================================================
    // --- the function values.
    initialise_all_functions<Xg, Yg>(functions_and_derivs, coord_transforms);
    initialise_2D_function(function_g);

    // --- the first derivatives computed from the function values.
    matrix_123.solve_deriv(functions_and_derivs);
    matrix_456.solve_deriv(functions_and_derivs);
    matrix_789.solve_deriv(functions_and_derivs);

    matrix_147.solve_deriv(functions_and_derivs);
    matrix_258.solve_deriv(functions_and_derivs);
    matrix_369.solve_deriv(functions_and_derivs);

    // --- the cross-derivatives computed from the first derivatives.
    matrix_147.solve_cross_deriv(functions_and_derivs);
    matrix_258.solve_cross_deriv(functions_and_derivs);
    matrix_369.solve_cross_deriv(functions_and_derivs);

    // Test the values of the derivatives ========================================================
    // --- Define an equivalent global spline.
    // Build global spline representation ---
    SplineRThetagBuilder builder_g(idx_range_xy_g);

    host_t<DFieldMem<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef_alloc(
            builder_g.batched_spline_domain(idx_range_xy_g));
    host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef
            = get_field(function_g_coef_alloc);

    builder_g(function_g_coef, get_const_field(function_g));

    // Global spline evaluator ---
    ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymin_g(yg_min);
    ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymax_g(yg_max);
    ddc::PeriodicExtrapolationRule<Xg> bc_x_g;
    SplineRThetagEvaluator evaluator_g(bc_x_g, bc_x_g, bc_ymin_g, bc_ymax_g);

    // Check each derivatives ---
    using PatchSeqLowerBound = ddc::detail::TypeSeq<Patch7, Patch8, Patch9>;
    using PatchSeqUpperBound = ddc::detail::TypeSeq<Patch1, Patch2, Patch3>;

    check_all_x_derivatives(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            coord_transforms,
            1e-4);
    check_all_y_derivatives<PatchSeqLowerBound, PatchSeqUpperBound>(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            coord_transforms,
            1e-4);
    check_all_xy_derivatives<PatchSeqLowerBound, PatchSeqUpperBound>(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            coord_transforms,
            1e-4);

    // Check the whole spline representations ---
    check_all_spline_representation_agreement<PatchSeqLowerBound, PatchSeqUpperBound>(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            coord_transforms,
            1e-7);
}
