// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "../test_utils.hpp"
#include "interface_derivatives/interface_derivatives_test_utils.hpp"
#include "interface_derivatives/interfaces_derivative_calculator_test_utils.hpp"

#include "9patches_2d_periodic_strips_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "interface.hpp"
#include "interface_derivative_coefficients.hpp"
#include "interfaces_derivative_calculator.hpp"
#include "interfaces_exact_derivative_calculator.hpp"
#include "linear_coord_transform.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field_mem.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "orthogonal_coord_transforms.hpp"
#include "types.hpp"
#include "view.hpp"


// Multi-patch tags ---
using namespace periodic_strips_non_uniform_2d_9patches;

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;
using DeviceExecSpace = Kokkos::DefaultExecutionSpace;

// Equivalent global mesh tags ---
struct Xg
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = true;
    using Dual = Xg;
};
struct Yg
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    using Dual = Yg;
};

struct GridXg : NonUniformGridBase<Xg>
{
};
struct GridYg : NonUniformGridBase<Yg>
{
};

struct BSplinesXg : ddc::NonUniformBSplines<Xg, 3, Xg::PERIODIC>
{
};
struct BSplinesYg : ddc::NonUniformBSplines<Yg, 3, Yg::PERIODIC>
{
};

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

// Interpolation points type for the patches.
template <std::size_t PatchIdx>
using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
        BSplinesX<PatchIdx>,
        ddc::SplineBuilderClosure::HERMITE,
        ddc::SplineBuilderClosure::HERMITE>;

template <
        std::size_t PatchIdx,
        ddc::SplineBuilderClosure SplineBuilderClosureMin,
        ddc::SplineBuilderClosure SplineBuilderClosureMax>
using SplineInterpPointsY = ddcHelper::NonUniformInterpolationPoints<
        BSplinesY<PatchIdx>,
        SplineBuilderClosureMin,
        SplineBuilderClosureMax>;

// Interpolation points type for the equivalent global spline.
using SplineInterpPointsXg = ddcHelper::NonUniformInterpolationPoints<
        BSplinesXg,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC>;
using SplineInterpPointsYg = ddcHelper::NonUniformInterpolationPoints<
        BSplinesYg,
        ddc::SplineBuilderClosure::GREVILLE,
        ddc::SplineBuilderClosure::GREVILLE>;


using MPIdxRange = MultipatchType<
        IdxRangeOnPatch,
        Patch1,
        Patch2,
        Patch3,
        Patch4,
        Patch5,
        Patch6,
        Patch7,
        Patch8,
        Patch9>;

using MPDerivFieldMem = MultipatchFieldMem<
        DerivFieldMemOnPatch_host,
        Patch1,
        Patch2,
        Patch3,
        Patch4,
        Patch5,
        Patch6,
        Patch7,
        Patch8,
        Patch9>;

using MPDerivField = MultipatchField<
        DerivFieldOnPatch_host,
        Patch1,
        Patch2,
        Patch3,
        Patch4,
        Patch5,
        Patch6,
        Patch7,
        Patch8,
        Patch9>;

        // Operators on the equivalent global spline.
        using GlobalSplineBuilder = ddc::SplineBuilder2D<
                HostExecSpace,
                typename HostExecSpace::memory_space,
                BSplinesXg,
                BSplinesYg,
                GridXg,
                GridYg,
                ddc::SplineBuilderClosure::PERIODIC,
                ddc::SplineBuilderClosure::PERIODIC,
                ddc::SplineBuilderClosure::GREVILLE,
                ddc::SplineBuilderClosure::GREVILLE,
                ddc::SplineSolver::LAPACK>;

using GlobalSplineEvaluator = ddc::SplineEvaluator2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::PeriodicExtrapolationRule<Xg>,
        ddc::PeriodicExtrapolationRule<Xg>,
        ddc::ConstantExtrapolationRule<Yg>,
        ddc::ConstantExtrapolationRule<Yg>>;


template <bool is_approximation, class Connectivity, class Grid1D, ddc::SplineBuilderClosure SBC>
using InterfacesDerivativeCalculatorType = std::conditional_t<
        is_approximation,
        InterfacesDerivativeCalculator<Connectivity, Grid1D>,
        InterfacesExactDerivativeCalculator<Connectivity, Grid1D, SBC, SBC>>;

template <bool is_approximation, class InterfaceType, class IdxRange1, class IdxRange2>
InterfaceDerivCoeffs<InterfaceType> build_interface(
        IdxRange1 idx_range_1,
        IdxRange2 idx_range_2,
        std::size_t nb_chosen_cells,
        ddc::SplineBuilderClosure BC1,
        ddc::SplineBuilderClosure BC2)
{
    if constexpr (is_approximation) {
        return InterfaceDerivCoeffs<InterfaceType>(idx_range_1, idx_range_2, nb_chosen_cells);
    } else {
        return InterfaceDerivCoeffs<InterfaceType>(idx_range_1, idx_range_2, BC1, BC2);
    }
}

template <std::size_t I>
struct CoordTransformGroup
{
    using XTransform = LinearCoordTransform<Xg, X<I>>;
    using YTransform = LinearCoordTransform<Yg, Y<I>>;
    using XYTransform = OrthogonalCoordTransforms<
            Coord<Yg, Xg>,
            Coord<X<I>, Y<I>>,
            Coord<Xg, Yg>,
            XTransform,
            YTransform>;
    XTransform x_transform;
    YTransform y_transform;
    XYTransform coord_transform;
    CoordTransformGroup()
        : x_transform(Coord<Xg>(0), Coord<X<I>>(0), 1.0)
        , y_transform(Coord<Yg>(0), Coord<Y<I>>(0), 1.0)
        , coord_transform(x_transform, y_transform)
    {
    }
};

template <class is_approximation_literal>
struct Multipatch9PatchesPeriodicStripsNonUniformTest : public ::testing::Test
{
    static constexpr bool is_approximation = is_approximation_literal::value;

protected:
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

    static constexpr ddc::SplineBuilderClosure BcH = ddc::SplineBuilderClosure::HERMITE;
    static constexpr ddc::SplineBuilderClosure BcG = ddc::SplineBuilderClosure::GREVILLE;

protected:
    // Index ranges
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

    const MPIdxRange idx_ranges;

    const IdxRange<GridXg> idx_range_xg;
    const IdxRange<GridYg> idx_range_yg;
    const IdxRange<GridXg, GridYg> idx_range_xy_g;

    // Coordinate transformations
    std::tuple<
            CoordTransformGroup<1>,
            CoordTransformGroup<2>,
            CoordTransformGroup<3>,
            CoordTransformGroup<4>,
            CoordTransformGroup<5>,
            CoordTransformGroup<6>,
            CoordTransformGroup<7>,
            CoordTransformGroup<8>,
            CoordTransformGroup<9>> const coord_transforms;

    // Global spline operators
    GlobalSplineBuilder const builder_g;
    ddc::ConstantExtrapolationRule<Yg> const bc_ymin_g;
    ddc::ConstantExtrapolationRule<Yg> const bc_ymax_g;
    ddc::PeriodicExtrapolationRule<Xg> const bc_x_g;
    GlobalSplineEvaluator const evaluator_g;

public:
    Multipatch9PatchesPeriodicStripsNonUniformTest()
        // Index ranges
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
        , idx_ranges(
                  idx_range_xy1,
                  idx_range_xy2,
                  idx_range_xy3,
                  idx_range_xy4,
                  idx_range_xy5,
                  idx_range_xy6,
                  idx_range_xy7,
                  idx_range_xy8,
                  idx_range_xy9)
        , idx_range_xg(SplineInterpPointsXg::get_domain<GridXg>())
        , idx_range_yg(SplineInterpPointsYg::get_domain<GridYg>())
        , idx_range_xy_g(idx_range_xg, idx_range_yg)
        // Coordinate transformations
        , coord_transforms(
                  CoordTransformGroup<1>(),
                  CoordTransformGroup<2>(),
                  CoordTransformGroup<3>(),
                  CoordTransformGroup<4>(),
                  CoordTransformGroup<5>(),
                  CoordTransformGroup<6>(),
                  CoordTransformGroup<7>(),
                  CoordTransformGroup<8>(),
                  CoordTransformGroup<9>())
        // Global spline operators
        , builder_g(idx_range_xy_g)
        , bc_ymin_g(yg_min)
        , bc_ymax_g(yg_max)
        , bc_x_g()
        , evaluator_g(bc_x_g, bc_x_g, bc_ymin_g, bc_ymax_g)
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

        interpolation_points_x147.pop_back();
        interpolation_points_x258.pop_back();
        interpolation_points_x369.pop_back();

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
        append(break_points_xg, break_points_x147);
        append(interpolation_points_xg, interpolation_points_x147);

        break_points_x258.pop_back();
        append(break_points_xg, break_points_x258);
        append(interpolation_points_xg, interpolation_points_x258);

        append(break_points_xg, break_points_x369);
        append(interpolation_points_xg, interpolation_points_x369);

        std::vector<Coord<Yg>> break_points_yg;
        std::vector<Coord<Yg>> interpolation_points_yg;

        break_points_y789.pop_back();
        interpolation_points_y789.pop_back();
        append(break_points_yg, break_points_y789);
        append(interpolation_points_yg, interpolation_points_y789);

        break_points_y456.pop_back();
        interpolation_points_y456.pop_back();
        append(break_points_yg, break_points_y456);
        append(interpolation_points_yg, interpolation_points_y456);

        append(break_points_yg, break_points_y123);
        append(interpolation_points_yg, interpolation_points_y123);

        ddc::init_discrete_space<BSplinesXg>(break_points_xg);
        ddc::init_discrete_space<BSplinesYg>(break_points_yg);

        ddc::init_discrete_space<GridXg>(
                SplineInterpPointsXg::get_sampling<GridXg>(interpolation_points_xg));
        ddc::init_discrete_space<GridYg>(
                SplineInterpPointsYg::get_sampling<GridYg>(interpolation_points_yg));
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

using is_approximation = std::tuple<std::true_type, std::false_type>;
using TestTypes = tuple_to_types_t<is_approximation>;
TYPED_TEST_SUITE(Multipatch9PatchesPeriodicStripsNonUniformTest, TestTypes);
