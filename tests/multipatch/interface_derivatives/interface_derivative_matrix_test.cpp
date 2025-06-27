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
#include "interface_derivatives_test_utils.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"


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


Coord<Xg, Yg> get_global_coord(Coord<X, Y> const& local_coord)
{
    double const x = ddc::select<X>(local_coord);
    double const y = ddc::select<Y>(local_coord);
    return Coord<Xg, Yg>(x, y);
}

/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y).
template <class Grid1, class Grid2>
void initialise_2D_function(host_t<DField<IdxRange<Grid1, Grid2>>> function)
{
    ddc::for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        double const xg = ddc::coordinate(Idx<Grid1>(idx));
        double const yg = ddc::coordinate(Idx<Grid2>(idx));
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI) * Kokkos::sin(yg);
    });
}


/// @brief Initialise the y-derivatives from the global spline.
template <class Grid1, class Grid2>
void initialise_2D_y_derivative(
        host_t<DField<IdxRange<Grid1, ddc::Deriv<typename Grid2::continuous_dimension_type>>>>
                deriv_y,
        Idx<Grid2> const& idx_y,
        SplineRThetagEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
{
    ddc::for_each(
            get_idx_range(deriv_y),
            [&](Idx<Grid1, ddc::Deriv<typename Grid2::continuous_dimension_type>> idx) {
                // Get the coordinate on the equivalent global mesh.
                Idx<Grid1, Grid2> idx_xy(Idx<Grid1>(idx), idx_y);
                Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx_xy)));
                deriv_y(idx) = evaluator_g.deriv_dim_2(interface_coord, function_g_coef);
            });
}


/// @brief Initialise the cross derivatives from the global spline.
template <class Grid1, class Grid2>
void initialise_2D_xy_derivative(
        host_t<DField<IdxRange<
                ddc::Deriv<typename Grid1::continuous_dimension_type>,
                ddc::Deriv<typename Grid2::continuous_dimension_type>>>> deriv_xy,
        Idx<Grid1> const& idx_x,
        Idx<Grid2> const& idx_y,
        SplineRThetagEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
{
    // Get the coordinate on the equivalent global mesh.
    Idx<Grid1, Grid2> idx(idx_x, idx_y);
    Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx)));
    deriv_xy(get_idx_range(deriv_xy).front())
            = evaluator_g.deriv_1_and_2(interface_coord, function_g_coef);
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

    const IdxRange<GridXg> idx_range_xg;
    const IdxRange<GridYg> idx_range_yg;
    const IdxRange<GridXg, GridYg> idx_range_xy_g;


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
        , idx_range_xg(SplineInterpPointsXg::get_domain<GridXg>())
        , idx_range_yg(SplineInterpPointsYg::get_domain<GridYg>())
        , idx_range_xy_g(idx_range_xg, idx_range_yg)
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
    template <class... Patches>
    void check_all_x_derivatives(
            MultipatchField<Deriv1_OnPatch_2D_host, Patches...> const& local_derivs_min,
            MultipatchField<Deriv1_OnPatch_2D_host, Patches...> const& local_derivs_max,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges)
    {
        (check_x_derivatives<Patches>(
                 local_derivs_min.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 typename Patches::IdxRange1(idx_ranges.template get<Patches>()).front()),
         ...);
        (check_x_derivatives<Patches>(
                 local_derivs_max.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 typename Patches::IdxRange1(idx_ranges.template get<Patches>()).back()),
         ...);
    }


    template <class Patch>
    void check_x_derivatives(
            Deriv1_OnPatch_2D_host<Patch> const& local_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::Idx1 const& idx_perp)
    {
        typename Patch::IdxRange2 idx_range_par(get_idx_range(local_derivs));
        IdxRange<ddc::Deriv<typename Patch::Dim1>> idx_range_deriv(get_idx_range(local_derivs));
        Idx<ddc::Deriv<typename Patch::Dim1>> idx_deriv = idx_range_deriv.front();

        ddc::for_each(idx_range_par, [&](typename Patch::Idx2 const& idx_par) {
            typename Patch::Idx12 idx(idx_perp, idx_par);
            Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx)));
            // std::cout << "interface coordinate: " << interface_coord << std::endl;

            double const local_deriv = local_derivs(idx_deriv, idx_par);
            double const global_deriv = evaluator_g.deriv_dim_1(interface_coord, function_g_coef);

            EXPECT_NEAR(local_deriv, global_deriv, 2e-14);
        });
    }

    template <class... Patches>
    void check_all_y_derivatives(
            MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& local_derivs_min,
            MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& local_derivs_max,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges)
    {
        (check_y_derivatives<Patches, true>(
                 local_derivs_min.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 typename Patches::IdxRange2(idx_ranges.template get<Patches>()).front()),
         ...);
        (check_y_derivatives<Patches, false>(
                 local_derivs_max.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 typename Patches::IdxRange2(idx_ranges.template get<Patches>()).back()),
         ...);
    }


    template <class Patch, bool is_deriv_min>
    void check_y_derivatives(
            Deriv2_OnPatch_2D_host<Patch> const& local_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::Idx2 const& idx_perp)
    {
        typename Patch::IdxRange1 idx_range_par(get_idx_range(local_derivs));
        IdxRange<ddc::Deriv<typename Patch::Dim2>> idx_range_deriv(get_idx_range(local_derivs));
        Idx<ddc::Deriv<typename Patch::Dim2>> idx_deriv = idx_range_deriv.front();

        ddc::for_each(idx_range_par, [&](typename Patch::Idx1 const& idx_par) {
            typename Patch::Idx12 idx(idx_par, idx_perp);
            Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx)));

            double const local_deriv = local_derivs(idx_par, idx_deriv);
            double const global_deriv = evaluator_g.deriv_dim_2(interface_coord, function_g_coef);

            EXPECT_NEAR(local_deriv, global_deriv, 2e-14);
        });
    }

    template <std::size_t Index>
    using PatchI = ddc::type_seq_element_t<Index, PatchOrdering>;

    template <class... Patches, std::size_t... I>
    void check_all_xy_derivatives(
            std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& local_derivs_min_min,
            std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& local_derivs_max_min,
            std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& local_derivs_min_max,
            std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& local_derivs_max_max,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            std::integer_sequence<std::size_t, I...>)
    {
        (check_xy_derivatives<PatchI<I>, true>(
                 std::get<I>(local_derivs_min_min),
                 evaluator_g,
                 function_g_coef,
                 typename PatchI<I>::IdxRange1(idx_ranges.template get<PatchI<I>>()).front(),
                 typename PatchI<I>::IdxRange2(idx_ranges.template get<PatchI<I>>()).front()),
         ...);
        (check_xy_derivatives<Patches, true>(
                 std::get<I>(local_derivs_max_min),
                 evaluator_g,
                 function_g_coef,
                 typename PatchI<I>::IdxRange1(idx_ranges.template get<PatchI<I>>()).back(),
                 typename PatchI<I>::IdxRange2(idx_ranges.template get<PatchI<I>>()).front()),
         ...);
        (check_xy_derivatives<Patches, false>(
                 std::get<I>(local_derivs_min_max),
                 evaluator_g,
                 function_g_coef,
                 typename PatchI<I>::IdxRange1(idx_ranges.template get<PatchI<I>>()).front(),
                 typename PatchI<I>::IdxRange2(idx_ranges.template get<PatchI<I>>()).back()),
         ...);
        (check_xy_derivatives<Patches, false>(
                 std::get<I>(local_derivs_max_max),
                 evaluator_g,
                 function_g_coef,
                 typename PatchI<I>::IdxRange1(idx_ranges.template get<PatchI<I>>()).back(),
                 typename PatchI<I>::IdxRange2(idx_ranges.template get<PatchI<I>>()).back()),
         ...);
    }


    template <class Patch, bool is_deriv2_min>
    void check_xy_derivatives(
            Deriv12_OnPatch_2D_host<Patch> const& local_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::Idx1 const& idx_1,
            typename Patch::Idx2 const& idx_2)
    {
        typename Patch::Idx12 idx(idx_1, idx_2);
        Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx)));

        Idx<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2>> idx_deriv(
                get_idx_range(local_derivs).front());

        double const local_deriv = local_derivs(idx_deriv);
        double const global_deriv = evaluator_g.deriv_1_and_2(interface_coord, function_g_coef);

        EXPECT_NEAR(local_deriv, global_deriv, 6e-14);
    }
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
    // Instantiate the derivatives calculators ---------------------------------------------------
    // SingleInterfaceDerivativesCalculators for interfaces along y (periodic).
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

    // SingleInterfaceDerivativesCalculators for interfaces along x.
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



    // Order in sequences ------------------------------------------------------------------------
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


    // Instantiate the matrix calculators --------------------------------------------------------
    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<1>,
            DConstFieldOnPatch_host,
            // Deriv1_OnPatch_2D_host,
            true,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            Kokkos::DefaultHostExecutionSpace,
            Patch1,
            Patch2,
            Patch3>
            matrix_123(idx_ranges_123, derivative_calculators_123);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<4>,
            DConstFieldOnPatch_host,
            // Deriv1_OnPatch_2D_host,
            true,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            Kokkos::DefaultHostExecutionSpace,
            Patch4,
            Patch5,
            Patch6>
            matrix_456(idx_ranges_456, derivative_calculators_456);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<7>,
            DConstFieldOnPatch_host,
            // Deriv1_OnPatch_2D_host,
            true,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            Kokkos::DefaultHostExecutionSpace,
            Patch7,
            Patch8,
            Patch9>
            matrix_789(idx_ranges_789, derivative_calculators_789);



    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<1>,
            DConstFieldOnPatch_host,
            // Deriv2_OnPatch_2D_host,
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
            DConstFieldOnPatch_host,
            // Deriv2_OnPatch_2D_host,
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
            DConstFieldOnPatch_host,
            // Deriv2_OnPatch_2D_host,
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

    // std::cout << boost::typeindex::type_id_with_cvr<
    //                      ddc::type_seq_element_t<0, interface_collection>>()
    //                      .pretty_name()
    //           << std::endl;
    // std::cout << boost::typeindex::type_id_with_cvr<
    //                      ddc::type_seq_element_t<1, interface_collection>>()
    //                      .pretty_name()
    //           << std::endl;
    // std::cout << boost::typeindex::type_id_with_cvr<
    //                      ddc::type_seq_element_t<2, interface_collection>>()
    //                      .pretty_name()
    //           << std::endl;
    // std::cout << boost::typeindex::type_id_with_cvr<ddc::type_seq_element_t<3,interface_collection>>().pretty_name() << std::endl;


    // Instantiate test function values ==========================================================
    // --- patch 1
    host_t<DFieldMem<Patch1::IdxRange12>> function_1_alloc(idx_range_xy1);
    host_t<DField<Patch1::IdxRange12>> function_1 = get_field(function_1_alloc);

    // --- patch 2
    host_t<DFieldMem<Patch2::IdxRange12>> function_2_alloc(idx_range_xy2);
    host_t<DField<Patch2::IdxRange12>> function_2 = get_field(function_2_alloc);

    // --- patch 3
    host_t<DFieldMem<Patch3::IdxRange12>> function_3_alloc(idx_range_xy3);
    host_t<DField<Patch3::IdxRange12>> function_3 = get_field(function_3_alloc);

    // --- patch 4
    host_t<DFieldMem<Patch4::IdxRange12>> function_4_alloc(idx_range_xy4);
    host_t<DField<Patch4::IdxRange12>> function_4 = get_field(function_4_alloc);

    // --- patch 5
    host_t<DFieldMem<Patch5::IdxRange12>> function_5_alloc(idx_range_xy5);
    host_t<DField<Patch5::IdxRange12>> function_5 = get_field(function_5_alloc);

    // --- patch 6
    host_t<DFieldMem<Patch6::IdxRange12>> function_6_alloc(idx_range_xy6);
    host_t<DField<Patch6::IdxRange12>> function_6 = get_field(function_6_alloc);

    // --- patch 7
    host_t<DFieldMem<Patch7::IdxRange12>> function_7_alloc(idx_range_xy7);
    host_t<DField<Patch7::IdxRange12>> function_7 = get_field(function_7_alloc);

    // --- patch 8
    host_t<DFieldMem<Patch8::IdxRange12>> function_8_alloc(idx_range_xy8);
    host_t<DField<Patch8::IdxRange12>> function_8 = get_field(function_8_alloc);

    // --- patch 9
    host_t<DFieldMem<Patch9::IdxRange12>> function_9_alloc(idx_range_xy9);
    host_t<DField<Patch9::IdxRange12>> function_9 = get_field(function_9_alloc);

    // --- global
    host_t<DFieldMem<IdxRange<GridXg, GridYg>>> function_g_alloc(idx_range_xy_g);
    host_t<DField<IdxRange<GridXg, GridYg>>> function_g = get_field(function_g_alloc);


    // Collect the function values.
    MultipatchField<DFieldOnPatch_host, Patch1, Patch2, Patch3>
            functions_123(function_1, function_2, function_3);
    MultipatchField<DFieldOnPatch_host, Patch4, Patch5, Patch6>
            functions_456(function_4, function_5, function_6);
    MultipatchField<DFieldOnPatch_host, Patch7, Patch8, Patch9>
            functions_789(function_7, function_8, function_9);

    MultipatchField<DFieldOnPatch_host, Patch1, Patch4, Patch7>
            functions_147(function_1, function_4, function_7);
    MultipatchField<DFieldOnPatch_host, Patch2, Patch5, Patch8>
            functions_258(function_2, function_5, function_8);
    MultipatchField<DFieldOnPatch_host, Patch3, Patch6, Patch9>
            functions_369(function_3, function_6, function_9);


    // Instantiate derivatives of the test functions ==============================================
    // --- derivatives along x.
    Idx<ddc::Deriv<X>> first_deriv_x(1);
    IdxStep<ddc::Deriv<X>> n_deriv_x(1);
    IdxRange<ddc::Deriv<X>> deriv_x_idx_range(first_deriv_x, n_deriv_x);

    IdxRange<ddc::Deriv<X>, GridY<1>> derivs_x_idx_range_1(deriv_x_idx_range, idx_range_y1);
    IdxRange<ddc::Deriv<X>, GridY<2>> derivs_x_idx_range_2(deriv_x_idx_range, idx_range_y2);
    IdxRange<ddc::Deriv<X>, GridY<3>> derivs_x_idx_range_3(deriv_x_idx_range, idx_range_y3);
    IdxRange<ddc::Deriv<X>, GridY<4>> derivs_x_idx_range_4(deriv_x_idx_range, idx_range_y4);
    IdxRange<ddc::Deriv<X>, GridY<5>> derivs_x_idx_range_5(deriv_x_idx_range, idx_range_y5);
    IdxRange<ddc::Deriv<X>, GridY<6>> derivs_x_idx_range_6(deriv_x_idx_range, idx_range_y6);
    IdxRange<ddc::Deriv<X>, GridY<7>> derivs_x_idx_range_7(deriv_x_idx_range, idx_range_y7);
    IdxRange<ddc::Deriv<X>, GridY<8>> derivs_x_idx_range_8(deriv_x_idx_range, idx_range_y8);
    IdxRange<ddc::Deriv<X>, GridY<9>> derivs_x_idx_range_9(deriv_x_idx_range, idx_range_y9);

    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<1>>>> derivs_xmin_1_alloc(derivs_x_idx_range_1);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<2>>>> derivs_xmin_2_alloc(derivs_x_idx_range_2);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<3>>>> derivs_xmin_3_alloc(derivs_x_idx_range_3);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<4>>>> derivs_xmin_4_alloc(derivs_x_idx_range_4);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<5>>>> derivs_xmin_5_alloc(derivs_x_idx_range_5);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<6>>>> derivs_xmin_6_alloc(derivs_x_idx_range_6);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<7>>>> derivs_xmin_7_alloc(derivs_x_idx_range_7);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<8>>>> derivs_xmin_8_alloc(derivs_x_idx_range_8);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<9>>>> derivs_xmin_9_alloc(derivs_x_idx_range_9);

    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<1>>>> derivs_xmin_1
            = get_field(derivs_xmin_1_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<2>>>> derivs_xmin_2
            = get_field(derivs_xmin_2_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<3>>>> derivs_xmin_3
            = get_field(derivs_xmin_3_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<4>>>> derivs_xmin_4
            = get_field(derivs_xmin_4_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<5>>>> derivs_xmin_5
            = get_field(derivs_xmin_5_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<6>>>> derivs_xmin_6
            = get_field(derivs_xmin_6_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<7>>>> derivs_xmin_7
            = get_field(derivs_xmin_7_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<8>>>> derivs_xmin_8
            = get_field(derivs_xmin_8_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<9>>>> derivs_xmin_9
            = get_field(derivs_xmin_9_alloc);


    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<1>>>> derivs_xmax_1_alloc(derivs_x_idx_range_1);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<2>>>> derivs_xmax_2_alloc(derivs_x_idx_range_2);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<3>>>> derivs_xmax_3_alloc(derivs_x_idx_range_3);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<4>>>> derivs_xmax_4_alloc(derivs_x_idx_range_4);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<5>>>> derivs_xmax_5_alloc(derivs_x_idx_range_5);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<6>>>> derivs_xmax_6_alloc(derivs_x_idx_range_6);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<7>>>> derivs_xmax_7_alloc(derivs_x_idx_range_7);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<8>>>> derivs_xmax_8_alloc(derivs_x_idx_range_8);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, GridY<9>>>> derivs_xmax_9_alloc(derivs_x_idx_range_9);

    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<1>>>> derivs_xmax_1
            = get_field(derivs_xmax_1_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<2>>>> derivs_xmax_2
            = get_field(derivs_xmax_2_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<3>>>> derivs_xmax_3
            = get_field(derivs_xmax_3_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<4>>>> derivs_xmax_4
            = get_field(derivs_xmax_4_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<5>>>> derivs_xmax_5
            = get_field(derivs_xmax_5_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<6>>>> derivs_xmax_6
            = get_field(derivs_xmax_6_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<7>>>> derivs_xmax_7
            = get_field(derivs_xmax_7_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<8>>>> derivs_xmax_8
            = get_field(derivs_xmax_8_alloc);
    host_t<DField<IdxRange<ddc::Deriv<X>, GridY<9>>>> derivs_xmax_9
            = get_field(derivs_xmax_9_alloc);


    // --- derivatives along y.
    Idx<ddc::Deriv<Y>> first_deriv_y(1);
    IdxStep<ddc::Deriv<Y>> n_deriv_y(1);
    IdxRange<ddc::Deriv<Y>> deriv_y_idx_range(first_deriv_y, n_deriv_y);

    IdxRange<GridX<1>, ddc::Deriv<Y>> derivs_y_idx_range_1(idx_range_x1, deriv_y_idx_range);
    IdxRange<GridX<2>, ddc::Deriv<Y>> derivs_y_idx_range_2(idx_range_x2, deriv_y_idx_range);
    IdxRange<GridX<3>, ddc::Deriv<Y>> derivs_y_idx_range_3(idx_range_x3, deriv_y_idx_range);
    IdxRange<GridX<4>, ddc::Deriv<Y>> derivs_y_idx_range_4(idx_range_x4, deriv_y_idx_range);
    IdxRange<GridX<5>, ddc::Deriv<Y>> derivs_y_idx_range_5(idx_range_x5, deriv_y_idx_range);
    IdxRange<GridX<6>, ddc::Deriv<Y>> derivs_y_idx_range_6(idx_range_x6, deriv_y_idx_range);
    IdxRange<GridX<7>, ddc::Deriv<Y>> derivs_y_idx_range_7(idx_range_x7, deriv_y_idx_range);
    IdxRange<GridX<8>, ddc::Deriv<Y>> derivs_y_idx_range_8(idx_range_x8, deriv_y_idx_range);
    IdxRange<GridX<9>, ddc::Deriv<Y>> derivs_y_idx_range_9(idx_range_x9, deriv_y_idx_range);

    host_t<DFieldMem<IdxRange<GridX<1>, ddc::Deriv<Y>>>> derivs_ymin_1_alloc(derivs_y_idx_range_1);
    host_t<DFieldMem<IdxRange<GridX<2>, ddc::Deriv<Y>>>> derivs_ymin_2_alloc(derivs_y_idx_range_2);
    host_t<DFieldMem<IdxRange<GridX<3>, ddc::Deriv<Y>>>> derivs_ymin_3_alloc(derivs_y_idx_range_3);
    host_t<DFieldMem<IdxRange<GridX<4>, ddc::Deriv<Y>>>> derivs_ymin_4_alloc(derivs_y_idx_range_4);
    host_t<DFieldMem<IdxRange<GridX<5>, ddc::Deriv<Y>>>> derivs_ymin_5_alloc(derivs_y_idx_range_5);
    host_t<DFieldMem<IdxRange<GridX<6>, ddc::Deriv<Y>>>> derivs_ymin_6_alloc(derivs_y_idx_range_6);
    host_t<DFieldMem<IdxRange<GridX<7>, ddc::Deriv<Y>>>> derivs_ymin_7_alloc(derivs_y_idx_range_7);
    host_t<DFieldMem<IdxRange<GridX<8>, ddc::Deriv<Y>>>> derivs_ymin_8_alloc(derivs_y_idx_range_8);
    host_t<DFieldMem<IdxRange<GridX<9>, ddc::Deriv<Y>>>> derivs_ymin_9_alloc(derivs_y_idx_range_9);

    host_t<DField<IdxRange<GridX<1>, ddc::Deriv<Y>>>> derivs_ymin_1
            = get_field(derivs_ymin_1_alloc);
    host_t<DField<IdxRange<GridX<2>, ddc::Deriv<Y>>>> derivs_ymin_2
            = get_field(derivs_ymin_2_alloc);
    host_t<DField<IdxRange<GridX<3>, ddc::Deriv<Y>>>> derivs_ymin_3
            = get_field(derivs_ymin_3_alloc);
    host_t<DField<IdxRange<GridX<4>, ddc::Deriv<Y>>>> derivs_ymin_4
            = get_field(derivs_ymin_4_alloc);
    host_t<DField<IdxRange<GridX<5>, ddc::Deriv<Y>>>> derivs_ymin_5
            = get_field(derivs_ymin_5_alloc);
    host_t<DField<IdxRange<GridX<6>, ddc::Deriv<Y>>>> derivs_ymin_6
            = get_field(derivs_ymin_6_alloc);
    host_t<DField<IdxRange<GridX<7>, ddc::Deriv<Y>>>> derivs_ymin_7
            = get_field(derivs_ymin_7_alloc);
    host_t<DField<IdxRange<GridX<8>, ddc::Deriv<Y>>>> derivs_ymin_8
            = get_field(derivs_ymin_8_alloc);
    host_t<DField<IdxRange<GridX<9>, ddc::Deriv<Y>>>> derivs_ymin_9
            = get_field(derivs_ymin_9_alloc);


    host_t<DFieldMem<IdxRange<GridX<1>, ddc::Deriv<Y>>>> derivs_ymax_1_alloc(derivs_y_idx_range_1);
    host_t<DFieldMem<IdxRange<GridX<2>, ddc::Deriv<Y>>>> derivs_ymax_2_alloc(derivs_y_idx_range_2);
    host_t<DFieldMem<IdxRange<GridX<3>, ddc::Deriv<Y>>>> derivs_ymax_3_alloc(derivs_y_idx_range_3);
    host_t<DFieldMem<IdxRange<GridX<4>, ddc::Deriv<Y>>>> derivs_ymax_4_alloc(derivs_y_idx_range_4);
    host_t<DFieldMem<IdxRange<GridX<5>, ddc::Deriv<Y>>>> derivs_ymax_5_alloc(derivs_y_idx_range_5);
    host_t<DFieldMem<IdxRange<GridX<6>, ddc::Deriv<Y>>>> derivs_ymax_6_alloc(derivs_y_idx_range_6);
    host_t<DFieldMem<IdxRange<GridX<7>, ddc::Deriv<Y>>>> derivs_ymax_7_alloc(derivs_y_idx_range_7);
    host_t<DFieldMem<IdxRange<GridX<8>, ddc::Deriv<Y>>>> derivs_ymax_8_alloc(derivs_y_idx_range_8);
    host_t<DFieldMem<IdxRange<GridX<9>, ddc::Deriv<Y>>>> derivs_ymax_9_alloc(derivs_y_idx_range_9);

    host_t<DField<IdxRange<GridX<1>, ddc::Deriv<Y>>>> derivs_ymax_1
            = get_field(derivs_ymax_1_alloc);
    host_t<DField<IdxRange<GridX<2>, ddc::Deriv<Y>>>> derivs_ymax_2
            = get_field(derivs_ymax_2_alloc);
    host_t<DField<IdxRange<GridX<3>, ddc::Deriv<Y>>>> derivs_ymax_3
            = get_field(derivs_ymax_3_alloc);
    host_t<DField<IdxRange<GridX<4>, ddc::Deriv<Y>>>> derivs_ymax_4
            = get_field(derivs_ymax_4_alloc);
    host_t<DField<IdxRange<GridX<5>, ddc::Deriv<Y>>>> derivs_ymax_5
            = get_field(derivs_ymax_5_alloc);
    host_t<DField<IdxRange<GridX<6>, ddc::Deriv<Y>>>> derivs_ymax_6
            = get_field(derivs_ymax_6_alloc);
    host_t<DField<IdxRange<GridX<7>, ddc::Deriv<Y>>>> derivs_ymax_7
            = get_field(derivs_ymax_7_alloc);
    host_t<DField<IdxRange<GridX<8>, ddc::Deriv<Y>>>> derivs_ymax_8
            = get_field(derivs_ymax_8_alloc);
    host_t<DField<IdxRange<GridX<9>, ddc::Deriv<Y>>>> derivs_ymax_9
            = get_field(derivs_ymax_9_alloc);


    // Collect the derivatives.
    MultipatchField<
            Deriv1_OnPatch_2D_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            derivs_xmin(
                    derivs_xmin_1,
                    derivs_xmin_2,
                    derivs_xmin_3,
                    derivs_xmin_4,
                    derivs_xmin_5,
                    derivs_xmin_6,
                    derivs_xmin_7,
                    derivs_xmin_8,
                    derivs_xmin_9);
    MultipatchField<
            Deriv1_OnPatch_2D_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            derivs_xmax(
                    derivs_xmax_1,
                    derivs_xmax_2,
                    derivs_xmax_3,
                    derivs_xmax_4,
                    derivs_xmax_5,
                    derivs_xmax_6,
                    derivs_xmax_7,
                    derivs_xmax_8,
                    derivs_xmax_9);

    MultipatchField<
            Deriv2_OnPatch_2D_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            derivs_ymin(
                    derivs_ymin_1,
                    derivs_ymin_2,
                    derivs_ymin_3,
                    derivs_ymin_4,
                    derivs_ymin_5,
                    derivs_ymin_6,
                    derivs_ymin_7,
                    derivs_ymin_8,
                    derivs_ymin_9);

    MultipatchField<
            Deriv2_OnPatch_2D_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            derivs_ymax(
                    derivs_ymax_1,
                    derivs_ymax_2,
                    derivs_ymax_3,
                    derivs_ymax_4,
                    derivs_ymax_5,
                    derivs_ymax_6,
                    derivs_ymax_7,
                    derivs_ymax_8,
                    derivs_ymax_9);


    MultipatchField<Deriv1_OnPatch_2D_host, Patch1, Patch2, Patch3>
            derivs_xmin_123(derivs_xmin_1, derivs_xmin_2, derivs_xmin_3);
    MultipatchField<Deriv1_OnPatch_2D_host, Patch4, Patch5, Patch6>
            derivs_xmin_456(derivs_xmin_4, derivs_xmin_5, derivs_xmin_6);
    MultipatchField<Deriv1_OnPatch_2D_host, Patch7, Patch8, Patch9>
            derivs_xmin_789(derivs_xmin_7, derivs_xmin_8, derivs_xmin_9);

    MultipatchField<Deriv1_OnPatch_2D_host, Patch1, Patch2, Patch3>
            derivs_xmax_123(derivs_xmax_1, derivs_xmax_2, derivs_xmax_3);
    MultipatchField<Deriv1_OnPatch_2D_host, Patch4, Patch5, Patch6>
            derivs_xmax_456(derivs_xmax_4, derivs_xmax_5, derivs_xmax_6);
    MultipatchField<Deriv1_OnPatch_2D_host, Patch7, Patch8, Patch9>
            derivs_xmax_789(derivs_xmax_7, derivs_xmax_8, derivs_xmax_9);


    MultipatchField<Deriv2_OnPatch_2D_host, Patch1, Patch4, Patch7>
            derivs_ymin_147(derivs_ymin_1, derivs_ymin_4, derivs_ymin_7);
    MultipatchField<Deriv2_OnPatch_2D_host, Patch2, Patch5, Patch8>
            derivs_ymin_258(derivs_ymin_2, derivs_ymin_5, derivs_ymin_8);
    MultipatchField<Deriv2_OnPatch_2D_host, Patch3, Patch6, Patch9>
            derivs_ymin_369(derivs_ymin_3, derivs_ymin_6, derivs_ymin_9);

    MultipatchField<Deriv2_OnPatch_2D_host, Patch1, Patch4, Patch7>
            derivs_ymax_147(derivs_ymax_1, derivs_ymax_4, derivs_ymax_7);
    MultipatchField<Deriv2_OnPatch_2D_host, Patch2, Patch5, Patch8>
            derivs_ymax_258(derivs_ymax_2, derivs_ymax_5, derivs_ymax_8);
    MultipatchField<Deriv2_OnPatch_2D_host, Patch3, Patch6, Patch9>
            derivs_ymax_369(derivs_ymax_3, derivs_ymax_6, derivs_ymax_9);



    MultipatchField<Deriv2_OnPatch_2D_host, Patch1, Patch2, Patch3>
            derivs_ymax_123(derivs_ymax_1, derivs_ymax_2, derivs_ymax_3);

    MultipatchField<Deriv2_OnPatch_2D_host, Patch4, Patch5, Patch6>
            derivs_ymax_456(derivs_ymax_4, derivs_ymax_5, derivs_ymax_6);

    MultipatchField<Deriv2_OnPatch_2D_host, Patch7, Patch8, Patch9>
            derivs_ymax_789(derivs_ymax_7, derivs_ymax_8, derivs_ymax_9);

    MultipatchField<Deriv2_OnPatch_2D_host, Patch1, Patch2, Patch3>
            derivs_ymin_123(derivs_ymin_1, derivs_ymin_2, derivs_ymin_3);

    MultipatchField<Deriv2_OnPatch_2D_host, Patch4, Patch5, Patch6>
            derivs_ymin_456(derivs_ymin_4, derivs_ymin_5, derivs_ymin_6);

    MultipatchField<Deriv2_OnPatch_2D_host, Patch7, Patch8, Patch9>
            derivs_ymin_789(derivs_ymin_7, derivs_ymin_8, derivs_ymin_9);


    // Instantiate cross-derivatives of the test functions ========================================
    IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>
            derivs_xy_idx_range(deriv_x_idx_range, deriv_y_idx_range);

    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_1_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_2_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_3_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_4_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_5_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_6_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_7_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_8_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_9_alloc(
            derivs_xy_idx_range);

    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_1(
            get_field(derivs_xy_min_min_1_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_2(
            get_field(derivs_xy_min_min_2_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_3(
            get_field(derivs_xy_min_min_3_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_4(
            get_field(derivs_xy_min_min_4_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_5(
            get_field(derivs_xy_min_min_5_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_6(
            get_field(derivs_xy_min_min_6_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_7(
            get_field(derivs_xy_min_min_7_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_8(
            get_field(derivs_xy_min_min_8_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_min_9(
            get_field(derivs_xy_min_min_9_alloc));


    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_1_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_2_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_3_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_4_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_5_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_6_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_7_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_8_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_9_alloc(
            derivs_xy_idx_range);

    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_1(
            get_field(derivs_xy_min_max_1_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_2(
            get_field(derivs_xy_min_max_2_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_3(
            get_field(derivs_xy_min_max_3_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_4(
            get_field(derivs_xy_min_max_4_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_5(
            get_field(derivs_xy_min_max_5_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_6(
            get_field(derivs_xy_min_max_6_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_7(
            get_field(derivs_xy_min_max_7_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_8(
            get_field(derivs_xy_min_max_8_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_min_max_9(
            get_field(derivs_xy_min_max_9_alloc));


    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_1_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_2_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_3_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_4_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_5_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_6_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_7_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_8_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_9_alloc(
            derivs_xy_idx_range);

    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_1(
            get_field(derivs_xy_max_min_1_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_2(
            get_field(derivs_xy_max_min_2_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_3(
            get_field(derivs_xy_max_min_3_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_4(
            get_field(derivs_xy_max_min_4_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_5(
            get_field(derivs_xy_max_min_5_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_6(
            get_field(derivs_xy_max_min_6_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_7(
            get_field(derivs_xy_max_min_7_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_8(
            get_field(derivs_xy_max_min_8_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_min_9(
            get_field(derivs_xy_max_min_9_alloc));


    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_1_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_2_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_3_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_4_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_5_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_6_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_7_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_8_alloc(
            derivs_xy_idx_range);
    host_t<DFieldMem<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_9_alloc(
            derivs_xy_idx_range);

    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_1(
            get_field(derivs_xy_max_max_1_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_2(
            get_field(derivs_xy_max_max_2_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_3(
            get_field(derivs_xy_max_max_3_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_4(
            get_field(derivs_xy_max_max_4_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_5(
            get_field(derivs_xy_max_max_5_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_6(
            get_field(derivs_xy_max_max_6_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_7(
            get_field(derivs_xy_max_max_7_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_8(
            get_field(derivs_xy_max_max_8_alloc));
    host_t<DField<IdxRange<ddc::Deriv<X>, ddc::Deriv<Y>>>> derivs_xy_max_max_9(
            get_field(derivs_xy_max_max_9_alloc));

    // Collect the cross-derivatives.
    std::tuple derivs_xy_min_min(
            derivs_xy_min_min_1,
            derivs_xy_min_min_2,
            derivs_xy_min_min_3,
            derivs_xy_min_min_4,
            derivs_xy_min_min_5,
            derivs_xy_min_min_6,
            derivs_xy_min_min_7,
            derivs_xy_min_min_8,
            derivs_xy_min_min_9);

    std::tuple derivs_xy_min_max(
            derivs_xy_min_max_1,
            derivs_xy_min_max_2,
            derivs_xy_min_max_3,
            derivs_xy_min_max_4,
            derivs_xy_min_max_5,
            derivs_xy_min_max_6,
            derivs_xy_min_max_7,
            derivs_xy_min_max_8,
            derivs_xy_min_max_9);

    std::tuple derivs_xy_max_min(
            derivs_xy_max_min_1,
            derivs_xy_max_min_2,
            derivs_xy_max_min_3,
            derivs_xy_max_min_4,
            derivs_xy_max_min_5,
            derivs_xy_max_min_6,
            derivs_xy_max_min_7,
            derivs_xy_max_min_8,
            derivs_xy_max_min_9);

    std::tuple derivs_xy_max_max(
            derivs_xy_max_max_1,
            derivs_xy_max_max_2,
            derivs_xy_max_max_3,
            derivs_xy_max_max_4,
            derivs_xy_max_max_5,
            derivs_xy_max_max_6,
            derivs_xy_max_max_7,
            derivs_xy_max_max_8,
            derivs_xy_max_max_9);

    std::tuple derivs_xy_min_max_123(derivs_xy_min_max_1, derivs_xy_min_max_2, derivs_xy_min_max_3);
    std::tuple derivs_xy_max_max_123(derivs_xy_max_max_1, derivs_xy_max_max_2, derivs_xy_max_max_3);

    std::tuple derivs_xy_min_max_456(derivs_xy_min_max_4, derivs_xy_min_max_5, derivs_xy_min_max_6);
    std::tuple derivs_xy_max_max_456(derivs_xy_max_max_4, derivs_xy_max_max_5, derivs_xy_max_max_6);

    std::tuple derivs_xy_min_max_789(derivs_xy_min_max_7, derivs_xy_min_max_8, derivs_xy_min_max_9);
    std::tuple derivs_xy_max_max_789(derivs_xy_max_max_7, derivs_xy_max_max_8, derivs_xy_max_max_9);


    std::tuple derivs_xy_min_min_123(derivs_xy_min_min_1, derivs_xy_min_min_2, derivs_xy_min_min_3);
    std::tuple derivs_xy_max_min_123(derivs_xy_max_min_1, derivs_xy_max_min_2, derivs_xy_max_min_3);

    std::tuple derivs_xy_min_min_456(derivs_xy_min_min_4, derivs_xy_min_min_5, derivs_xy_min_min_6);
    std::tuple derivs_xy_max_min_456(derivs_xy_max_min_4, derivs_xy_max_min_5, derivs_xy_max_min_6);

    std::tuple derivs_xy_min_min_789(derivs_xy_min_min_7, derivs_xy_min_min_8, derivs_xy_min_min_9);
    std::tuple derivs_xy_max_min_789(derivs_xy_max_min_7, derivs_xy_max_min_8, derivs_xy_max_min_9);


    // Initialise the data =======================================================================

    initialise_2D_function<GridX<1>, GridY<1>>(function_1);
    initialise_2D_function<GridX<2>, GridY<2>>(function_2);
    initialise_2D_function<GridX<3>, GridY<3>>(function_3);
    initialise_2D_function<GridX<4>, GridY<4>>(function_4);
    initialise_2D_function<GridX<5>, GridY<5>>(function_5);
    initialise_2D_function<GridX<6>, GridY<6>>(function_6);
    initialise_2D_function<GridX<7>, GridY<7>>(function_7);
    initialise_2D_function<GridX<8>, GridY<8>>(function_8);
    initialise_2D_function<GridX<9>, GridY<9>>(function_9);
    initialise_2D_function<GridXg, GridYg>(function_g);

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


    // Initiliase the boundary derivatives.
    initialise_2D_y_derivative<GridX<7>, GridY<7>>(
            derivs_ymin_7,
            idx_range_y7.front(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_y_derivative<GridX<8>, GridY<8>>(
            derivs_ymin_8,
            idx_range_y8.front(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_y_derivative<GridX<9>, GridY<9>>(
            derivs_ymin_9,
            idx_range_y9.front(),
            evaluator_g,
            get_const_field(function_g_coef));

    initialise_2D_y_derivative<GridX<1>, GridY<1>>(
            derivs_ymax_1,
            idx_range_y1.back(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_y_derivative<GridX<2>, GridY<2>>(
            derivs_ymax_2,
            idx_range_y2.back(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_y_derivative<GridX<3>, GridY<3>>(
            derivs_ymax_3,
            idx_range_y3.back(),
            evaluator_g,
            get_const_field(function_g_coef));


    initialise_2D_xy_derivative<GridX<1>, GridY<1>>(
            derivs_xy_min_max_1,
            idx_range_x1.front(),
            idx_range_y1.back(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_xy_derivative<GridX<1>, GridY<1>>(
            derivs_xy_max_max_1,
            idx_range_x1.back(),
            idx_range_y1.back(),
            evaluator_g,
            get_const_field(function_g_coef));

    initialise_2D_xy_derivative<GridX<2>, GridY<2>>(
            derivs_xy_min_max_2,
            idx_range_x2.front(),
            idx_range_y2.back(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_xy_derivative<GridX<2>, GridY<2>>(
            derivs_xy_max_max_2,
            idx_range_x2.back(),
            idx_range_y2.back(),
            evaluator_g,
            get_const_field(function_g_coef));

    initialise_2D_xy_derivative<GridX<3>, GridY<3>>(
            derivs_xy_min_max_3,
            idx_range_x3.front(),
            idx_range_y3.back(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_xy_derivative<GridX<3>, GridY<3>>(
            derivs_xy_max_max_3,
            idx_range_x3.back(),
            idx_range_y3.back(),
            evaluator_g,
            get_const_field(function_g_coef));

    initialise_2D_xy_derivative<GridX<7>, GridY<7>>(
            derivs_xy_min_min_7,
            idx_range_x7.front(),
            idx_range_y7.front(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_xy_derivative<GridX<7>, GridY<7>>(
            derivs_xy_max_min_7,
            idx_range_x7.back(),
            idx_range_y7.front(),
            evaluator_g,
            get_const_field(function_g_coef));

    initialise_2D_xy_derivative<GridX<8>, GridY<8>>(
            derivs_xy_min_min_8,
            idx_range_x8.front(),
            idx_range_y8.front(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_xy_derivative<GridX<8>, GridY<8>>(
            derivs_xy_max_min_8,
            idx_range_x8.back(),
            idx_range_y8.front(),
            evaluator_g,
            get_const_field(function_g_coef));

    initialise_2D_xy_derivative<GridX<9>, GridY<9>>(
            derivs_xy_min_min_9,
            idx_range_x9.front(),
            idx_range_y9.front(),
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_2D_xy_derivative<GridX<9>, GridY<9>>(
            derivs_xy_max_min_9,
            idx_range_x9.back(),
            idx_range_y9.front(),
            evaluator_g,
            get_const_field(function_g_coef));



    // Solve each matrix system ==================================================================
    // TODO: modify order: call derivatives first then functions.
    matrix_123.solve(functions_123, derivs_xmin_123, derivs_xmax_123);
    matrix_456.solve(functions_456, derivs_xmin_456, derivs_xmax_456);
    matrix_789.solve(functions_789, derivs_xmin_789, derivs_xmax_789);

    matrix_147.solve(functions_147, derivs_ymin_147, derivs_ymax_147);
    matrix_258.solve(functions_258, derivs_ymin_258, derivs_ymax_258);
    matrix_369.solve(functions_369, derivs_ymin_369, derivs_ymax_369);

    matrix_456.solve(derivs_ymax_456, derivs_xy_min_max_456, derivs_xy_max_max_456);
    matrix_789.solve(derivs_ymax_789, derivs_xy_min_max_789, derivs_xy_max_max_789);

    matrix_123.solve(derivs_ymin_123, derivs_xy_min_min_123, derivs_xy_max_min_123);
    matrix_456.solve(derivs_ymin_456, derivs_xy_min_min_456, derivs_xy_max_min_456);

    // Test the values of the derivatives ========================================================
    // // Build global spline representation ---
    // SplineRThetagBuilder builder_g(idx_range_xy_g);

    // host_t<DFieldMem<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef_alloc(
    //         builder_g.batched_spline_domain(idx_range_xy_g));
    // host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef
    //         = get_field(function_g_coef_alloc);

    // builder_g(function_g_coef, get_const_field(function_g));

    // // Global spline evaluator ---
    // ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymin_g(yg_min);
    // ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymax_g(yg_max);
    // ddc::PeriodicExtrapolationRule<Xg> bc_x_g;
    // SplineRThetagEvaluator evaluator_g(bc_x_g, bc_x_g, bc_ymin_g, bc_ymax_g);

    // Check only the derivatives values ---
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

    check_all_x_derivatives<Patch1, Patch2, Patch3, Patch4, Patch5, Patch6, Patch7, Patch8, Patch9>(
            derivs_xmin,
            derivs_xmax,
            evaluator_g,
            get_const_field(function_g_coef),
            idx_ranges);

    check_all_y_derivatives<Patch1, Patch2, Patch3, Patch4, Patch5, Patch6, Patch7, Patch8, Patch9>(
            derivs_ymin,
            derivs_ymax,
            evaluator_g,
            get_const_field(function_g_coef),
            idx_ranges);

    check_all_xy_derivatives<
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>(
            derivs_xy_min_min,
            derivs_xy_max_min,
            derivs_xy_min_max,
            derivs_xy_max_max,
            evaluator_g,
            get_const_field(function_g_coef),
            idx_ranges,
            std::make_integer_sequence<std::size_t, 9> {});
}
