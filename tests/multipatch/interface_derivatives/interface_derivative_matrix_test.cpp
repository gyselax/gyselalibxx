// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "9patches_2d_periodic_strips_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "interface.hpp"
#include "interface_derivative_matrix.hpp"
#include "interface_derivatives_test_utils.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"


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

template <std::size_t Index>
using PatchI = ddc::type_seq_element_t<Index, PatchOrdering>;

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


// USEFUL FUNCTIONS ==============================================================================
/// @brief Convert a local coordinate into a global coordinate.
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

/// @brief Initialise all the functions defined on the patches.
template <class... Patches>
void initialise_all_functions(MultipatchField<DFieldOnPatch_host, Patches...> const& functions)
{
    (initialise_2D_function<typename Patches::Grid1, typename Patches::Grid2>(
             functions.template get<Patches>()),
     ...);
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
                Idx<Grid1, Grid2> idx_xy(Idx<Grid1>(idx), idx_y);
                Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx_xy)));
                deriv_y(idx) = evaluator_g.deriv_dim_2(interface_coord, function_g_coef);
            });
}

/// @brief Initialise all the y-derivatives on the lower/left bounds.
template <class... Patches>
void initialise_all_y_derivatives_min(
        MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& derivs_ymin,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        SplineRThetagEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
{
    (initialise_2D_y_derivative<typename Patches::Grid1, typename Patches::Grid2>(
             derivs_ymin.template get<Patches>(),
             typename Patches::IdxRange2(idx_ranges.template get<Patches>()).front(),
             evaluator_g,
             function_g_coef),
     ...);
}

/// @brief Initialise all the y-derivatives on the upper/right bounds.
template <class... Patches>
void initialise_all_y_derivatives_max(
        MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& derivs_ymax,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        SplineRThetagEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
{
    (initialise_2D_y_derivative<typename Patches::Grid1, typename Patches::Grid2>(
             derivs_ymax.template get<Patches>(),
             typename Patches::IdxRange2(idx_ranges.template get<Patches>()).back(),
             evaluator_g,
             function_g_coef),
     ...);
}


/// @brief Initialise the cross-derivatives from the global spline.
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
    Idx<Grid1, Grid2> idx(idx_x, idx_y);
    Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx)));
    deriv_xy(get_idx_range(deriv_xy).front())
            = evaluator_g.deriv_1_and_2(interface_coord, function_g_coef);
}

/// @brief Initialise all the cross-derivatives on the lower/left - lower/left corners.
template <std::size_t... I, class... Patches>
void initialise_all_xy_derivatives_min_min(
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_min_min,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        SplineRThetagEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::integer_sequence<std::size_t, I...>)
{
    (initialise_2D_xy_derivative<typename Patches::Grid1, typename Patches::Grid2>(
             std::get<I>(derivs_min_min),
             typename Patches::IdxRange1(idx_ranges.template get<Patches>()).front(),
             typename Patches::IdxRange2(idx_ranges.template get<Patches>()).front(),
             evaluator_g,
             function_g_coef),
     ...);
}

/// @brief Initialise all the cross-derivatives on the upper/right - lower/left corners.
template <std::size_t... I, class... Patches>
void initialise_all_xy_derivatives_max_min(
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_max_min,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        SplineRThetagEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::integer_sequence<std::size_t, I...>)
{
    (initialise_2D_xy_derivative<typename Patches::Grid1, typename Patches::Grid2>(
             std::get<I>(derivs_max_min),
             typename Patches::IdxRange1(idx_ranges.template get<Patches>()).back(),
             typename Patches::IdxRange2(idx_ranges.template get<Patches>()).front(),
             evaluator_g,
             function_g_coef),
     ...);
}

/// @brief Initialise all the cross-derivatives on the lower/left - upper/right corners.
template <std::size_t... I, class... Patches>
void initialise_all_xy_derivatives_min_max(
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_min_max,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        SplineRThetagEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::integer_sequence<std::size_t, I...>)
{
    (initialise_2D_xy_derivative<typename Patches::Grid1, typename Patches::Grid2>(
             std::get<I>(derivs_min_max),
             typename Patches::IdxRange1(idx_ranges.template get<Patches>()).front(),
             typename Patches::IdxRange2(idx_ranges.template get<Patches>()).back(),
             evaluator_g,
             function_g_coef),
     ...);
}

/// @brief Initialise all the cross-derivatives on the upper/right - lower/left corners.
template <std::size_t... I, class... Patches>
void initialise_all_xy_derivatives_max_max(
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_max_max,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        SplineRThetagEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::integer_sequence<std::size_t, I...>)
{
    (initialise_2D_xy_derivative<typename Patches::Grid1, typename Patches::Grid2>(
             std::get<I>(derivs_max_max),
             typename Patches::IdxRange1(idx_ranges.template get<Patches>()).back(),
             typename Patches::IdxRange2(idx_ranges.template get<Patches>()).back(),
             evaluator_g,
             function_g_coef),
     ...);
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



    // TEST OPERATORS ============================================================================
    /** @brief Check that the local grids and the equivalent global grid match together for a 
     * given patch. The integers x_shift and y_shift correspond to a shift in the indices to match
     * with the correct index on the global grid. 
     */
    template <class Patch>
    void check_interpolation_grids(
            typename Patch::IdxRange12 const& idx_range,
            int const x_shift,
            int const y_shift)
    {
        ddc::for_each(idx_range, [&](typename Patch::Idx12 const& idx) {
            typename Patch::IdxStep1 idx_x(
                    typename Patch::Idx1(idx) - typename Patch::IdxRange1(idx_range).front());
            typename Patch::IdxStep2 idx_y(
                    typename Patch::Idx2(idx) - typename Patch::IdxRange2(idx_range).front());
            Idx<GridXg, GridYg> idx_g(idx_x.value() + x_shift, idx_y.value() + y_shift);
            EXPECT_NEAR(
                    ddc::coordinate(typename Patch::Idx1(idx)),
                    ddc::coordinate(Idx<GridXg>(idx_g)),
                    1e-15);
            EXPECT_NEAR(
                    ddc::coordinate(typename Patch::Idx2(idx)),
                    ddc::coordinate(Idx<GridYg>(idx_g)),
                    1e-15);
        });
    }


    /** @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
     * the interfaces. 
     */
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

    /** @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
     * the interfaces for a given patch.
     */
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

            double const local_deriv = local_derivs(idx_deriv, idx_par);
            double const global_deriv = evaluator_g.deriv_dim_1(interface_coord, function_g_coef);

            EXPECT_NEAR(local_deriv, global_deriv, 2e-14);
        });
    }


    /** @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
     * the interfaces.
     */
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

    /** @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
     * the interfaces for a given patch.
     */
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

    /** @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
     * the interfaces.
     */
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

    /** @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
     * the interfaces for a given patch.
     */
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


    /** @brief Check agreement between the local splines defined with the computed derivatives 
     * and the global spline. 
     */
    template <class... Patches, std::size_t... I>
    void check_all_spline_representation_conformity(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            MultipatchField<DFieldOnPatch_host, Patches...> const& functions,
            MultipatchField<Deriv1_OnPatch_2D_host, Patches...> const& local_derivs_xmin,
            MultipatchField<Deriv1_OnPatch_2D_host, Patches...> const& local_derivs_xmax,
            MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& local_derivs_ymin,
            MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& local_derivs_ymax,
            std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& local_derivs_min_min,
            std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& local_derivs_max_min,
            std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& local_derivs_min_max,
            std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& local_derivs_max_max,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            std::integer_sequence<std::size_t, I...>)
    {
        (check_spline_representation_conformity<PatchI<I>, I>(
                 idx_ranges.template get<PatchI<I>>(),
                 functions.template get<PatchI<I>>(),
                 local_derivs_xmin.template get<PatchI<I>>(),
                 local_derivs_xmax.template get<PatchI<I>>(),
                 local_derivs_ymin.template get<PatchI<I>>(),
                 local_derivs_ymax.template get<PatchI<I>>(),
                 std::get<I>(local_derivs_min_min),
                 std::get<I>(local_derivs_max_min),
                 std::get<I>(local_derivs_min_max),
                 std::get<I>(local_derivs_max_max),
                 evaluator_g,
                 get_const_field(function_g_coef)),
         ...);
    };

    /** @brief Check agreement between the local splines defined with the computed derivatives 
     * and the global spline for a given patch.  
     */
    template <class Patch, std::size_t I>
    void check_spline_representation_conformity(
            typename Patch::IdxRange12 const& idx_range_xy,
            host_t<DFieldOnPatch<Patch>> const& function,
            Deriv1_OnPatch_2D_host<Patch> const& deriv_xmin,
            Deriv1_OnPatch_2D_host<Patch> const& deriv_xmax,
            Deriv2_OnPatch_2D_host<Patch> const& deriv_ymin,
            Deriv2_OnPatch_2D_host<Patch> const& deriv_ymax,
            Deriv12_OnPatch_2D_host<Patch> const& deriv_xy_min_min,
            Deriv12_OnPatch_2D_host<Patch> const& deriv_xy_max_min,
            Deriv12_OnPatch_2D_host<Patch> const& deriv_xy_min_max,
            Deriv12_OnPatch_2D_host<Patch> const& deriv_xy_max_max,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
    {
        // For Patch7, Patch8, Patch9, the local lower Y-boundary is with the outside.
        const ddc::BoundCond BoundCondY1
                = (I >= 6) ? ddc::BoundCond::GREVILLE : ddc::BoundCond::HERMITE;
        // For Patch1, Patch2, Patch3, the local upper Y-boundary is with the outside.
        const ddc::BoundCond BoundCondY2
                = (I <= 2) ? ddc::BoundCond::GREVILLE : ddc::BoundCond::HERMITE;

        // Build the local spline representation -------------------------------------------------
        ddc::SplineBuilder2D<
                HostExecSpace,
                typename HostExecSpace::memory_space,
                typename Patch::BSplines1,
                typename Patch::BSplines2,
                typename Patch::Grid1,
                typename Patch::Grid2,
                ddc::BoundCond::HERMITE,
                ddc::BoundCond::HERMITE,
                BoundCondY1,
                BoundCondY2,
                ddc::SplineSolver::LAPACK>
                builder(idx_range_xy);

        SplineCoeffMemOnPatch_2D_host<Patch> function_coef_alloc(
                builder.batched_spline_domain(idx_range_xy));
        SplineCoeffOnPatch_2D_host<Patch> function_coef = get_field(function_coef_alloc);

        // If the boundary is not a ddc::BoundCond::HERMITE, we don't use derivatives.
        if constexpr (
                (BoundCondY1 == ddc::BoundCond::HERMITE)
                && (BoundCondY2 == ddc::BoundCond::HERMITE)) {
            builder(function_coef,
                    get_const_field(function),
                    std::optional(get_const_field(deriv_xmin)),
                    std::optional(get_const_field(deriv_xmax)),
                    std::optional(get_const_field(deriv_ymin)),
                    std::optional(get_const_field(deriv_ymax)),
                    std::optional(get_const_field(deriv_xy_min_min)),
                    std::optional(get_const_field(deriv_xy_max_min)),
                    std::optional(get_const_field(deriv_xy_min_max)),
                    std::optional(get_const_field(deriv_xy_max_max)));
        } else if constexpr (
                (BoundCondY1 == ddc::BoundCond::HERMITE)
                && !(BoundCondY2 == ddc::BoundCond::HERMITE)) {
            builder(function_coef,
                    get_const_field(function),
                    std::optional(get_const_field(deriv_xmin)),
                    std::optional(get_const_field(deriv_xmax)),
                    std::optional(get_const_field(deriv_ymin)),
                    std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional(get_const_field(deriv_xy_min_min)),
                    std::optional(get_const_field(deriv_xy_max_min)),
                    std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt});
        } else if constexpr (
                !(BoundCondY1 == ddc::BoundCond::HERMITE)
                && (BoundCondY2 == ddc::BoundCond::HERMITE)) {
            builder(function_coef,
                    get_const_field(function),
                    std::optional(get_const_field(deriv_xmin)),
                    std::optional(get_const_field(deriv_xmax)),
                    std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional(get_const_field(deriv_ymax)),
                    std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional(get_const_field(deriv_xy_min_max)),
                    std::optional(get_const_field(deriv_xy_max_max)));
        } else {
            builder(function_coef,
                    get_const_field(function),
                    std::optional(get_const_field(deriv_xmin)),
                    std::optional(get_const_field(deriv_xmax)),
                    std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                    std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt});
        }

        // Define local spline evaluator ---------------------------------------------------------
        Coord<X> const x_min(ddc::discrete_space<typename Patch::BSplines1>().rmin());
        Coord<X> const x_max(ddc::discrete_space<typename Patch::BSplines1>().rmax());
        Coord<Y> const y_min(ddc::discrete_space<typename Patch::BSplines2>().rmin());
        Coord<Y> const y_max(ddc::discrete_space<typename Patch::BSplines2>().rmax());
        ddc::ConstantExtrapolationRule<X, Y> bc_xmin(x_min, y_min, y_max);
        ddc::ConstantExtrapolationRule<X, Y> bc_xmax(x_max, y_min, y_max);
        ddc::ConstantExtrapolationRule<Y, X> bc_ymin(y_min, x_min, x_max);
        ddc::ConstantExtrapolationRule<Y, X> bc_ymax(y_max, x_min, x_max);
        ddc::SplineEvaluator2D<
                HostExecSpace,
                typename HostExecSpace::memory_space,
                typename Patch::BSplines1,
                typename Patch::BSplines2,
                typename Patch::Grid1,
                typename Patch::Grid2,
                ddc::ConstantExtrapolationRule<X, Y>,
                ddc::ConstantExtrapolationRule<X, Y>,
                ddc::ConstantExtrapolationRule<Y, X>,
                ddc::ConstantExtrapolationRule<Y, X>>
                evaluator(bc_xmin, bc_xmax, bc_ymin, bc_ymax);

        // Define evaluation points at the centre of the cells -----------------------------------
        host_t<CoordFieldMemOnPatch<Patch>> eval_points_alloc(idx_range_xy);
        host_t<CoordFieldOnPatch<Patch>> eval_points(eval_points_alloc);

        ddc::for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
            Coord<X, Y> const mesh_point(ddc::coordinate(idx));
            typename Patch::Idx1 idx_1(idx);
            typename Patch::Idx2 idx_2(idx);
            Coord<X> dx = (idx_1 != typename Patch::IdxRange1(idx_range_xy).back())
                          * distance_at_right(idx_1);
            Coord<Y> dy = (idx_2 != typename Patch::IdxRange2(idx_range_xy).back())
                          * distance_at_right(idx_2);
            eval_points(idx) = mesh_point + Coord<X, Y>(dx, dy);
        });

        // Evaluate and compare the local and global spline representations ----------------------
        ddc::for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
            Coord<Xg, Yg> const eval_point_g(get_global_coord(eval_points(idx)));

            double const local_spline = evaluator(eval_points(idx), get_const_field(function_coef));
            double const global_spline
                    = evaluator_g(eval_point_g, get_const_field(function_g_coef));

            EXPECT_NEAR(local_spline, global_spline, 1e-14);
        });
    }
};

} // end namespace



// Check that the local grids and the equivalent global grid match together.
TEST_F(InterfaceDerivativeMatrixTest, InterpolationPointsCheck)
{
    int const x_shift1 = x1_ncells.value();
    int const x_shift2 = x1_ncells.value() + x2_ncells.value();

    int const y_shift1 = y4_ncells.value() + 1;
    int const y_shift2 = y7_ncells.value() + y4_ncells.value() + 1;

    check_interpolation_grids<Patch1>(idx_range_xy1, 0, y_shift2);
    check_interpolation_grids<Patch2>(idx_range_xy2, x_shift1, y_shift2);
    check_interpolation_grids<Patch3>(idx_range_xy3, x_shift2, y_shift2);
    check_interpolation_grids<Patch4>(idx_range_xy4, 0, y_shift1);
    check_interpolation_grids<Patch5>(idx_range_xy5, x_shift1, y_shift1);
    check_interpolation_grids<Patch6>(idx_range_xy6, x_shift2, y_shift1);
    check_interpolation_grids<Patch7>(idx_range_xy7, 0, 0);
    check_interpolation_grids<Patch8>(idx_range_xy8, x_shift1, 0);
    check_interpolation_grids<Patch9>(idx_range_xy9, x_shift2, 0);
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


    // TODO: develop a MultipatchType for SingleInterfaceDerivativesCalculator.
    // MultipatchType<SingleInterfaceDerivativesCalculatorOnInterface, Interface_1_2, Interface_2_3, Interface_3_1>
    //         derivative_calculators_123_(
    //                 derivatives_calculator_1_2,
    //                 derivatives_calculator_2_3,
    //                 derivatives_calculator_3_1);

    MultipatchType<IdxRangeOnPatch, Patch1, Patch2, Patch3>
            idx_ranges_123(idx_range_xy1, idx_range_xy2, idx_range_xy3);
    MultipatchType<IdxRangeOnPatch, Patch4, Patch5, Patch6>
            idx_ranges_456(idx_range_xy4, idx_range_xy5, idx_range_xy6);
    MultipatchType<IdxRangeOnPatch, Patch7, Patch8, Patch9>
            idx_ranges_789(idx_range_xy7, idx_range_xy8, idx_range_xy9);

    // For 1|4|7, we artificially add another patch to test if the method will modify the correct patches.
    MultipatchType<IdxRangeOnPatch, Patch1, Patch4, Patch7, Patch2>
            idx_ranges_147(idx_range_xy1, idx_range_xy4, idx_range_xy7, idx_range_xy2);
    MultipatchType<IdxRangeOnPatch, Patch2, Patch5, Patch8>
            idx_ranges_258(idx_range_xy2, idx_range_xy5, idx_range_xy8);
    MultipatchType<IdxRangeOnPatch, Patch3, Patch6, Patch9>
            idx_ranges_369(idx_range_xy3, idx_range_xy6, idx_range_xy9);


    // Instantiate the matrix calculators --------------------------------------------------------
    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<1>,
            DConstFieldOnPatch_host,
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
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            Kokkos::DefaultHostExecutionSpace,
            Patch7,
            Patch8,
            Patch9>
            matrix_789(idx_ranges_789, derivative_calculators_789);


    // Test with an exact patch to check it will only take the needed patches.
    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<1>,
            DConstFieldOnPatch_host,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            Kokkos::DefaultHostExecutionSpace,
            Patch1,
            Patch4,
            Patch7,
            Patch2>
            matrix_147(idx_ranges_147, derivative_calculators_147);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<2>,
            DConstFieldOnPatch_host,
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
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            Kokkos::DefaultHostExecutionSpace,
            Patch3,
            Patch6,
            Patch9>
            matrix_369(idx_ranges_369, derivative_calculators_369);

    // Instantiate test function values ==========================================================
    host_t<DFieldMemOnPatch<Patch1>> function_1_alloc(idx_range_xy1);
    host_t<DFieldMemOnPatch<Patch2>> function_2_alloc(idx_range_xy2);
    host_t<DFieldMemOnPatch<Patch3>> function_3_alloc(idx_range_xy3);
    host_t<DFieldMemOnPatch<Patch4>> function_4_alloc(idx_range_xy4);
    host_t<DFieldMemOnPatch<Patch5>> function_5_alloc(idx_range_xy5);
    host_t<DFieldMemOnPatch<Patch6>> function_6_alloc(idx_range_xy6);
    host_t<DFieldMemOnPatch<Patch7>> function_7_alloc(idx_range_xy7);
    host_t<DFieldMemOnPatch<Patch8>> function_8_alloc(idx_range_xy8);
    host_t<DFieldMemOnPatch<Patch9>> function_9_alloc(idx_range_xy9);

    host_t<DFieldOnPatch<Patch1>> function_1(get_field(function_1_alloc));
    host_t<DFieldOnPatch<Patch2>> function_2(get_field(function_2_alloc));
    host_t<DFieldOnPatch<Patch3>> function_3(get_field(function_3_alloc));
    host_t<DFieldOnPatch<Patch4>> function_4(get_field(function_4_alloc));
    host_t<DFieldOnPatch<Patch5>> function_5(get_field(function_5_alloc));
    host_t<DFieldOnPatch<Patch6>> function_6(get_field(function_6_alloc));
    host_t<DFieldOnPatch<Patch7>> function_7(get_field(function_7_alloc));
    host_t<DFieldOnPatch<Patch8>> function_8(get_field(function_8_alloc));
    host_t<DFieldOnPatch<Patch9>> function_9(get_field(function_9_alloc));

    host_t<DFieldMem<IdxRange<GridXg, GridYg>>> function_g_alloc(idx_range_xy_g);
    host_t<DField<IdxRange<GridXg, GridYg>>> function_g = get_field(function_g_alloc);

    // Collect the function values.
    MultipatchField<
            DFieldOnPatch_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            functions(
                    function_1,
                    function_2,
                    function_3,
                    function_4,
                    function_5,
                    function_6,
                    function_7,
                    function_8,
                    function_9);

    MultipatchField<DFieldOnPatch_host, Patch1, Patch2, Patch3>
            functions_123(function_1, function_2, function_3);
    MultipatchField<DFieldOnPatch_host, Patch4, Patch5, Patch6>
            functions_456(function_4, function_5, function_6);
    MultipatchField<DFieldOnPatch_host, Patch7, Patch8, Patch9>
            functions_789(function_7, function_8, function_9);

    MultipatchField<DFieldOnPatch_host, Patch1, Patch4, Patch7, Patch2>
            functions_147(function_1, function_4, function_7, function_2);
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

    DerivMem1_OnPatch_2D_host<Patch1> derivs_xmin_1_alloc(derivs_x_idx_range_1);
    DerivMem1_OnPatch_2D_host<Patch2> derivs_xmin_2_alloc(derivs_x_idx_range_2);
    DerivMem1_OnPatch_2D_host<Patch3> derivs_xmin_3_alloc(derivs_x_idx_range_3);
    DerivMem1_OnPatch_2D_host<Patch4> derivs_xmin_4_alloc(derivs_x_idx_range_4);
    DerivMem1_OnPatch_2D_host<Patch5> derivs_xmin_5_alloc(derivs_x_idx_range_5);
    DerivMem1_OnPatch_2D_host<Patch6> derivs_xmin_6_alloc(derivs_x_idx_range_6);
    DerivMem1_OnPatch_2D_host<Patch7> derivs_xmin_7_alloc(derivs_x_idx_range_7);
    DerivMem1_OnPatch_2D_host<Patch8> derivs_xmin_8_alloc(derivs_x_idx_range_8);
    DerivMem1_OnPatch_2D_host<Patch9> derivs_xmin_9_alloc(derivs_x_idx_range_9);

    Deriv1_OnPatch_2D_host<Patch1> derivs_xmin_1 = get_field(derivs_xmin_1_alloc);
    Deriv1_OnPatch_2D_host<Patch2> derivs_xmin_2 = get_field(derivs_xmin_2_alloc);
    Deriv1_OnPatch_2D_host<Patch3> derivs_xmin_3 = get_field(derivs_xmin_3_alloc);
    Deriv1_OnPatch_2D_host<Patch4> derivs_xmin_4 = get_field(derivs_xmin_4_alloc);
    Deriv1_OnPatch_2D_host<Patch5> derivs_xmin_5 = get_field(derivs_xmin_5_alloc);
    Deriv1_OnPatch_2D_host<Patch6> derivs_xmin_6 = get_field(derivs_xmin_6_alloc);
    Deriv1_OnPatch_2D_host<Patch7> derivs_xmin_7 = get_field(derivs_xmin_7_alloc);
    Deriv1_OnPatch_2D_host<Patch8> derivs_xmin_8 = get_field(derivs_xmin_8_alloc);
    Deriv1_OnPatch_2D_host<Patch9> derivs_xmin_9 = get_field(derivs_xmin_9_alloc);


    DerivMem1_OnPatch_2D_host<Patch1> derivs_xmax_1_alloc(derivs_x_idx_range_1);
    DerivMem1_OnPatch_2D_host<Patch2> derivs_xmax_2_alloc(derivs_x_idx_range_2);
    DerivMem1_OnPatch_2D_host<Patch3> derivs_xmax_3_alloc(derivs_x_idx_range_3);
    DerivMem1_OnPatch_2D_host<Patch4> derivs_xmax_4_alloc(derivs_x_idx_range_4);
    DerivMem1_OnPatch_2D_host<Patch5> derivs_xmax_5_alloc(derivs_x_idx_range_5);
    DerivMem1_OnPatch_2D_host<Patch6> derivs_xmax_6_alloc(derivs_x_idx_range_6);
    DerivMem1_OnPatch_2D_host<Patch7> derivs_xmax_7_alloc(derivs_x_idx_range_7);
    DerivMem1_OnPatch_2D_host<Patch8> derivs_xmax_8_alloc(derivs_x_idx_range_8);
    DerivMem1_OnPatch_2D_host<Patch9> derivs_xmax_9_alloc(derivs_x_idx_range_9);

    Deriv1_OnPatch_2D_host<Patch1> derivs_xmax_1 = get_field(derivs_xmax_1_alloc);
    Deriv1_OnPatch_2D_host<Patch2> derivs_xmax_2 = get_field(derivs_xmax_2_alloc);
    Deriv1_OnPatch_2D_host<Patch3> derivs_xmax_3 = get_field(derivs_xmax_3_alloc);
    Deriv1_OnPatch_2D_host<Patch4> derivs_xmax_4 = get_field(derivs_xmax_4_alloc);
    Deriv1_OnPatch_2D_host<Patch5> derivs_xmax_5 = get_field(derivs_xmax_5_alloc);
    Deriv1_OnPatch_2D_host<Patch6> derivs_xmax_6 = get_field(derivs_xmax_6_alloc);
    Deriv1_OnPatch_2D_host<Patch7> derivs_xmax_7 = get_field(derivs_xmax_7_alloc);
    Deriv1_OnPatch_2D_host<Patch8> derivs_xmax_8 = get_field(derivs_xmax_8_alloc);
    Deriv1_OnPatch_2D_host<Patch9> derivs_xmax_9 = get_field(derivs_xmax_9_alloc);


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

    DerivMem2_OnPatch_2D_host<Patch1> derivs_ymin_1_alloc(derivs_y_idx_range_1);
    DerivMem2_OnPatch_2D_host<Patch2> derivs_ymin_2_alloc(derivs_y_idx_range_2);
    DerivMem2_OnPatch_2D_host<Patch3> derivs_ymin_3_alloc(derivs_y_idx_range_3);
    DerivMem2_OnPatch_2D_host<Patch4> derivs_ymin_4_alloc(derivs_y_idx_range_4);
    DerivMem2_OnPatch_2D_host<Patch5> derivs_ymin_5_alloc(derivs_y_idx_range_5);
    DerivMem2_OnPatch_2D_host<Patch6> derivs_ymin_6_alloc(derivs_y_idx_range_6);
    DerivMem2_OnPatch_2D_host<Patch7> derivs_ymin_7_alloc(derivs_y_idx_range_7);
    DerivMem2_OnPatch_2D_host<Patch8> derivs_ymin_8_alloc(derivs_y_idx_range_8);
    DerivMem2_OnPatch_2D_host<Patch9> derivs_ymin_9_alloc(derivs_y_idx_range_9);

    Deriv2_OnPatch_2D_host<Patch1> derivs_ymin_1 = get_field(derivs_ymin_1_alloc);
    Deriv2_OnPatch_2D_host<Patch2> derivs_ymin_2 = get_field(derivs_ymin_2_alloc);
    Deriv2_OnPatch_2D_host<Patch3> derivs_ymin_3 = get_field(derivs_ymin_3_alloc);
    Deriv2_OnPatch_2D_host<Patch4> derivs_ymin_4 = get_field(derivs_ymin_4_alloc);
    Deriv2_OnPatch_2D_host<Patch5> derivs_ymin_5 = get_field(derivs_ymin_5_alloc);
    Deriv2_OnPatch_2D_host<Patch6> derivs_ymin_6 = get_field(derivs_ymin_6_alloc);
    Deriv2_OnPatch_2D_host<Patch7> derivs_ymin_7 = get_field(derivs_ymin_7_alloc);
    Deriv2_OnPatch_2D_host<Patch8> derivs_ymin_8 = get_field(derivs_ymin_8_alloc);
    Deriv2_OnPatch_2D_host<Patch9> derivs_ymin_9 = get_field(derivs_ymin_9_alloc);


    DerivMem2_OnPatch_2D_host<Patch1> derivs_ymax_1_alloc(derivs_y_idx_range_1);
    DerivMem2_OnPatch_2D_host<Patch2> derivs_ymax_2_alloc(derivs_y_idx_range_2);
    DerivMem2_OnPatch_2D_host<Patch3> derivs_ymax_3_alloc(derivs_y_idx_range_3);
    DerivMem2_OnPatch_2D_host<Patch4> derivs_ymax_4_alloc(derivs_y_idx_range_4);
    DerivMem2_OnPatch_2D_host<Patch5> derivs_ymax_5_alloc(derivs_y_idx_range_5);
    DerivMem2_OnPatch_2D_host<Patch6> derivs_ymax_6_alloc(derivs_y_idx_range_6);
    DerivMem2_OnPatch_2D_host<Patch7> derivs_ymax_7_alloc(derivs_y_idx_range_7);
    DerivMem2_OnPatch_2D_host<Patch8> derivs_ymax_8_alloc(derivs_y_idx_range_8);
    DerivMem2_OnPatch_2D_host<Patch9> derivs_ymax_9_alloc(derivs_y_idx_range_9);

    Deriv2_OnPatch_2D_host<Patch1> derivs_ymax_1 = get_field(derivs_ymax_1_alloc);
    Deriv2_OnPatch_2D_host<Patch2> derivs_ymax_2 = get_field(derivs_ymax_2_alloc);
    Deriv2_OnPatch_2D_host<Patch3> derivs_ymax_3 = get_field(derivs_ymax_3_alloc);
    Deriv2_OnPatch_2D_host<Patch4> derivs_ymax_4 = get_field(derivs_ymax_4_alloc);
    Deriv2_OnPatch_2D_host<Patch5> derivs_ymax_5 = get_field(derivs_ymax_5_alloc);
    Deriv2_OnPatch_2D_host<Patch6> derivs_ymax_6 = get_field(derivs_ymax_6_alloc);
    Deriv2_OnPatch_2D_host<Patch7> derivs_ymax_7 = get_field(derivs_ymax_7_alloc);
    Deriv2_OnPatch_2D_host<Patch8> derivs_ymax_8 = get_field(derivs_ymax_8_alloc);
    Deriv2_OnPatch_2D_host<Patch9> derivs_ymax_9 = get_field(derivs_ymax_9_alloc);


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


    MultipatchField<Deriv2_OnPatch_2D_host, Patch1, Patch4, Patch7, Patch2>
            derivs_ymin_147(derivs_ymin_1, derivs_ymin_4, derivs_ymin_7, derivs_ymin_2);
    MultipatchField<Deriv2_OnPatch_2D_host, Patch2, Patch5, Patch8>
            derivs_ymin_258(derivs_ymin_2, derivs_ymin_5, derivs_ymin_8);
    MultipatchField<Deriv2_OnPatch_2D_host, Patch3, Patch6, Patch9>
            derivs_ymin_369(derivs_ymin_3, derivs_ymin_6, derivs_ymin_9);

    MultipatchField<Deriv2_OnPatch_2D_host, Patch1, Patch4, Patch7, Patch2>
            derivs_ymax_147(derivs_ymax_1, derivs_ymax_4, derivs_ymax_7, derivs_ymax_2);
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

    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_1_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_2_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_3_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_4_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_5_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_6_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_7_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_8_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_min_9_alloc(derivs_xy_idx_range);

    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_1(get_field(derivs_xy_min_min_1_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_2(get_field(derivs_xy_min_min_2_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_3(get_field(derivs_xy_min_min_3_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_4(get_field(derivs_xy_min_min_4_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_5(get_field(derivs_xy_min_min_5_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_6(get_field(derivs_xy_min_min_6_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_7(get_field(derivs_xy_min_min_7_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_8(get_field(derivs_xy_min_min_8_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_min_9(get_field(derivs_xy_min_min_9_alloc));


    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_1_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_2_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_3_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_4_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_5_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_6_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_7_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_8_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_min_max_9_alloc(derivs_xy_idx_range);

    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_1(get_field(derivs_xy_min_max_1_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_2(get_field(derivs_xy_min_max_2_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_3(get_field(derivs_xy_min_max_3_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_4(get_field(derivs_xy_min_max_4_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_5(get_field(derivs_xy_min_max_5_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_6(get_field(derivs_xy_min_max_6_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_7(get_field(derivs_xy_min_max_7_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_8(get_field(derivs_xy_min_max_8_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_min_max_9(get_field(derivs_xy_min_max_9_alloc));


    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_1_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_2_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_3_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_4_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_5_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_6_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_7_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_8_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_min_9_alloc(derivs_xy_idx_range);

    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_1(get_field(derivs_xy_max_min_1_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_2(get_field(derivs_xy_max_min_2_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_3(get_field(derivs_xy_max_min_3_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_4(get_field(derivs_xy_max_min_4_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_5(get_field(derivs_xy_max_min_5_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_6(get_field(derivs_xy_max_min_6_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_7(get_field(derivs_xy_max_min_7_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_8(get_field(derivs_xy_max_min_8_alloc));
    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_min_9(get_field(derivs_xy_max_min_9_alloc));


    DerivMem12_OnPatch_2D_host<Patch1> derivs_xy_max_max_1_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch2> derivs_xy_max_max_2_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch3> derivs_xy_max_max_3_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch4> derivs_xy_max_max_4_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch5> derivs_xy_max_max_5_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch6> derivs_xy_max_max_6_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch7> derivs_xy_max_max_7_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch8> derivs_xy_max_max_8_alloc(derivs_xy_idx_range);
    DerivMem12_OnPatch_2D_host<Patch9> derivs_xy_max_max_9_alloc(derivs_xy_idx_range);

    Deriv12_OnPatch_2D_host<Patch1> derivs_xy_max_max_1(get_field(derivs_xy_max_max_1_alloc));
    Deriv12_OnPatch_2D_host<Patch2> derivs_xy_max_max_2(get_field(derivs_xy_max_max_2_alloc));
    Deriv12_OnPatch_2D_host<Patch3> derivs_xy_max_max_3(get_field(derivs_xy_max_max_3_alloc));
    Deriv12_OnPatch_2D_host<Patch4> derivs_xy_max_max_4(get_field(derivs_xy_max_max_4_alloc));
    Deriv12_OnPatch_2D_host<Patch5> derivs_xy_max_max_5(get_field(derivs_xy_max_max_5_alloc));
    Deriv12_OnPatch_2D_host<Patch6> derivs_xy_max_max_6(get_field(derivs_xy_max_max_6_alloc));
    Deriv12_OnPatch_2D_host<Patch7> derivs_xy_max_max_7(get_field(derivs_xy_max_max_7_alloc));
    Deriv12_OnPatch_2D_host<Patch8> derivs_xy_max_max_8(get_field(derivs_xy_max_max_8_alloc));
    Deriv12_OnPatch_2D_host<Patch9> derivs_xy_max_max_9(get_field(derivs_xy_max_max_9_alloc));

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
    initialise_all_functions(functions);
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


    // Initialise the boundary derivatives.
    initialise_all_y_derivatives_min(
            derivs_ymin_789,
            idx_ranges_789,
            evaluator_g,
            get_const_field(function_g_coef));
    initialise_all_y_derivatives_max(
            derivs_ymax_123,
            idx_ranges_123,
            evaluator_g,
            get_const_field(function_g_coef));

    // Initialise the boundary cross-derivatives.
    initialise_all_xy_derivatives_min_min(
            derivs_xy_min_min_789,
            idx_ranges_789,
            evaluator_g,
            get_const_field(function_g_coef),
            std::make_integer_sequence<std::size_t, 3> {});

    initialise_all_xy_derivatives_max_min(
            derivs_xy_max_min_789,
            idx_ranges_789,
            evaluator_g,
            get_const_field(function_g_coef),
            std::make_integer_sequence<std::size_t, 3> {});

    initialise_all_xy_derivatives_min_max(
            derivs_xy_min_max_123,
            idx_ranges_123,
            evaluator_g,
            get_const_field(function_g_coef),
            std::make_integer_sequence<std::size_t, 3> {});

    initialise_all_xy_derivatives_max_max(
            derivs_xy_max_max_123,
            idx_ranges_123,
            evaluator_g,
            get_const_field(function_g_coef),
            std::make_integer_sequence<std::size_t, 3> {});


    // Solve each matrix system ==================================================================
    matrix_123.solve(derivs_xmin_123, derivs_xmax_123, get_const_field(functions_123));
    matrix_456.solve(derivs_xmin_456, derivs_xmax_456, get_const_field(functions_456));
    matrix_789.solve(derivs_xmin_789, derivs_xmax_789, get_const_field(functions_789));


    matrix_258.solve(derivs_ymin_258, derivs_ymax_258, get_const_field(functions_258));
    matrix_369.solve(derivs_ymin_369, derivs_ymax_369, get_const_field(functions_369));
    matrix_147.solve(derivs_ymin_147, derivs_ymax_147, get_const_field(functions_147));

    matrix_456
            .solve(derivs_xy_min_max_456, derivs_xy_max_max_456, get_const_field(derivs_ymax_456));
    matrix_789
            .solve(derivs_xy_min_max_789, derivs_xy_max_max_789, get_const_field(derivs_ymax_789));

    matrix_123
            .solve(derivs_xy_min_min_123, derivs_xy_max_min_123, get_const_field(derivs_ymin_123));
    matrix_456
            .solve(derivs_xy_min_min_456, derivs_xy_max_min_456, get_const_field(derivs_ymin_456));

    // Test the values of the derivatives ========================================================
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

    check_all_x_derivatives(
            derivs_xmin,
            derivs_xmax,
            evaluator_g,
            get_const_field(function_g_coef),
            idx_ranges);

    check_all_y_derivatives(
            derivs_ymin,
            derivs_ymax,
            evaluator_g,
            get_const_field(function_g_coef),
            idx_ranges);

    check_all_xy_derivatives(
            derivs_xy_min_min,
            derivs_xy_max_min,
            derivs_xy_min_max,
            derivs_xy_max_max,
            evaluator_g,
            get_const_field(function_g_coef),
            idx_ranges,
            std::make_integer_sequence<std::size_t, 9> {});


    // Test the spline representations ===========================================================
    check_all_spline_representation_conformity(
            idx_ranges,
            functions,
            derivs_xmin,
            derivs_xmax,
            derivs_ymin,
            derivs_ymax,
            derivs_xy_min_min,
            derivs_xy_max_min,
            derivs_xy_min_max,
            derivs_xy_max_max,
            evaluator_g,
            get_const_field(function_g_coef),
            std::make_integer_sequence<std::size_t, 9> {});
}
