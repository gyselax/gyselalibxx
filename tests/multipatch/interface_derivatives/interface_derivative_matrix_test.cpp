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
#include "interface_derivatives_test_utils.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "view.hpp"


/*
    Test InterfaceDerivativeMatrix on the following geometry:

        |  1  |  2  |  3  |  1 ...
        -------------------
        |  4  |  5  |  6  |  4 ...
        -------------------b
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


template <class Patch>
using DerivFieldMemOnPatch = DerivFieldMem<
        double,
        IdxRange<
                ddc::Deriv<typename Patch::Dim1>,
                typename Patch::Grid1,
                ddc::Deriv<typename Patch::Dim2>,
                typename Patch::Grid2>,
        1>;

template <class Patch>
using DerivFieldMemOnPatch_host = host_t<DerivFieldMemOnPatch<Patch>>;

template <class Patch>
using DerivFieldOnPatch = DerivField<
        double,
        IdxRange<
                ddc::Deriv<typename Patch::Dim1>,
                typename Patch::Grid1,
                ddc::Deriv<typename Patch::Dim2>,
                typename Patch::Grid2>>;

template <class Patch>
using DerivFieldOnPatch_host = host_t<DerivFieldOnPatch<Patch>>;

template <class Patch>
using IdxRange1SliceOnPatch = IdxRangeSlice<typename Patch::Grid1>;

template <class Patch>
using IdxRange2SliceOnPatch = IdxRangeSlice<typename Patch::Grid2>;


// USEFUL FUNCTIONS ==============================================================================
/// @brief Convert a local coordinate into a global coordinate.
Coord<Xg, Yg> get_global_coord(Coord<X, Y> const& local_coord)
{
    double const x = ddc::select<X>(local_coord);
    double const y = ddc::select<Y>(local_coord);
    return Coord<Xg, Yg>(x, y);
}

/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y).
template <class Grid1, class Grid2, class Layout>
void initialise_2D_function(DField<IdxRange<Grid1, Grid2>, Kokkos::HostSpace, Layout> function)
{
    ddc::for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        double const xg = ddc::coordinate(Idx<Grid1>(idx));
        double const yg = ddc::coordinate(Idx<Grid2>(idx));
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI) * Kokkos::sin(yg);
    });
}

/// @brief Initialise all the functions defined on the patches.
// template <class... Patches>
// void initialise_all_functions(MultipatchField<DFieldOnPatch_host, Patches...> const& functions)
// {
//     (initialise_2D_function<typename Patches::Grid1, typename Patches::Grid2>(
//              functions.template get<Patches>()),
//      ...);
// }

template <class... Patches>
void initialise_all_functions(
        MultipatchField<DerivFieldOnPatch_host, Patches...> const& functions_and_derivs)
{
    (initialise_2D_function<typename Patches::Grid1, typename Patches::Grid2>(
             (functions_and_derivs.template get<Patches>()).get_values_field()),
     ...);
}

// /// @brief Initialise the y-derivatives from the global spline.
// template <class Grid1, class Grid2>
// void initialise_2D_y_derivative(
//         host_t<DField<IdxRange<Grid1, ddc::Deriv<typename Grid2::continuous_dimension_type>>>>
//                 deriv_y,
//         Idx<Grid2> const& idx_y,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
// {
//     ddc::for_each(
//             get_idx_range(deriv_y),
//             [&](Idx<Grid1, ddc::Deriv<typename Grid2::continuous_dimension_type>> idx) {
//                 Idx<Grid1, Grid2> idx_xy(Idx<Grid1>(idx), idx_y);
//                 Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx_xy)));
//                 deriv_y(idx) = evaluator_g.deriv_dim_2(interface_coord, function_g_coef);
//             });
// }

// /// @brief Initialise all the y-derivatives on the lower/left bounds.
// template <class... Patches>
// void initialise_all_y_derivatives_min(
//         MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& derivs_ymin,
//         MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
// {
//     (initialise_2D_y_derivative<typename Patches::Grid1, typename Patches::Grid2>(
//              derivs_ymin.template get<Patches>(),
//              typename Patches::IdxRange2(idx_ranges.template get<Patches>()).front(),
//              evaluator_g,
//              function_g_coef),
//      ...);
// }

// /// @brief Initialise all the y-derivatives on the upper/right bounds.
// template <class... Patches>
// void initialise_all_y_derivatives_max(
//         MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& derivs_ymax,
//         MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
// {
//     (initialise_2D_y_derivative<typename Patches::Grid1, typename Patches::Grid2>(
//              derivs_ymax.template get<Patches>(),
//              typename Patches::IdxRange2(idx_ranges.template get<Patches>()).back(),
//              evaluator_g,
//              function_g_coef),
//      ...);
// }


// /// @brief Initialise the cross-derivatives from the global spline.
// template <class Grid1, class Grid2>
// void initialise_2D_xy_derivative(
//         host_t<DField<IdxRange<
//                 ddc::Deriv<typename Grid1::continuous_dimension_type>,
//                 ddc::Deriv<typename Grid2::continuous_dimension_type>>>> deriv_xy,
//         Idx<Grid1> const& idx_x,
//         Idx<Grid2> const& idx_y,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
// {
//     Idx<Grid1, Grid2> idx(idx_x, idx_y);
//     Coord<Xg, Yg> interface_coord(get_global_coord(ddc::coordinate(idx)));
//     deriv_xy(get_idx_range(deriv_xy).front())
//             = evaluator_g.deriv_1_and_2(interface_coord, function_g_coef);
// }

// /// @brief Initialise all the cross-derivatives on the lower/left - lower/left corners.
// template <std::size_t... I, class... Patches>
// void initialise_all_xy_derivatives_min_min(
//         std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_min_min,
//         MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
//         std::integer_sequence<std::size_t, I...>)
// {
//     (initialise_2D_xy_derivative<typename Patches::Grid1, typename Patches::Grid2>(
//              std::get<I>(derivs_min_min),
//              typename Patches::IdxRange1(idx_ranges.template get<Patches>()).front(),
//              typename Patches::IdxRange2(idx_ranges.template get<Patches>()).front(),
//              evaluator_g,
//              function_g_coef),
//      ...);
// }

// /// @brief Initialise all the cross-derivatives on the upper/right - lower/left corners.
// template <std::size_t... I, class... Patches>
// void initialise_all_xy_derivatives_max_min(
//         std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_max_min,
//         MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
//         std::integer_sequence<std::size_t, I...>)
// {
//     (initialise_2D_xy_derivative<typename Patches::Grid1, typename Patches::Grid2>(
//              std::get<I>(derivs_max_min),
//              typename Patches::IdxRange1(idx_ranges.template get<Patches>()).back(),
//              typename Patches::IdxRange2(idx_ranges.template get<Patches>()).front(),
//              evaluator_g,
//              function_g_coef),
//      ...);
// }

// /// @brief Initialise all the cross-derivatives on the lower/left - upper/right corners.
// template <std::size_t... I, class... Patches>
// void initialise_all_xy_derivatives_min_max(
//         std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_min_max,
//         MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
//         std::integer_sequence<std::size_t, I...>)
// {
//     (initialise_2D_xy_derivative<typename Patches::Grid1, typename Patches::Grid2>(
//              std::get<I>(derivs_min_max),
//              typename Patches::IdxRange1(idx_ranges.template get<Patches>()).front(),
//              typename Patches::IdxRange2(idx_ranges.template get<Patches>()).back(),
//              evaluator_g,
//              function_g_coef),
//      ...);
// }

// /// @brief Initialise all the cross-derivatives on the upper/right - lower/left corners.
// template <std::size_t... I, class... Patches>
// void initialise_all_xy_derivatives_max_max(
//         std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_max_max,
//         MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
//         std::integer_sequence<std::size_t, I...>)
// {
//     (initialise_2D_xy_derivative<typename Patches::Grid1, typename Patches::Grid2>(
//              std::get<I>(derivs_max_max),
//              typename Patches::IdxRange1(idx_ranges.template get<Patches>()).back(),
//              typename Patches::IdxRange2(idx_ranges.template get<Patches>()).back(),
//              evaluator_g,
//              function_g_coef),
//      ...);
// }


template <class Patch>
void copy_function_values_in_deriv_field_on_patch(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        DFieldOnPatch_host<Patch> const& function)
{
    DField<typename Patch::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_extracted
            = function_and_derivs.get_values_field();
    ddc::for_each(get_idx_range(function), [&](typename Patch::Idx12 const idx) {
        function_extracted(idx) = function(idx);
    });
}

template <class... Patches>
void copy_function_values_in_deriv_field(
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        MultipatchField<DFieldOnPatch_host, Patches...> const& functions)
{
    (copy_function_values_in_deriv_field_on_patch<Patches>(
             functions_and_derivs.template get<Patches>(),
             functions.template get<Patches>()),
     ...);
}


template <class Patch>
void copy_deriv_x_in_deriv_field_on_patch(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        Deriv1_OnPatch_2D_host<Patch> const& deriv_xmin,
        Deriv1_OnPatch_2D_host<Patch> const& deriv_xmax,
        IdxRangeSlice<typename Patch::Grid1> const& idx_range_slice_dx)
{
    using dX = typename ddc::Deriv<typename Patch::Grid1>;
    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    Idx<dX, typename Patch::Grid1> idx_slice_deriv_xmin(Idx<dX>(1), idx_range_slice_dx.front());
    Idx<dX, typename Patch::Grid1> idx_slice_deriv_xmax(Idx<dX>(1), idx_range_slice_dx.back());

    DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_xmin_extracted = function_and_derivs[idx_slice_deriv_xmin];
    DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_xmax_extracted = function_and_derivs[idx_slice_deriv_xmax];

    ddc::for_each(get_idx_range(derivs_xmin_extracted), [&](typename Patch::Idx2 const idx_y) {
        derivs_xmin_extracted(idx_y) = deriv_xmin(Idx<DerivX>(1), idx_y);
        derivs_xmax_extracted(idx_y) = deriv_xmax(Idx<DerivX>(1), idx_y);
    });
}

template <class... Patches>
void copy_derivs_x_in_deriv_field(
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        MultipatchField<Deriv1_OnPatch_2D_host, Patches...> const& derivs_xmin,
        MultipatchField<Deriv1_OnPatch_2D_host, Patches...> const& derivs_xmax,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_ranges_slice_dx)
{
    (copy_deriv_x_in_deriv_field_on_patch<Patches>(
             functions_and_derivs.template get<Patches>(),
             derivs_xmin.template get<Patches>(),
             derivs_xmax.template get<Patches>(),
             idx_ranges_slice_dx.template get<Patches>()),
     ...);
}



template <class Patch>
void copy_deriv_y_in_deriv_field_on_patch(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        Deriv2_OnPatch_2D_host<Patch> const& deriv_ymin,
        Deriv2_OnPatch_2D_host<Patch> const& deriv_ymax,
        IdxRangeSlice<typename Patch::Grid2> const& idx_range_slice_dy)
{
    using dY = typename ddc::Deriv<typename Patch::Grid2>;
    using DerivY = typename ddc::Deriv<typename Patch::Dim2>;
    Idx<dY, typename Patch::Grid2> idx_slice_deriv_ymin(Idx<dY>(1), idx_range_slice_dy.front());
    Idx<dY, typename Patch::Grid2> idx_slice_deriv_ymax(Idx<dY>(1), idx_range_slice_dy.back());

    DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymin_extracted = function_and_derivs[idx_slice_deriv_ymin];
    DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymax_extracted = function_and_derivs[idx_slice_deriv_ymax];

    ddc::for_each(get_idx_range(derivs_ymin_extracted), [&](typename Patch::Idx1 const idx_x) {
        derivs_ymin_extracted(idx_x) = deriv_ymin(Idx<DerivY>(1), idx_x);
        derivs_ymax_extracted(idx_x) = deriv_ymax(Idx<DerivY>(1), idx_x);
    });
}

template <class... Patches>
void copy_derivs_y_in_deriv_field(
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& derivs_ymin,
        MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& derivs_ymax,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_ranges_slice_dy)
{
    (copy_deriv_y_in_deriv_field_on_patch<Patches>(
             functions_and_derivs.template get<Patches>(),
             derivs_ymin.template get<Patches>(),
             derivs_ymax.template get<Patches>(),
             idx_ranges_slice_dy.template get<Patches>()),
     ...);
}


template <class Patch>
void copy_deriv_xy_in_deriv_field_on_patch(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        Deriv12_OnPatch_2D_host<Patch> const& derivs_xy_min_min,
        Deriv12_OnPatch_2D_host<Patch> const& derivs_xy_max_min,
        Deriv12_OnPatch_2D_host<Patch> const& derivs_xy_min_max,
        Deriv12_OnPatch_2D_host<Patch> const& derivs_xy_max_max,
        IdxRangeSlice<typename Patch::Grid1> const& idx_range_slice_dx,
        IdxRangeSlice<typename Patch::Grid2> const& idx_range_slice_dy)
{
    using dX = typename ddc::Deriv<typename Patch::Grid1>;
    using dY = typename ddc::Deriv<typename Patch::Grid2>;
    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    using DerivY = typename ddc::Deriv<typename Patch::Dim2>;

    IdxRange<dX> deriv_block_x(Idx<dX>(1), IdxStep<dX>(2));
    IdxRange<dY> deriv_block_y(Idx<dY>(1), IdxStep<dY>(2));
    IdxRange<dX, dY> deriv_block_xy(deriv_block_x, deriv_block_y);

    detail::ViewNDMaker<4, double, false>::type derivs_xy_extracted
            = function_and_derivs.get_mdspan(deriv_block_xy);

    Idx<DerivX, DerivY> idx(1, 1);

    derivs_xy_extracted(IdxStep<dX>(0), IdxStep<dY>(0), 0, 0) = derivs_xy_min_min(idx);
    derivs_xy_extracted(IdxStep<dX>(1), IdxStep<dY>(0), 1, 0) = derivs_xy_max_min(idx);
    derivs_xy_extracted(IdxStep<dX>(0), IdxStep<dY>(1), 0, 1) = derivs_xy_min_max(idx);
    derivs_xy_extracted(IdxStep<dX>(1), IdxStep<dY>(1), 1, 1) = derivs_xy_max_max(idx);
}

template <class... Patches>
void copy_derivs_xy_in_deriv_field(
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_xy_min_min,
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_xy_max_min,
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_xy_min_max,
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_xy_max_max,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_ranges_slice_dx,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_ranges_slice_dy)
{
    (copy_deriv_xy_in_deriv_field_on_patch<Patches>(
             functions_and_derivs.template get<Patches>(),
             std::get<0>(derivs_xy_min_min),
             std::get<0>(derivs_xy_max_min),
             std::get<0>(derivs_xy_min_max),
             std::get<0>(derivs_xy_max_max),
             idx_ranges_slice_dx.template get<Patches>(),
             idx_ranges_slice_dy.template get<Patches>()),
     ...);
}



template <class Patch>
void copy_deriv_field_in_function_values_on_patch(
        DFieldOnPatch_host<Patch> const& function,
        DerivFieldOnPatch_host<Patch> function_and_derivs)
{
    DField<typename Patch::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_extracted
            = function_and_derivs.get_values_field();
    ddc::for_each(get_idx_range(function), [&](typename Patch::Idx12 const idx) {
        function(idx) = function_extracted(idx);
    });
}

template <class... Patches>
void copy_deriv_field_in_function_values(
        MultipatchField<DFieldOnPatch_host, Patches...> const& functions,
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
{
    (copy_deriv_field_in_function_values_on_patch<Patches>(
             functions.template get<Patches>(),
             functions_and_derivs.template get<Patches>()),
     ...);
}


template <class Patch>
void copy_deriv_field_in_derivs_x_on_patch(
        Deriv1_OnPatch_2D_host<Patch> const& deriv_xmin,
        Deriv1_OnPatch_2D_host<Patch> const& deriv_xmax,
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRangeSlice<typename Patch::Grid1> const& idx_range_slice_dx)
{
    using dX = typename ddc::Deriv<typename Patch::Grid1>;
    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    Idx<dX, typename Patch::Grid1> idx_slice_deriv_xmin(Idx<dX>(1), idx_range_slice_dx.front());
    Idx<dX, typename Patch::Grid1> idx_slice_deriv_xmax(Idx<dX>(1), idx_range_slice_dx.back());

    DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_xmin_extracted = function_and_derivs[idx_slice_deriv_xmin];
    DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_xmax_extracted = function_and_derivs[idx_slice_deriv_xmax];

    ddc::for_each(get_idx_range(derivs_xmin_extracted), [&](typename Patch::Idx2 const idx_y) {
        deriv_xmin(Idx<DerivX>(1), idx_y) = derivs_xmin_extracted(idx_y);
        deriv_xmax(Idx<DerivX>(1), idx_y) = derivs_xmax_extracted(idx_y);
    });
}

template <class... Patches>
void copy_deriv_field_in_derivs_x(
        MultipatchField<Deriv1_OnPatch_2D_host, Patches...> const& derivs_xmin,
        MultipatchField<Deriv1_OnPatch_2D_host, Patches...> const& derivs_xmax,
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_ranges_slice_dx)
{
    (copy_deriv_field_in_derivs_x_on_patch<Patches>(
             derivs_xmin.template get<Patches>(),
             derivs_xmax.template get<Patches>(),
             functions_and_derivs.template get<Patches>(),
             idx_ranges_slice_dx.template get<Patches>()),
     ...);
}



template <class Patch>
void copy_deriv_field_in_derivs_y_on_patch(
        Deriv2_OnPatch_2D_host<Patch> const& deriv_ymin,
        Deriv2_OnPatch_2D_host<Patch> const& deriv_ymax,
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRangeSlice<typename Patch::Grid2> const& idx_range_slice_dy)
{
    using dY = typename ddc::Deriv<typename Patch::Grid2>;
    using DerivY = typename ddc::Deriv<typename Patch::Dim2>;
    Idx<dY, typename Patch::Grid2> idx_slice_deriv_ymin(Idx<dY>(1), idx_range_slice_dy.front());
    Idx<dY, typename Patch::Grid2> idx_slice_deriv_ymax(Idx<dY>(1), idx_range_slice_dy.back());

    DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymin_extracted = function_and_derivs[idx_slice_deriv_ymin];
    DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymax_extracted = function_and_derivs[idx_slice_deriv_ymax];

    ddc::for_each(get_idx_range(derivs_ymin_extracted), [&](typename Patch::Idx1 const idx_x) {
        deriv_ymin(Idx<DerivY>(1), idx_x) = derivs_ymin_extracted(idx_x);
        deriv_ymax(Idx<DerivY>(1), idx_x) = derivs_ymax_extracted(idx_x);
    });
}

template <class... Patches>
void copy_deriv_field_in_derivs_y(
        MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& derivs_ymin,
        MultipatchField<Deriv2_OnPatch_2D_host, Patches...> const& derivs_ymax,
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_ranges_slice_dy)
{
    (copy_deriv_field_in_derivs_y_on_patch<Patches>(
             derivs_ymin.template get<Patches>(),
             derivs_ymax.template get<Patches>(),
             functions_and_derivs.template get<Patches>(),
             idx_ranges_slice_dy.template get<Patches>()),
     ...);
}


template <class Patch>
void copy_deriv_field_in_deriv_xy_on_patch(
        Deriv12_OnPatch_2D_host<Patch> const& derivs_xy_min_min,
        Deriv12_OnPatch_2D_host<Patch> const& derivs_xy_max_min,
        Deriv12_OnPatch_2D_host<Patch> const& derivs_xy_min_max,
        Deriv12_OnPatch_2D_host<Patch> const& derivs_xy_max_max,
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRangeSlice<typename Patch::Grid1> const& idx_range_slice_dx,
        IdxRangeSlice<typename Patch::Grid2> const& idx_range_slice_dy)
{
    using dX = typename ddc::Deriv<typename Patch::Grid1>;
    using dY = typename ddc::Deriv<typename Patch::Grid2>;
    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    using DerivY = typename ddc::Deriv<typename Patch::Dim2>;

    IdxRange<dX> deriv_block_x(Idx<dX>(1), IdxStep<dX>(2));
    IdxRange<dY> deriv_block_y(Idx<dY>(1), IdxStep<dY>(2));
    IdxRange<dX, dY> deriv_block_xy(deriv_block_x, deriv_block_y);

    detail::ViewNDMaker<4, double, false>::type derivs_xy_extracted
            = function_and_derivs.get_mdspan(deriv_block_xy);

    Idx<DerivX, DerivY> idx(1, 1);

    derivs_xy_min_min(idx) = derivs_xy_extracted(IdxStep<dX>(0), IdxStep<dY>(0), 0, 0);
    derivs_xy_max_min(idx) = derivs_xy_extracted(IdxStep<dX>(1), IdxStep<dY>(0), 1, 0);
    derivs_xy_min_max(idx) = derivs_xy_extracted(IdxStep<dX>(0), IdxStep<dY>(1), 0, 1);
    derivs_xy_max_max(idx) = derivs_xy_extracted(IdxStep<dX>(1), IdxStep<dY>(1), 1, 1);
}

template <class... Patches>
void copy_deriv_field_in_derivs_xy(
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_xy_min_min,
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_xy_max_min,
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_xy_min_max,
        std::tuple<Deriv12_OnPatch_2D_host<Patches>...> const& derivs_xy_max_max,
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_ranges_slice_dx,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_ranges_slice_dy)
{
    (copy_deriv_field_in_deriv_xy_on_patch<Patches>(
             std::get<0>(derivs_xy_min_min),
             std::get<0>(derivs_xy_max_min),
             std::get<0>(derivs_xy_min_max),
             std::get<0>(derivs_xy_max_max),
             functions_and_derivs.template get<Patches>(),
             idx_ranges_slice_dx.template get<Patches>(),
             idx_ranges_slice_dy.template get<Patches>()),
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
            MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_dx)
    {
        (check_x_derivatives<Patches>(
                 functions_and_derivs.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 typename Patches::IdxRange1(idx_ranges.template get<Patches>()),
                 idx_range_slices_dx.template get<Patches>()),
         ...);
    }

    /** @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
     * the interfaces for a given patch.
     */
    template <class Patch>
    void check_x_derivatives(
            DerivFieldOnPatch_host<Patch> function_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::IdxRange1 const& idx_range_perp,
            IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx)
    {
        using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
        Idx<DerivX> idx_deriv(1);

        Idx<DerivX, typename Patch::Grid1> idx_slice_min(idx_deriv, idx_range_slice_dx.front());
        Idx<DerivX, typename Patch::Grid1> idx_slice_max(idx_deriv, idx_range_slice_dx.back());

        DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_xmin_extracted = function_and_derivs[idx_slice_min];
        DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_xmax_extracted = function_and_derivs[idx_slice_max];

        typename Patch::IdxRange2 idx_range_par(get_idx_range(derivs_xmin_extracted));
        ddc::for_each(idx_range_par, [&](typename Patch::Idx2 const& idx_par) {
            typename Patch::Idx12 idx_min(idx_range_perp.front(), idx_par);
            typename Patch::Idx12 idx_max(idx_range_perp.back(), idx_par);
            Coord<Xg, Yg> interface_coord_min(get_global_coord(ddc::coordinate(idx_min)));
            Coord<Xg, Yg> interface_coord_max(get_global_coord(ddc::coordinate(idx_max)));

            double const global_deriv_min
                    = evaluator_g.deriv_dim_1(interface_coord_min, function_g_coef);
            double const global_deriv_max
                    = evaluator_g.deriv_dim_1(interface_coord_max, function_g_coef);

            EXPECT_NEAR(derivs_xmin_extracted(idx_par), global_deriv_min, 5e-14);
            EXPECT_NEAR(derivs_xmax_extracted(idx_par), global_deriv_max, 5e-14);
        });
    }


    /** @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
     * the interfaces.
     */
    template <class... Patches>
    void check_all_y_derivatives(
            MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_dy)
    {
        (check_y_derivatives<Patches>(
                 functions_and_derivs.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 typename Patches::IdxRange2(idx_ranges.template get<Patches>()),
                 idx_range_slices_dy.template get<Patches>()),
         ...);
    }

    /** @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
     * the interfaces for a given patch.
     */
    template <class Patch>
    void check_y_derivatives(
            DerivFieldOnPatch_host<Patch> function_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::IdxRange2 const& idx_range_perp,
            IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy)
    {
        using DerivY = typename ddc::Deriv<typename Patch::Dim2>;
        Idx<DerivY> idx_deriv(1);

        Idx<DerivY, typename Patch::Grid2> idx_slice_min(idx_deriv, idx_range_slice_dy.front());
        Idx<DerivY, typename Patch::Grid2> idx_slice_max(idx_deriv, idx_range_slice_dy.back());

        DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_ymin_extracted = function_and_derivs[idx_slice_min];
        DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_ymax_extracted = function_and_derivs[idx_slice_max];


        typename Patch::IdxRange1 idx_range_par(get_idx_range(derivs_ymin_extracted));
        ddc::for_each(idx_range_par, [&](typename Patch::Idx1 const& idx_par) {
            typename Patch::Idx12 idx_min(idx_par, idx_range_perp.front());
            typename Patch::Idx12 idx_max(idx_par, idx_range_perp.back());
            Coord<Xg, Yg> interface_coord_min(get_global_coord(ddc::coordinate(idx_min)));
            Coord<Xg, Yg> interface_coord_max(get_global_coord(ddc::coordinate(idx_max)));

            double const global_deriv_min
                    = evaluator_g.deriv_dim_2(interface_coord_min, function_g_coef);
            double const global_deriv_max
                    = evaluator_g.deriv_dim_2(interface_coord_max, function_g_coef);

            if constexpr (!ddc::in_tags_v<Patch, ddc::detail::TypeSeq<Patch1, Patch2, Patch3>>) {
                EXPECT_NEAR(derivs_ymax_extracted(idx_par), global_deriv_max, 5e-14);
            }
            if constexpr (!ddc::in_tags_v<Patch, ddc::detail::TypeSeq<Patch7, Patch8, Patch9>>) {
                EXPECT_NEAR(derivs_ymin_extracted(idx_par), global_deriv_min, 5e-14);
            }
        });
    }


    /** @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
     * the interfaces.
     */
    template <class... Patches>
    void check_all_xy_derivatives(
            MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_dx,
            MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_dy)


    {
        (check_xy_derivatives<Patches>(
                 functions_and_derivs.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 idx_ranges.template get<Patches>(),
                 idx_range_slices_dx.template get<Patches>(),
                 idx_range_slices_dy.template get<Patches>()),
         ...);
    }

    /** @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
     * the interfaces for a given patch.
     */
    template <class Patch>
    void check_xy_derivatives(
            DerivFieldOnPatch_host<Patch> function_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::IdxRange12 const& idx_range,
            IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx,
            IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy)
    {
        using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
        using DerivY = typename ddc::Deriv<typename Patch::Dim2>;
        Idx<DerivX> idx_deriv_x(1);
        Idx<DerivY> idx_deriv_y(1);
        IdxStep<DerivX> nelems_deriv_x(1);
        IdxStep<DerivY> nelems_deriv_y(1);
        IdxRange<DerivX> idx_range_deriv_x(idx_deriv_x, nelems_deriv_x);
        IdxRange<DerivY> idx_range_deriv_y(idx_deriv_y, nelems_deriv_y);
        IdxRange<DerivX, DerivY> idx_range_deriv_xy(idx_range_deriv_x, idx_range_deriv_y);

        DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_xy_min_min_extracted
                = function_and_derivs[idx_range_deriv_xy][idx_range_slice_dx.front()]
                                     [idx_range_slice_dy.front()];
        DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_xy_max_min_extracted
                = function_and_derivs[idx_range_deriv_xy][idx_range_slice_dx.back()]
                                     [idx_range_slice_dy.front()];
        DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_xy_min_max_extracted
                = function_and_derivs[idx_range_deriv_xy][idx_range_slice_dx.front()]
                                     [idx_range_slice_dy.back()];
        DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_xy_max_max_extracted
                = function_and_derivs[idx_range_deriv_xy][idx_range_slice_dx.back()]
                                     [idx_range_slice_dy.back()];


        typename Patch::Idx1 idx_1min(idx_range.front());
        typename Patch::Idx1 idx_1max(idx_range.back());
        typename Patch::Idx2 idx_2min(idx_range.front());
        typename Patch::Idx2 idx_2max(idx_range.back());

        typename Patch::Idx12 idx_min_min(idx_1min, idx_2min);
        typename Patch::Idx12 idx_max_min(idx_1max, idx_2min);
        typename Patch::Idx12 idx_min_max(idx_1min, idx_2max);
        typename Patch::Idx12 idx_max_max(idx_1max, idx_2max);

        Coord<Xg, Yg> interface_coord_min_min(get_global_coord(ddc::coordinate(idx_min_min)));
        Coord<Xg, Yg> interface_coord_max_min(get_global_coord(ddc::coordinate(idx_max_min)));
        Coord<Xg, Yg> interface_coord_min_max(get_global_coord(ddc::coordinate(idx_min_max)));
        Coord<Xg, Yg> interface_coord_max_max(get_global_coord(ddc::coordinate(idx_max_max)));

        double const global_deriv_min_min
                = evaluator_g.deriv_1_and_2(interface_coord_min_min, function_g_coef);
        double const global_deriv_max_min
                = evaluator_g.deriv_1_and_2(interface_coord_max_min, function_g_coef);
        double const global_deriv_min_max
                = evaluator_g.deriv_1_and_2(interface_coord_min_max, function_g_coef);
        double const global_deriv_max_max
                = evaluator_g.deriv_1_and_2(interface_coord_max_max, function_g_coef);

        EXPECT_NEAR(
                derivs_xy_min_min_extracted(idx_deriv_x, idx_deriv_y),
                global_deriv_min_min,
                6e-14);
        EXPECT_NEAR(
                derivs_xy_max_min_extracted(idx_deriv_x, idx_deriv_y),
                global_deriv_max_min,
                6e-14);
        EXPECT_NEAR(
                derivs_xy_min_max_extracted(idx_deriv_x, idx_deriv_y),
                global_deriv_min_max,
                6e-14);
        EXPECT_NEAR(
                derivs_xy_max_max_extracted(idx_deriv_x, idx_deriv_y),
                global_deriv_max_max,
                6e-14);
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
            decltype(deriv_calculators_collect_123),
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            Kokkos::DefaultHostExecutionSpace,
            Patch1,
            Patch2,
            Patch3>
            matrix_123(idx_ranges_123, deriv_calculators_collect_123);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<4>,
            decltype(deriv_calculators_collect_456),
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            Kokkos::DefaultHostExecutionSpace,
            Patch4,
            Patch5,
            Patch6>
            matrix_456(idx_ranges_456, deriv_calculators_collect_456);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<7>,
            decltype(deriv_calculators_collect_789),
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            Kokkos::DefaultHostExecutionSpace,
            Patch7,
            Patch8,
            Patch9>
            matrix_789(idx_ranges_789, deriv_calculators_collect_789);


    // Test with an exact patch to check it will only take the needed patches.
    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<1>,
            decltype(deriv_calculators_collect_147),
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            Kokkos::DefaultHostExecutionSpace,
            Patch1,
            Patch4,
            Patch7,
            Patch2>
            matrix_147(idx_ranges_147, deriv_calculators_collect_147);

    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<2>,
            decltype(deriv_calculators_collect_258),
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            Kokkos::DefaultHostExecutionSpace,
            Patch2,
            Patch5,
            Patch8>
            matrix_258(idx_ranges_258, deriv_calculators_collect_258);


    InterfaceDerivativeMatrix<
            Connectivity,
            GridY<3>,
            decltype(deriv_calculators_collect_369),
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            Kokkos::DefaultHostExecutionSpace,
            Patch3,
            Patch6,
            Patch9>
            matrix_369(idx_ranges_369, deriv_calculators_collect_369);

    // Instantiate test function values ==========================================================

    host_t<DFieldMem<IdxRange<GridXg, GridYg>>> function_g_alloc(idx_range_xy_g);
    host_t<DField<IdxRange<GridXg, GridYg>>> function_g = get_field(function_g_alloc);


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


    // Try to use DerivField =====================================================================
    // Instantiate index range slices ------------------------------------------------------------
    IdxRangeSlice<GridX<1>>
            idx_range_slice_dx1(idx_range_x1.front(), IdxStep<GridX<1>>(2), idx_range_x1.extents());
    IdxRangeSlice<GridX<2>>
            idx_range_slice_dx2(idx_range_x2.front(), IdxStep<GridX<2>>(2), idx_range_x2.extents());
    IdxRangeSlice<GridX<3>>
            idx_range_slice_dx3(idx_range_x3.front(), IdxStep<GridX<3>>(2), idx_range_x3.extents());
    IdxRangeSlice<GridX<4>>
            idx_range_slice_dx4(idx_range_x4.front(), IdxStep<GridX<4>>(2), idx_range_x4.extents());
    IdxRangeSlice<GridX<5>>
            idx_range_slice_dx5(idx_range_x5.front(), IdxStep<GridX<5>>(2), idx_range_x5.extents());
    IdxRangeSlice<GridX<6>>
            idx_range_slice_dx6(idx_range_x6.front(), IdxStep<GridX<6>>(2), idx_range_x6.extents());
    IdxRangeSlice<GridX<7>>
            idx_range_slice_dx7(idx_range_x7.front(), IdxStep<GridX<7>>(2), idx_range_x7.extents());
    IdxRangeSlice<GridX<8>>
            idx_range_slice_dx8(idx_range_x8.front(), IdxStep<GridX<8>>(2), idx_range_x8.extents());
    IdxRangeSlice<GridX<9>>
            idx_range_slice_dx9(idx_range_x9.front(), IdxStep<GridX<9>>(2), idx_range_x9.extents());

    IdxRangeSlice<GridY<1>>
            idx_range_slice_dy1(idx_range_y1.front(), IdxStep<GridY<1>>(2), idx_range_y1.extents());
    IdxRangeSlice<GridY<2>>
            idx_range_slice_dy2(idx_range_y2.front(), IdxStep<GridY<2>>(2), idx_range_y2.extents());
    IdxRangeSlice<GridY<3>>
            idx_range_slice_dy3(idx_range_y3.front(), IdxStep<GridY<3>>(2), idx_range_y3.extents());
    IdxRangeSlice<GridY<4>>
            idx_range_slice_dy4(idx_range_y4.front(), IdxStep<GridY<4>>(2), idx_range_y4.extents());
    IdxRangeSlice<GridY<5>>
            idx_range_slice_dy5(idx_range_y5.front(), IdxStep<GridY<5>>(2), idx_range_y5.extents());
    IdxRangeSlice<GridY<6>>
            idx_range_slice_dy6(idx_range_y6.front(), IdxStep<GridY<6>>(2), idx_range_y6.extents());
    IdxRangeSlice<GridY<7>>
            idx_range_slice_dy7(idx_range_y7.front(), IdxStep<GridY<7>>(2), idx_range_y7.extents());
    IdxRangeSlice<GridY<8>>
            idx_range_slice_dy8(idx_range_y8.front(), IdxStep<GridY<8>>(2), idx_range_y8.extents());
    IdxRangeSlice<GridY<9>>
            idx_range_slice_dy9(idx_range_y9.front(), IdxStep<GridY<9>>(2), idx_range_y9.extents());

    MultipatchType<
            IdxRange1SliceOnPatch,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            idx_ranges_slice_dx(
                    idx_range_slice_dx1,
                    idx_range_slice_dx2,
                    idx_range_slice_dx3,
                    idx_range_slice_dx4,
                    idx_range_slice_dx5,
                    idx_range_slice_dx6,
                    idx_range_slice_dx7,
                    idx_range_slice_dx8,
                    idx_range_slice_dx9);

    MultipatchType<
            IdxRange2SliceOnPatch,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            idx_ranges_slice_dy(
                    idx_range_slice_dy1,
                    idx_range_slice_dy2,
                    idx_range_slice_dy3,
                    idx_range_slice_dy4,
                    idx_range_slice_dy5,
                    idx_range_slice_dy6,
                    idx_range_slice_dy7,
                    idx_range_slice_dy8,
                    idx_range_slice_dy9);

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


    // Initialise the data =======================================================================
    initialise_all_functions(functions_and_derivs);
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


    // // Initialise the boundary derivatives.
    // initialise_all_y_derivatives_min(
    //         derivs_ymin_789,
    //         idx_ranges_789,
    //         evaluator_g,
    //         get_const_field(function_g_coef));
    // initialise_all_y_derivatives_max(
    //         derivs_ymax_123,
    //         idx_ranges_123,
    //         evaluator_g,
    //         get_const_field(function_g_coef));

    // // Initialise the boundary cross-derivatives.
    // initialise_all_xy_derivatives_min_min(
    //         derivs_xy_min_min_789,
    //         idx_ranges_789,
    //         evaluator_g,
    //         get_const_field(function_g_coef),
    //         std::make_integer_sequence<std::size_t, 3> {});

    // initialise_all_xy_derivatives_max_min(
    //         derivs_xy_max_min_789,
    //         idx_ranges_789,
    //         evaluator_g,
    //         get_const_field(function_g_coef),
    //         std::make_integer_sequence<std::size_t, 3> {});

    // initialise_all_xy_derivatives_min_max(
    //         derivs_xy_min_max_123,
    //         idx_ranges_123,
    //         evaluator_g,
    //         get_const_field(function_g_coef),
    //         std::make_integer_sequence<std::size_t, 3> {});

    // initialise_all_xy_derivatives_max_max(
    //         derivs_xy_max_max_123,
    //         idx_ranges_123,
    //         evaluator_g,
    //         get_const_field(function_g_coef),
    //         std::make_integer_sequence<std::size_t, 3> {});

    // Use the InterfaceDerivativeMatrix ---------------------------------------------------------
    // MultipatchField<
    //         DerivFieldOnPatch_host,
    //         Patch1,
    //         Patch2,
    //         Patch3,
    //         Patch4,
    //         Patch5,
    //         Patch6,
    //         Patch7,
    //         Patch8,
    //         Patch9>
    //         functions_and_derivs(
    //                 function_and_derivs_1,
    //                 function_and_derivs_2,
    //                 function_and_derivs_3,
    //                 function_and_derivs_4,
    //                 function_and_derivs_5,
    //                 function_and_derivs_6,
    //                 function_and_derivs_7,
    //                 function_and_derivs_8,
    //                 function_and_derivs_9);

    MultipatchField<DerivFieldOnPatch_host, Patch1, Patch2, Patch3> functions_and_derivs_123(
            function_and_derivs_1,
            function_and_derivs_2,
            function_and_derivs_3);

    MultipatchField<DerivFieldOnPatch_host, Patch4, Patch5, Patch6> functions_and_derivs_456(
            function_and_derivs_4,
            function_and_derivs_5,
            function_and_derivs_6);

    MultipatchField<DerivFieldOnPatch_host, Patch7, Patch8, Patch9> functions_and_derivs_789(
            function_and_derivs_7,
            function_and_derivs_8,
            function_and_derivs_9);

    MultipatchField<DerivFieldOnPatch_host, Patch1, Patch4, Patch7, Patch2>
            functions_and_derivs_147(
                    function_and_derivs_1,
                    function_and_derivs_4,
                    function_and_derivs_7,
                    function_and_derivs_2);

    MultipatchField<DerivFieldOnPatch_host, Patch2, Patch5, Patch8> functions_and_derivs_258(
            function_and_derivs_2,
            function_and_derivs_5,
            function_and_derivs_8);

    MultipatchField<DerivFieldOnPatch_host, Patch3, Patch6, Patch9> functions_and_derivs_369(
            function_and_derivs_3,
            function_and_derivs_6,
            function_and_derivs_9);

    std::cout << "before solving." << std::endl;
    matrix_123.solve_deriv(functions_and_derivs_123);
    matrix_456.solve_deriv(functions_and_derivs_456);
    matrix_789.solve_deriv(functions_and_derivs_789);

    matrix_147.solve_deriv(functions_and_derivs_147);
    matrix_258.solve_deriv(functions_and_derivs_258);
    matrix_369.solve_deriv(functions_and_derivs_369);

    matrix_456.solve_cross_deriv(functions_and_derivs_456);
    matrix_789.solve_cross_deriv(functions_and_derivs_789);

    matrix_123.solve_cross_deriv(functions_and_derivs_123);
    matrix_456.solve_cross_deriv(functions_and_derivs_456);
    std::cout << "after solving." << std::endl;

    // Test the values of the derivatives ========================================================
    check_all_x_derivatives(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            idx_ranges,
            idx_ranges_slice_dx);

    check_all_y_derivatives(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            idx_ranges,
            idx_ranges_slice_dy);

    // check_all_xy_derivatives(
    //         functions_and_derivs,
    //         evaluator_g,
    //         get_const_field(function_g_coef),
    //         idx_ranges,
    //         idx_ranges_slice_dx,
    //         idx_ranges_slice_dy);


    // check_all_spline_representation_conformity(
    //         idx_ranges,
    //         functions,
    //         derivs_xmin,
    //         derivs_xmax,
    //         derivs_ymin,
    //         derivs_ymax,
    //         derivs_xy_min_min,
    //         derivs_xy_max_min,
    //         derivs_xy_min_max,
    //         derivs_xy_max_max,
    //         evaluator_g,
    //         get_const_field(function_g_coef),
    //         std::make_integer_sequence<std::size_t, 9> {});
}
