// SPDX-License-Identifier: MIT

#pragma once

#include <typeinfo>

#include <ddc/ddc.hpp>

// #include "connectivity_details.hpp"
// #include "ddc_aliases.hpp"
// #include "edge.hpp"
// #include "edge_transformation.hpp"
// #include "interface.hpp"
// #include "interface_derivative_coefficients_collection.hpp"
#include "multipatch_field.hpp"
// #include "multipatch_type.hpp"
// #include "types.hpp"
// #include "view.hpp"
#include "multipatch_spline_builder_2d.hpp"
#include "multipatch_spline_evaluator_2d.hpp"


template <
        class ExecSpace,
        class MemorySpace,
        template <typename P>
        typename BSpline1OnPatch,
        template <typename P>
        typename BSpline2OnPatch,
        template <typename P>
        typename Grid1OnPatch,
        template <typename P>
        typename Grid2OnPatch,
        ddc::SplineBuilderClosure SBCLower1,
        ddc::SplineBuilderClosure SBCUpper1,
        ddc::SplineBuilderClosure SBCLower2,
        ddc::SplineBuilderClosure SBCUpper2,
        ddc::SplineSolver Solver,
        class ExtrapolationRule,
        class PatchLocator,
        class Connectivity,
        // class InterfacesDerivCaculatorCollection,
        template <typename P>
        typename ValuesOnPatch,
        class PatchesSeq,
        class InterfacesDerivCaculatorsSeq>
class MultipatchInterpolator;

template <
        class ExecSpace,
        class MemorySpace,
        template <typename P>
        typename BSpline1OnPatch,
        template <typename P>
        typename BSpline2OnPatch,
        template <typename P>
        typename Grid1OnPatch,
        template <typename P>
        typename Grid2OnPatch,
        ddc::SplineBuilderClosure SBCLower1,
        ddc::SplineBuilderClosure SBCUpper1,
        ddc::SplineBuilderClosure SBCLower2,
        ddc::SplineBuilderClosure SBCUpper2,
        ddc::SplineSolver Solver,
        class ExtrapolationRule,
        class PatchLocator,
        class Connectivity,
        // class InterfacesDerivCaculatorCollection,
        template <typename P>
        typename ValuesOnPatch,
        class... Patches,
        class... InterfacesDerivCaculators>
class MultipatchInterpolator<
        ExecSpace,
        MemorySpace,
        BSpline1OnPatch,
        BSpline2OnPatch,
        Grid1OnPatch,
        Grid2OnPatch,
        SBCLower1,
        SBCUpper1,
        SBCLower2,
        SBCUpper2,
        Solver,
        ExtrapolationRule,
        PatchLocator,
        Connectivity,
        ValuesOnPatch,
        ddc::detail::TypeSeq<Patches...>,
        ddc::detail::TypeSeq<InterfacesDerivCaculators...>>
{
    // TODO: check InterfacesDerivCaculatorCollection contains InterfacesDerivCaculator.
    // TODO: assert that all the grids are covered.

    using MultipatchSplineBuilder = MultipatchSplineBuilder2D<
            ExecSpace,
            MemorySpace,
            BSpline1OnPatch,
            BSpline2OnPatch,
            Grid1OnPatch,
            Grid2OnPatch,
            SBCLower1,
            SBCUpper1,
            SBCLower2,
            SBCUpper2,
            ddc::SplineBuilderClosure::HERMITE,
            Connectivity,
            Solver,
            ValuesOnPatch,
            Patches...>;

    using MultipatchSplineEvaluator = MultipatchSplineEvaluator2D<
            ExecSpace,
            MemorySpace,
            BSpline1OnPatch,
            BSpline2OnPatch,
            Grid1OnPatch,
            Grid2OnPatch,
            ExtrapolationRule,
            ValuesOnPatch,
            PatchLocator,
            Patches...>;

    using InterfacesDerivCaculatorCollection = std::tuple<InterfacesDerivCaculators const&...>;

    // InterfacesDerivativeCalculator<Connectivity, Grid1D, PatchSeq, DerivCoeffCollection>;

    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;
    using MultipatchCoordConstField = MultipatchField<ConstCoordFieldOnPatch, Patches...>;

private:
    MultipatchSplineBuilder const& m_builders;

    MultipatchSplineEvaluator const& m_evaluators;

    InterfacesDerivCaculatorCollection const& m_interfaces_deriv_calculator_collect;


public:
    MultipatchInterpolator(
            MultipatchSplineBuilder const& builders,
            MultipatchSplineEvaluator const& evaluators,
            InterfacesDerivCaculatorCollection const& interfaces_deriv_calculator_collect)
        : m_builders(builders)
        , m_evaluators(evaluators)
        , m_interfaces_deriv_calculator_collect(interfaces_deriv_calculator_collect)
    {
    }

    void operator()(
            MultipatchValues& functions_and_derivs,
            MultipatchCoordConstField const& eval_coords);
    {
        MultipatchType<IdxRangeBatchedOnPatch, Patches...> idx_ranges;
        MultipachFieldMem<SplineFieldMemOnPatch, Patches...> splines_alloc(
                m_builders.batched_spline_idx_ranges(idx_ranges));
        MultipachField<SplineFieldOnPatch, Patches...> splines(splines_alloc);

        m_builders(splines, get_const_field(functions_and_derivs));

        // Evaluate the function values at the given coordinates.
        m_evaluators(functions_and_derivs, eval_coords, get_const_field(splines));

        // Reconstructe the derivatives.
        (reconstructe_interfaces_derivatives(
                 functions_and_derivs,
                 std::get<InterfacesDerivCaculators const&...>(
                         m_interfaces_deriv_calculator_collect)),
         ...);

        (reconstructe_corner_cross_derivatives(
                 functions_and_derivs,
                 std::get<InterfacesDerivCaculators const&...>(
                         m_interfaces_deriv_calculator_collect)),
         ...);
    }

private:
    template <class InterfacesDerivCaculatorType>
    void reconstructe_interfaces_derivatives(
            MultipatchValues& functions_and_derivs,
            InterfacesDerivCaculatorType const& interfaces_deriv_calculator)
    {
        interfaces_deriv_calculator.solve_deriv(functions_and_derivs);
    }

    template <class InterfacesDerivCaculatorType>
    void reconstructe_corner_cross_derivatives(
            MultipatchValues& functions_and_derivs,
            InterfacesDerivCaculatorType const& interfaces_deriv_calculator)
    {
        // TODO: Select when it is needed to avoid duplication.
        interfaces_deriv_calculator.solve_cross_deriv(functions_and_derivs);
    }
}