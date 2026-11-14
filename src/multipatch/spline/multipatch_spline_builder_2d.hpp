// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <tuple>
#include <utility>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"
#include "spline_builder_2d_collection.hpp"
#include "types.hpp"

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
        ddc::SplineBuilderClosure BcTransition,
        class Connectivity,
        ddc::SplineSolver Solver,
        template <typename P>
        typename ValuesOnPatch,
        class DerivativeCalculatorCollect,
        class PatchesSeq>
class MultipatchSplineBuilder2D;


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
        ddc::SplineBuilderClosure BcTransition,
        class Connectivity,
        ddc::SplineSolver Solver,
        template <typename P>
        typename ValuesOnPatch,
        class DerivativeCalculatorCollect,
        class... Patches>
class MultipatchSplineBuilder2D<
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
        BcTransition,
        Connectivity,
        Solver,
        ValuesOnPatch,
        DerivativeCalculatorCollect,
        ddc::TypeSeq<Patches...>>
{
    using SplineBuilderCollection = SplineBuilder2DCollection<
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
            BcTransition,
            Connectivity,
            Solver,
            ValuesOnPatch,
            Patches...>;

    // using DerivativeCalculator = InterfacesDerivativeCalculator<
    //         Connectivity,
    //         Grid1D,
    //         Patches...,
    //         DerivCoeffCollection>;

    using MultipatchSplineCoeffs = typename SplineBuilderCollection::multipatch_spline_type;

    // using MultipatchDerivConstField = typename SplineBuilderCollection::multipatch_deriv_field_type;

    // template <class P>
    // using DerivConstFieldOnPatch = typename MultipatchDerivConstField::template TypeOnPatch<P>;

    template <class P>
    using DerivFieldOnPatch = DerivField<
            // std::remove_const<typename DerivConstFieldOnPatch<P>::element_type>,
            // typename DerivConstFieldOnPatch<P>::idx_range_type,
            // typename DerivConstFieldOnPatch<P>::memory_space>;
            double,
            IdxRange<
                    ddc::Deriv<typename P::Dim1>,
                    typename P::Grid1,
                    ddc::Deriv<typename P::Dim2>,
                    typename P::Grid2>,
            MemorySpace>;

    using MultipatchDerivField = MultipatchField<DerivFieldOnPatch, Patches...>;


    // template <std::size_t IOp>
    // using DerivativeCalculator
    //         = ddc::type_seq_element_t<IOp, ddc::TypeSeq<DerivativeCalculators...>>;

    // template <std::size_t IOp>
    // using Indice = ddc::type_seq_element_t<IOp, ddc::TypeSeq<Indices...>>;


    // TODO: Check that given DerivativeCalculators are InterfacesDerivativeCalculator or InterfacesExactDerivativeCalculator
    // TODO: Check that given Indices are index on the patches/matching with the correct DerivativeCalculators.

private:
    SplineBuilderCollection const& m_builders;

    DerivativeCalculatorCollect& m_deriv_calculator_collect;



public:
    MultipatchSplineBuilder2D(
            SplineBuilderCollection const& builders,
            DerivativeCalculatorCollect& deriv_calculator_collect)
        : m_builders(builders)
        , m_deriv_calculator_collect(deriv_calculator_collect) {};


    void operator()(MultipatchSplineCoeffs splines, MultipatchDerivField functions_and_derivs)
    {
        set_all_interface_derivatives(
                std::make_integer_sequence<std::size_t, DerivativeCalculatorCollect::size> {},
                functions_and_derivs);

        m_builders(splines, get_const_field(functions_and_derivs));
    };

private:
    template <std::size_t... IOp>
    void set_all_interface_derivatives(
            std::integer_sequence<std::size_t, IOp...>,
            MultipatchDerivField functions_and_derivs)
    {
        (set_interface_derivatives<IOp>(functions_and_derivs), ...);
    };

    template <std::size_t IOp>
    void set_interface_derivatives(MultipatchDerivField functions_and_derivs)
    {
        m_deriv_calculator_collect.template apply<IOp>(functions_and_derivs);
        // std::get<IOp>(m_deriv_calculator_collect)(std::get<IOp>(m_indices), functions_and_derivs);
    };
};