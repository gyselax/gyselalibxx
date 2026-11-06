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
        class... DerivativeCalculators,
        class... Patches>
class MultipatchSplineBuilder2D
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
            DerivCoeffCollection,
            Patches...>;

    // using DerivativeCalculator = InterfacesDerivativeCalculator<
    //         Connectivity,
    //         Grid1D,
    //         Patches...,
    //         DerivCoeffCollection>;

    using MultipatchSplineCoeffs = typename SplineBuilderCollection::multipatch_spline_type;

    using MultipatchDerivConstField = typename SplineBuilderCollection::multipatch_deriv_field_type;

    template <class P>
    using DerivConstFieldOnPatch = typename MultipatchDerivConstField::template TypeOnPatch<P>;

    template <class P>
    using DerivFieldOnPatch = DerivField<
            std::remove_const<typename DerivConstFieldOnPatch<P>::element_type>,
            typename DerivConstFieldOnPatch<P>::idx_range_type,
            typename DerivConstFieldOnPatch<P>::memory_space>;

    using MultipatchDerivField = MultipatchField<DerivFieldOnPatch, Patches...>;


private:
    SplineBuilderCollection const& m_builders;

public:
    MultipatchSplineBuilder2D(SplineBuilderCollection const& builders) : m_builders(builders) {};


    void operator()(MultipatchSplineCoeffs splines, MultipatchDerivField fonctions_and_derivs)
    {
        // TODO: compute the interfaces derivatives.

        deriv_calculator.solve_deriv(fonctions_and_derivs); 
        // ...

        deriv_calculator.solve_cross_deriv(fonctions_and_derivs); 
        // ... 

        m_builders(splines, get_const_field(fonctions_and_derivs));
    };
};