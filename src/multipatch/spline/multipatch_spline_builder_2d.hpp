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

/**
 * @brief Tag to indicate to compute a perpendicular derivative with 
 * a derivative calculator.
 */
template <class Grid1D>
struct DerivPerp
{
    /// @brief Grid 1D given to the class indicating the direction of the mesh line along which the derivatives are computed.
    using grid_direction = Grid1D;
    /// @brief Index of the perpendicular derivative.
    template <class Calculator>
    using idx_deriv = typename Calculator::idx_deriv_perp;
};

/**
 * @brief Tag to indicate to compute a cross-derivative with 
 * a derivative calculator.
 */
template <class Grid1D>
struct CrossDeriv
{
    /// @brief Grid 1D given to the class indicating the direction of the mesh line along which the derivatives are computed.
    using grid_direction = Grid1D;
    /// @brief Index of the cross-derivative.
    template <class Calculator>
    using idx_deriv = typename Calculator::idx_cross_deriv;
};

/**
 * @brief Sorted list of DerivPerp and/or CrossDeriv.
 */
template <class... DirectionDerivs>
struct DerivList
{
};



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
        class DerivListType,
        class DerivativeCalculatorSeq,
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
        class... DirectionDerivs,
        class... DerivativeCalculators,
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
        DerivList<DirectionDerivs...>,
        ddc::TypeSeq<DerivativeCalculators...>,
        ddc::TypeSeq<Patches...>>
{
    static_assert(
            sizeof...(DerivativeCalculators) <= sizeof...(DirectionDerivs),
            "Please provide enough indices for the derivative calculators.");

    static_assert(
            ((is_deriv_calculator_v<DerivativeCalculators>)&&...),
            "Please provide derivative calculators.");

    template <class Grid1D>
    static constexpr std::size_t number_of_calculators_on
            = (std::size_t(0) + ...
               + std::is_same_v<Grid1D, typename DerivativeCalculators::grid_direction>);
    static_assert(
            ((number_of_calculators_on<typename DerivativeCalculators::grid_direction> == 1)
             && ...),
            "Two derivative calculators share the same grid_direction.");

    using grid_calculator_seq = ddc::TypeSeq<typename DerivativeCalculators::grid_direction...>;

    static_assert(
            (ddc::type_seq_contains_v<
                     ddc::TypeSeq<typename DirectionDerivs::grid_direction>,
                     grid_calculator_seq> && ...),
            "No derivative calculator found for one of the DirectionDeriv grid_direction.");


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

    // InterfacesDerivativeCalculator<Connectivity, GridX<1>

    using MultipatchSplineCoeffs = typename SplineBuilderCollection::multipatch_spline_type;

    template <class P>
    using DerivFieldOnPatch = DerivField<
            double,
            IdxRange<
                    ddc::Deriv<typename P::Dim1>,
                    typename P::Grid1,
                    ddc::Deriv<typename P::Dim2>,
                    typename P::Grid2>,
            MemorySpace>;

    using MultipatchDerivField = MultipatchField<DerivFieldOnPatch, Patches...>;

    using DerivativeCalculatorCollection = std::tuple<DerivativeCalculators&...>;

private:
    /// @brief Sorted list of DerivPerp and/or CrossDeriv.
    using direction_deriv_list = DerivList<DirectionDerivs...>;

    /// @brief Type sequence of DerivPerp and/or CrossDeriv.
    using direction_deriv_seq = ddc::TypeSeq<DirectionDerivs...>;

    /// @brief Number of steps in the algorithm (i.e. of given derivative indices).
    static constexpr std::size_t size = sizeof...(DirectionDerivs);

    /// @brief Type sequence of derivative calculators.
    using derivative_calculator_seq = ddc::TypeSeq<DerivativeCalculators...>;

    /// @brief 1D grid indicating the direction to follow at the step IOpth.
    template <std::size_t IOp>
    using grid_direction =
            typename ddc::type_seq_element_t<IOp, direction_deriv_seq>::grid_direction;

    /// @brief Derivative calculator to use at the step IOpth.
    template <std::size_t IOp>
    using derivative_calculator_type = ddc::type_seq_element_t<
            ddc::type_seq_rank_v<grid_direction<IOp>, grid_calculator_seq>,
            derivative_calculator_seq>;

    /// @brief Derivative index to apply at the step IOpth.
    template <std::size_t IOp>
    using index_type =
            typename ddc::type_seq_element_t<IOp, direction_deriv_seq>::template idx_deriv<
                    derivative_calculator_type<IOp>>;

private:
    SplineBuilderCollection const& m_builders_collection;

    DerivativeCalculatorCollection m_deriv_calculator_collection;

public:
    MultipatchSplineBuilder2D(
            SplineBuilderCollection const& builders_collection,
            DerivativeCalculators&... deriv_calculators)
        : m_builders_collection(builders_collection)
        , m_deriv_calculator_collection(deriv_calculators...) {};


    void operator()(MultipatchSplineCoeffs splines, MultipatchDerivField functions_and_derivs)
    {
        set_all_interface_derivatives(std::make_index_sequence<size> {}, functions_and_derivs);

        m_builders_collection(splines, get_const_field(functions_and_derivs));
    };

private:
    template <std::size_t... IOp>
    void set_all_interface_derivatives(
            std::integer_sequence<std::size_t, IOp...>,
            MultipatchDerivField functions_and_derivs)
    {
        (set_interface_derivatives<IOp>(functions_and_derivs), ...);
    };

     /**
     * @brief Apply the correct derivative calculator 
     * to set the correct derivatives of the given MultipatchDerivField
     * at the step IOpth of the algorithm. 
     */
    template <std::size_t IOp>
    void set_interface_derivatives(MultipatchDerivField functions_and_derivs)
    {
        this->get<IOp>()(this->get_index<IOp>(), functions_and_derivs);
    };

    /**
     * @brief Get the derivative index at the step IOpth of the algorithm. 
     */
    template <std::size_t IOp>
    static index_type<IOp> get_index()
    {
        static_assert((IOp < size), "Index outside of the range.");
        return order_one<index_type<IOp>>::value;
    }

    /**
     * @brief Get a reference to the derivative calculator needed
     * at the step IOpth of the algorithm. 
     */
    template <std::size_t IOp>
    derivative_calculator_type<IOp>& get()
    {
        static_assert((IOp < size), "Index outside of the range.");
        return std::get<derivative_calculator_type<IOp>&>(m_deriv_calculator_collection);
    }

    template <class IdxDeriv>
    struct order_one;

    template <class... DerivDims>
    struct order_one<Idx<DerivDims...>>
    {
        static constexpr Idx<DerivDims...> value = Idx<DerivDims...>(Idx<DerivDims>(1)...);
    };
};