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


/**
 * @brief A class to set all the needed interface derivatives and 
 * to call all the builders of all the patches once.
 *
 * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
 * @tparam MemorySpace The space (CPU/GPU) where the coefficients and values are stored.
 * @tparam BSpline1OnPatch A type alias which provides the first BSpline type along which the splines are built.
 * @tparam BSpline2OnPatch A type alias which provides the second BSpline type along which the splines are built.
 * @tparam Grid1OnPatch A type alias which provides the first Grid type along which the interpolation points of the splines are found.
 * @tparam Grid2OnPatch A type alias which provides the second Grid type along which the interpolation points of the splines are found.
 * @tparam SBCLower1 The lower spline closure on the first dimension.
 * @tparam SBCUpper1 The upper spline closure on the first dimension.
 * @tparam SBCLower2 The lower spline closure on the second dimension.
 * @tparam SBCUpper2 The upper spline closure on the second dimension.
 * @tparam BcTransition The boundary condition used at the interface between 2 patches.
 * @tparam Connectivity A MultipatchConnectivity object describing the interfaces between patches.
 * @tparam Solver The SplineSolver giving the backend used to perform the spline approximation. See DDC for more details.
 * @tparam ValuesOnPatch A type alias which provides the field type which will be used to pass the values of the function
 *                      at the interpolation points. The index range of this field type should contain any batch dimensions.
 * @tparam DerivListType A DerivList type that stores the sorted DerivPerp and/or CrossDeriv.
 *                      Their order acts as an indication on the algorithm to follow to compute all the derivatives. 
 * @tparam DerivativeCalculatorSeq  A TypeSeq of DerivativeCalculator derivative calculators (non sorted).
 * @tparam PatchesSeq A TypeSeq of patches.
 */
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

/**
 * @brief A class to set all the needed interface derivatives and 
 * to call all the builders of all the patches once.
 *
 * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
 * @tparam MemorySpace The space (CPU/GPU) where the coefficients and values are stored.
 * @tparam BSpline1OnPatch A type alias which provides the first BSpline type along which the splines are built.
 * @tparam BSpline2OnPatch A type alias which provides the second BSpline type along which the splines are built.
 * @tparam Grid1OnPatch A type alias which provides the first Grid type along which the interpolation points of the splines are found.
 * @tparam Grid2OnPatch A type alias which provides the second Grid type along which the interpolation points of the splines are found.
 * @tparam SBCLower1 The lower spline closure on the first dimension.
 * @tparam SBCUpper1 The upper spline closure on the first dimension.
 * @tparam SBCLower2 The lower spline closure on the second dimension.
 * @tparam SBCUpper2 The upper spline closure on the second dimension.
 * @tparam BcTransition The boundary condition used at the interface between 2 patches.
 * @tparam Connectivity A MultipatchConnectivity object describing the interfaces between patches.
 * @tparam Solver The SplineSolver giving the backend used to perform the spline approximation. See DDC for more details.
 * @tparam ValuesOnPatch A type alias which provides the field type which will be used to pass the values of the function
 *                      at the interpolation points. The index range of this field type should contain any batch dimensions.
 * @tparam DirectionDerivs Sorted DerivPerp and/or CrossDeriv.
 *                      Their order acts as an indication on the algorithm to follow to compute all the derivatives. 
 * @tparam DerivativeCalculators  Derivative calculators (non sorted).
 * @tparam PatchesSeq Tags for patches.
 */
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

    using MultipatchIdxRangeBS = typename SplineBuilderCollection::coeff_idx_range_type;

    template <class Patch, class BatchedIdxRange>
    using batched_spline_domain_type_on_patch
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<BatchedIdxRange>,
                    ddc::TypeSeq<typename Patch::Grid1, typename Patch::Grid2>,
                    ddc::TypeSeq<typename Patch::BSplines1, typename Patch::BSplines2>>>;

    template <class MultipatchIdxRangeBatched>
    using MultipatchIdxRangeBatchedBS = detail::MultipatchType<
            typename MultipatchIdxRangeBatched::PatchOrdering,
            ddc::TypeSeq<batched_spline_domain_type_on_patch<
                    Patches,
                    typename MultipatchIdxRangeBatched::template TypeOnPatch<Patches>>...>>;


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
    SplineBuilderCollection const m_builders_collection;

    DerivativeCalculatorCollection m_deriv_calculator_collection;

public:
    /**
     * @brief Instantiate a MultipatchSplineBuilder2D.
     * 
     * @param builders_collection A SplineBuilder2DCollection type that collects the spline builders.
     * @param deriv_calculators DerivativeCalculatorCommon types that compute the interface derivatives 
     *      in a given direction.
     */
    MultipatchSplineBuilder2D(
            SplineBuilderCollection const& builders_collection,
            DerivativeCalculators&... deriv_calculators)
        : m_builders_collection(builders_collection)
        , m_deriv_calculator_collection(deriv_calculators...) {};

    /**
     * @brief Instantiate a MultipatchSplineBuilder2D.
     * 
     * @param builders A SplineBuilder2D types.
     * @param deriv_calculators DerivativeCalculatorCommon types that compute the interface derivatives 
     *      in a given direction.
     */
    MultipatchSplineBuilder2D(
            typename SplineBuilderCollection::template builder_type_on_patch<
                    Patches> const&... builders,
            DerivativeCalculators&... deriv_calculators)
        : m_builders_collection(builders...)
        , m_deriv_calculator_collection(deriv_calculators...) {};


    /**
     * @brief Get the spline coefficients of the given function. 
     * The interface derivatives in the given DerivFields are computed using the derivative calculators
     * of the class and the builders are called to build the splines. 
     * @param[out] splines A MultipatchField of splines on each patch.
     * @param[in, out] functions_and_derivs A MultipatchField of DerivFields on each patch.
     *      The interface derivatives are modified. 
     */
    void operator()(MultipatchSplineCoeffs splines, MultipatchDerivField functions_and_derivs)
    {
        set_all_interface_derivatives(std::make_index_sequence<size> {}, functions_and_derivs);
        m_builders_collection(splines, get_const_field(functions_and_derivs));
    };


    /**
     * @brief Get a MultipatchType collecting the index ranges typed on the B-splines dimensions
     * of the patches. 
     * These index ranges are adapted to allocate splines on the multi-patch domain. 
     * @return a MultipatchType collecting the index ranges typed on the B-splines dimensions
     * of the patches. 
     */
    MultipatchIdxRangeBS const spline_idx_ranges() const
    {
        return m_builders_collection.spline_idx_ranges();
    };

    /**
     * @brief Get a MultipatchType collecting the index ranges typed on the whole B-splines dimensions.
     * These index ranges are adapted to allocate splines on the multi-patch domain. 
     * @param[in] idx_ranges MultipatchType collecting the index ranges on the grids of the patches. 
     * @return a MultipatchType collecting the index ranges typed on the whole B-splines dimensions.
     */
    template <class MultipatchIdxRangeBatched>
    MultipatchIdxRangeBatchedBS<MultipatchIdxRangeBatched> const batched_basis_idx_ranges(
            MultipatchIdxRangeBatched idx_ranges) const
    {
        return m_builders_collection.batched_basis_idx_ranges(idx_ranges);
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