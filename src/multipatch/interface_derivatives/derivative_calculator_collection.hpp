// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "derivative_calculator_common.hpp"

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


template <class T>
inline constexpr bool enable_derivative_calculator_collection = false;

template <class T>
inline constexpr bool is_derivative_calculator_collection_v
        = enable_derivative_calculator_collection<std::remove_const_t<std::remove_reference_t<T>>>;


template <class DerivListType, class... DerivativeCalculator>
class DerivativeCalculatorCollection;

/**
 * @brief Collection of derivative calculators.
 * 
 * The DirectionDerivs can be seen as steps of the algorithm determining all 
 * the derivatives at the patch interfaces. 
 * 
 * @tparam DirectionDerivs Sorted DerivPerp and/or CrossDeriv.
 *      Their order acts as an indication on the algorithm to follow to compute all the derivatives. 
 * @tparam DerivativeCalculator Derivative calculators (non sorted).
 */
template <class... DirectionDerivs, class... DerivativeCalculators>
class DerivativeCalculatorCollection<DerivList<DirectionDerivs...>, DerivativeCalculators...>
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

public:
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
    std::tuple<DerivativeCalculators&...> m_deriv_calculator_collection;

public:
    /**
     * @brief Collect the derivative calculators and the DerivPerp and/or CrossDeriv
     * tags.
     * @param deriv_calculators Derivative calculators (non sorted).
     */
    DerivativeCalculatorCollection(
            DerivList<DirectionDerivs...>,
            DerivativeCalculators&... deriv_calculators)
        : m_deriv_calculator_collection(deriv_calculators...)
    {
    }

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

    /**
     * @brief Apply the correct derivative calculator 
     * to set the correct derivatives of the given MultipatchDerivField
     * at the step IOpth of the algorithm. 
     */
    template <std::size_t IOp, class MultipatchDerivFieldType>
    void apply(MultipatchDerivFieldType field)
    {
        this->get<IOp>()(this->get_index<IOp>(), field);
    }

private:
    template <class IdxDeriv>
    struct order_one;

    template <class... DerivDims>
    struct order_one<Idx<DerivDims...>>
    {
        static constexpr Idx<DerivDims...> value = Idx<DerivDims...>(Idx<DerivDims>(1)...);
    };
};

// To help template deduction
template <class... DirectionDerivs, class... DerivativeCalculators>
DerivativeCalculatorCollection(DerivList<DirectionDerivs...>, DerivativeCalculators&...)
        -> DerivativeCalculatorCollection<DerivList<DirectionDerivs...>, DerivativeCalculators...>;
