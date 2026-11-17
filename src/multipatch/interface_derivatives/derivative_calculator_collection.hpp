// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "derivative_calculator_common.hpp"

template <class Grid1D>
struct DerivPerp
{
    using grid_direction = Grid1D;
    template <class Calculator>
    using idx_deriv = typename Calculator::idx_deriv_perp;
};

template <class Grid1D>
struct CrossDeriv
{
    using grid_direction = Grid1D;
    template <class Calculator>
    using idx_deriv = typename Calculator::idx_cross_deriv;
};

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
    using direction_deriv_list = DerivList<DirectionDerivs...>;

    using direction_deriv_seq = ddc::TypeSeq<DirectionDerivs...>;
    static constexpr std::size_t size = sizeof...(DirectionDerivs);

    using derivative_calculator_seq = ddc::TypeSeq<DerivativeCalculators...>;

    template <std::size_t IOp>
    using grid_direction =
            typename ddc::type_seq_element_t<IOp, direction_deriv_seq>::grid_direction;

    template <std::size_t IOp>
    using derivative_calculator_type = ddc::type_seq_element_t<
            ddc::type_seq_rank_v<grid_direction<IOp>, grid_calculator_seq>,
            derivative_calculator_seq>;

    template <std::size_t IOp>
    using index_type =
            typename ddc::type_seq_element_t<IOp, direction_deriv_seq>::template idx_deriv<
                    derivative_calculator_type<IOp>>;

private:
    std::tuple<DerivativeCalculators&...> m_deriv_calculator_collection;

public:
    DerivativeCalculatorCollection(
            DerivList<DirectionDerivs...>,
            DerivativeCalculators&... deriv_calculators)
        : m_deriv_calculator_collection(deriv_calculators...)
    {
    }

    template <std::size_t IOp>
    index_type<IOp> get_index() const
    {
        static_assert((IOp < size), "Index outside of the range.");
        return order_one<index_type<IOp>>::value;
    }

    template <std::size_t IOp>
    derivative_calculator_type<IOp>& get()
    {
        static_assert((IOp < size), "Index outside of the range.");
        return std::get<derivative_calculator_type<IOp>&>(m_deriv_calculator_collection);
    }

    template <std::size_t IOp, class FieldType>
    void apply(FieldType field)
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

template <class... DirectionDerivs, class... DerivativeCalculators>
DerivativeCalculatorCollection(DerivList<DirectionDerivs...>, DerivativeCalculators&...)
        -> DerivativeCalculatorCollection<DerivList<DirectionDerivs...>, DerivativeCalculators...>;


// /**
//  * @brief Helper to build a DerivativeCalculatorCollection.
//  * The DirectionDerivs are given explicitly, the DerivativeCalculators are deduced.
//  * (CTAD cannot be used as it does not work with partially specified template arguments.)
//  *
//  * Usage: make_derivative_calculator_collection<DirectionDeriv<...>, ...>(calc_1, calc_2, ...);
//  */
// template <class... DirectionDerivs, class... DerivativeCalculators>
// DerivativeCalculatorCollection<
//         ddc::TypeSeq<DirectionDerivs...>,
//         ddc::TypeSeq<DerivativeCalculators...>>
// make_derivative_calculator_collection(DerivativeCalculators&... deriv_calculators)
// {
//     return DerivativeCalculatorCollection<
//             ddc::TypeSeq<DirectionDerivs...>,
//             ddc::TypeSeq<DerivativeCalculators...>>(std::tie(deriv_calculators...));
// }
