// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "derivative_calculator_common.hpp"


template <class T>
inline constexpr bool enable_derivative_calculator_collection = false;

template <class T>
inline constexpr bool is_derivative_calculator_collection_v
        = enable_derivative_calculator_collection<std::remove_const_t<std::remove_reference_t<T>>>;


template <class DerivativeCalculatorSeq, class IndiceSeq>
class DerivativeCalculatorCollection;

template <class... DerivativeCalculators, class... Indices>
class DerivativeCalculatorCollection<
        ddc::TypeSeq<DerivativeCalculators...>,
        ddc::TypeSeq<Indices...>>
{
    static_assert(
            sizeof...(DerivativeCalculators) == sizeof...(Indices),
            "Please provide as many indices as derivative calculators.");


public:
    template <std::size_t IOp>
    using derivative_calculator_type
            = ddc::type_seq_element_t<IOp, ddc::TypeSeq<DerivativeCalculators...>>;

    template <std::size_t IOp>
    using index_type = ddc::type_seq_element_t<IOp, ddc::TypeSeq<Indices...>>;

    static constexpr std::size_t size = sizeof...(Indices);

private:
    static_assert(
            ((is_deriv_calculator_v<DerivativeCalculators>)&&...),
            "Please provide derivative calculators.");

    static_assert((ddc::is_discrete_element_v<Indices> && ...), "Please provide Idx types.");

    using IndexSeq = ddc::TypeSeq<Indices...>;

    template <std::size_t IOp>
    struct index_matches_calculator
    {
        using patch_test = typename derivative_calculator_type<IOp>::patch_direction;
        using idx_admissible_seq = ddc::TypeSeq<
                typename derivative_calculator_type<IOp>::idx_deriv_perp,
                Idx<ddc::Deriv<typename patch_test::Dim1>, ddc::Deriv<typename patch_test::Dim2>>,
                Idx<ddc::Deriv<typename patch_test::Dim2>, ddc::Deriv<typename patch_test::Dim1>>>;

        static constexpr bool value
                = (ddc::type_seq_contains_v<ddc::TypeSeq<index_type<IOp>>, idx_admissible_seq>);
    };

    static_assert(
            (index_matches_calculator<ddc::type_seq_rank_v<Indices, IndexSeq>>::value && ...),
            "One of the indice types do not match with the operator.");

private:
    std::tuple<DerivativeCalculators&...> m_deriv_calculator_collection;

    std::tuple<Indices...> m_indice_collection;

public:
    DerivativeCalculatorCollection(
            std::tuple<DerivativeCalculators&...> deriv_calculators,
            std::tuple<Indices...> indices)
        : m_deriv_calculator_collection(deriv_calculators)
        , m_indice_collection(indices)
    {
    }

    template <std::size_t IOp>
    index_type<IOp>& get_index() const
    {
        static_assert((IOp < size), "Index outside of the range.");
        return std::get<IOp>(m_indice_collection);
    }


    template <std::size_t IOp>
    derivative_calculator_type<IOp>& get()
    {
        static_assert((IOp < size), "Index outside of the range.");
        return std::get<IOp>(m_deriv_calculator_collection);
    }

    template <std::size_t IOp, class FieldType>
    void apply(FieldType field)
    {
        this->get<IOp>()(this->get_index<IOp>(), field);
    }
};


// To help the template deduction.
template <class... DerivativeCalculators, class... Indices>
DerivativeCalculatorCollection(std::tuple<DerivativeCalculators&...>, std::tuple<Indices...>)
        -> DerivativeCalculatorCollection<
                ddc::TypeSeq<DerivativeCalculators...>,
                ddc::TypeSeq<Indices...>>;
