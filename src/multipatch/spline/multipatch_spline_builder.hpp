// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <tuple>
#include <utility>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


/**
 * @brief A class to call all the builders of all the patches once. 
 * 
 * We need to instantiate all the builders for all the pacthes in the main code.
 * We process the same way for the Field containing the spline coefficients and the 
 * values of the function on each patch, and we store them in std::tuple. 
 * This class is instantiated with all the builders. 
 * The operator() allows to call all the builders stored in the member of this class in 
 * one single line. 
 * 
 * This function is useful to avoid calling individually all the builders, especially in 
 * multipatch geometry with several patches. 
 *  
 * @tparam Builders Type of builders for each patch. 
 */
template <class... Builders>
class MultipatchSplineBuilder
{
    using BuilderTuple = std::tuple<Builders const&...>;
    using SplineTuple = std::tuple<
            DField<typename Builders::batched_spline_domain_type,
                   std::experimental::layout_right,
                   typename Builders::memory_space>...>;
    // For PERIODIC or GREVILLE boundary conditions
    using ValuesTuple = std::tuple<
            DField<typename Builders::batched_interpolation_domain_type,
                   std::experimental::layout_right,
                   typename Builders::memory_space>...>;


    BuilderTuple const m_builders;

public:
    /**
     * @brief Instantiate the MultipatchSplineBuilder from a std::tuple 
     * of all the builder on each patch. 
     * 
     * @warning The builders have to be sorted in the same order as the patches
     * in the tuple. 
     * 
     * @param builders Spline builders for each patch.
     */
    MultipatchSplineBuilder(Builders const&... builders) : m_builders(std::tie(builders...)) {};


    ~MultipatchSplineBuilder() = default;

    /**
     * @brief Build the spline representation of each given function.
     * 
     * @warning The splines and the values std::tuples have to be sorted 
     * in the same order as the patches in the tuple. 
     * 
     * @param splines Tuple of all the Fields pointing to the spline representations. 
     * @param values Tuple of all the Fields pointing to the function values. 
     */
    void operator()(SplineTuple const& splines, ValuesTuple const& values)
    {
        apply_builders(splines, values, std::make_index_sequence<sizeof...(Builders)> {});
    };

private:
    /// @brief Apply the builders to the splines and values elements.
    template <std::size_t... I>
    void apply_builders(
            SplineTuple const& splines,
            ValuesTuple const& values,
            std::index_sequence<I...>)
    {
        (std::get<I>(m_builders)(std::get<I>(splines), get_const_field(std::get<I>(values))), ...);
    };
};
