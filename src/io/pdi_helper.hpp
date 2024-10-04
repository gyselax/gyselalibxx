// SPDX-License-Identifier: MIT
#pragma once

#include <pdi.h>

namespace detail {

template <class NameTuple, class OutputTuple, size_t... I>
void PDI_get_array(
        std::string event_name,
        NameTuple names,
        OutputTuple args,
        std::integer_sequence<size_t, I...>)
{
    std::array<size_t, sizeof...(I)> sizes;
    // Put sizes into PDI
    ((PDI_share((std::get<I>(names) + "_extents").c_str(), &sizes[I], PDI_INOUT)), ...);
    // Collect sizes from file
    PDI_event((event_name + "_extents").c_str());
    // Collect sizes from PDI
    ((PDI_reclaim((std::get<I>(names) + "_extents").c_str())), ...);

    // Update the size of the vectors
    ((std::get<I>(args).resize(sizes[I])), ...);
    // Put vector into PDI
    ((PDI_share(std::get<I>(names).c_str(), std::get<I>(args).data(), PDI_INOUT)), ...);
    // Collect vector from file
    PDI_event(event_name.c_str());
    // Collect vector from PDI
    ((PDI_reclaim(std::get<I>(names).c_str())), ...);
}

template <class TupleType, size_t... I>
auto get_name_tuple(TupleType input_args, std::integer_sequence<size_t, I...>)
{
    return std::make_tuple(std::string(std::get<I * 2>(input_args))...);
}

template <class TupleType, size_t... I>
auto get_vector_tuple(TupleType input_args, std::integer_sequence<size_t, I...>)
{
    return std::tie(std::get<I * 2 + 1>(input_args)...);
}

} // namespace detail

/**
 * A helper function to read an unknown number of arrays from a file using PDI.
 */
template <class T, class... Args>
void PDI_get_arrays(
        std::string event_name,
        std::string name,
        std::vector<T>& out_vector,
        Args&... input_args)
{
    auto idx_sequence = std::make_index_sequence<sizeof...(Args) / 2 + 1> {};
    auto arg_tuple = std::tie(name, out_vector, input_args...);
    auto names = detail::get_name_tuple(arg_tuple, idx_sequence);
    auto out_vectors = detail::get_vector_tuple(arg_tuple, idx_sequence);
    detail::PDI_get_array(event_name, names, out_vectors, idx_sequence);
}
