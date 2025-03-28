

# File pdi\_helper.hpp

[**File List**](files.md) **>** [**io**](dir_c184e51c84f2c3f0345bbc8a0d75d3e1.md) **>** [**pdi\_helper.hpp**](pdi__helper_8hpp.md)

[Go to the documentation of this file](pdi__helper_8hpp.md)


```C++
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

template <class T, class... Args>
void PDI_get_arrays(
        std::string const& event_name,
        std::string const& name,
        std::vector<T>& out_vector,
        Args&... input_args)
{
    std::integer_sequence idx_sequence = std::make_index_sequence<sizeof...(Args) / 2 + 1> {};
    std::tuple arg_tuple = std::tie(name, out_vector, input_args...);
    auto names = detail::get_name_tuple(arg_tuple, idx_sequence);
    auto out_vectors = detail::get_vector_tuple(arg_tuple, idx_sequence);
    detail::PDI_get_array(event_name, names, out_vectors, idx_sequence);
}

template <class... Grids>
void PDI_expose_idx_range(IdxRange<Grids...> index_range, std::string name)
{
    IdxStep<Grids...> extents = index_range.extents();
    Idx<Grids...> local_starts = index_range.front();
    Idx<Grids...> global_starts(Idx<Grids> {0}...);
    // TODO: Ghosts?
    IdxStep<Grids...> starts = local_starts - global_starts;
    int constexpr n_grids = sizeof...(Grids);
    std::array<ddc::DiscreteVectorElement, n_grids> starts_arr = ddc::detail::array(starts);
    std::array<ddc::DiscreteVectorElement, n_grids> extents_arr = ddc::detail::array(extents);
    std::array<std::size_t, n_grids> starts_s_arr;
    std::array<std::size_t, n_grids> extents_s_arr;
    for (int i(0); i < n_grids; ++i) {
        starts_s_arr[i] = std::size_t(starts_arr[i]);
        extents_s_arr[i] = std::size_t(extents_arr[i]);
    }
    PDI_expose((name + "_starts").c_str(), starts_s_arr.data(), PDI_OUT);
    PDI_expose((name + "_extents").c_str(), extents_s_arr.data(), PDI_OUT);
}
```


