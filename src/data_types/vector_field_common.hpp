// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "vector_index_tools.hpp"

template <class T>
inline constexpr bool enable_vector_field = false;

template <class T>
inline constexpr bool enable_borrowed_vector_field = false;

template <class T>
inline constexpr bool is_vector_field_v
        = enable_vector_field<std::remove_const_t<std::remove_reference_t<T>>>;

template <class T>
inline constexpr bool is_borrowed_vector_field_v
        = is_vector_field_v<
                  T> && (std::is_lvalue_reference_v<T> || enable_borrowed_vector_field<std::remove_cv_t<std::remove_reference_t<T>>>);


namespace detail {
template<class VectorIndexSetType>
class VecIdxSetWrapper {
    static_assert(is_vector_index_set_v<VectorIndexSetType>);
    using my_tag = VectorIndexSetType;
};

template <class VectorIndexSetType>
constexpr IdxRange<VecIdxSetWrapper<VectorIndexSetType>> make_idx_range()
{
    Idx<VecIdxSetWrapper<VectorIndexSetType>> start(0);
    IdxStep<VecIdxSetWrapper<VectorIndexSetType>> len(ddc::type_seq_size_v<VectorIndexSetType>);
    return IdxRange<VecIdxSetWrapper<VectorIndexSetType>>(start, len);
}

} // namespace detail

