// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>
#include <utility>

namespace ddc {

template <class T>
inline constexpr bool enable_borrowed_chunk = false;

template <class T>
inline constexpr bool enable_chunk = false;

template <class T>
inline constexpr bool is_chunk_v = enable_chunk<std::remove_const_t<std::remove_reference_t<T>>>;

template <class T>
inline constexpr bool is_borrowed_chunk_v
        = is_chunk_v<
                  T> && (std::is_lvalue_reference_v<T> || enable_borrowed_chunk<std::remove_cv_t<std::remove_reference_t<T>>>);

template <class T>
struct chunk_traits
{
    static_assert(is_chunk_v<T>);
    using value_type
            = std::remove_cv_t<std::remove_pointer_t<decltype(std::declval<T>().data_handle())>>;
    using pointer_type = decltype(std::declval<T>().data_handle());
    using reference_type = decltype(*std::declval<T>().data_handle());
};

template <class T>
using chunk_value_t = typename chunk_traits<T>::value_type;

template <class T>
using chunk_pointer_t = typename chunk_traits<T>::pointer_type;

template <class T>
using chunk_reference_t = typename chunk_traits<T>::reference_type;

template <class T>
inline constexpr bool is_writable_chunk_v
        = !std::is_const_v<std::remove_pointer_t<chunk_pointer_t<T>>>;

} // namespace ddc
