// SPDX-License-Identifier: MIT
#pragma once

#include <type_traits>

namespace detail {
template <class ExecSpace, class Type>
struct MappingAccessibility : std::false_type
{
};
} // namespace detail

template <class ExecSpace, class Type>
static constexpr bool is_accessible_v = detail::MappingAccessibility<ExecSpace, Type>::value;
