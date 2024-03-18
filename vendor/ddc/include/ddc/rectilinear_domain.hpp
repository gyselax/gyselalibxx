// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>

#include "ddc/discrete_domain.hpp"
#include "ddc/non_uniform_point_sampling.hpp"
#include "ddc/uniform_point_sampling.hpp"

namespace ddc {

template <class T>
struct is_rectilinear_domain : std::false_type
{
};

template <class... DDims>
struct is_rectilinear_domain<DiscreteDomain<DDims...>>
    : std::conditional_t<
              ((is_uniform_sampling_v<DDims> || is_non_uniform_sampling_v<DDims>)&&...),
              std::true_type,
              std::false_type>
{
};

template <class T>
constexpr bool is_rectilinear_domain_v = is_rectilinear_domain<T>::value;

} // namespace ddc
