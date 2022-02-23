#pragma once

#include <type_traits>

#include "ddc/discrete_domain.hpp"
#include "ddc/uniform_discretization.hpp"

template <class T>
struct is_uniform_domain : std::false_type
{
};

template <class... DDims>
struct is_uniform_domain<DiscreteDomain<DDims...>>
    : std::conditional_t<
              (is_uniform_disretization_v<DDims> && ...),
              std::true_type,
              std::false_type>
{
};

template <class T>
constexpr bool is_uniform_domain_v = is_uniform_domain<T>::value;
