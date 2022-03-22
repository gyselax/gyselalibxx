// SPDX-License-Identifier: MIT

#pragma once

#include <optional>
#include <stdexcept>

#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_domain.hpp"

namespace detail {

// For now, in the future, this should be specialized by tag
template <class IDim>
struct Discretization
{
    static std::optional<IDim> s_disc;
};

template <class IDim>
std::optional<IDim> Discretization<IDim>::s_disc;

template <class Tuple, std::size_t... Ids>
auto extract_after(Tuple&& t, std::index_sequence<Ids...>)
{
    return std::make_tuple(std::move(std::get<Ids + 1>(t))...);
}

} // namespace detail

template <class D, class Arg>
Arg init_discretization(std::tuple<D, Arg>&& a)
{
    if (detail::Discretization<std::remove_cv_t<std::remove_reference_t<D>>>::s_disc) {
        throw std::runtime_error("Discretization function already initialized.");
    }
    detail::Discretization<std::remove_cv_t<std::remove_reference_t<D>>>::s_disc.emplace(
            std::move(std::get<0>(a)));
    return std::get<1>(a);
}

template <class D, class... Args>
std::enable_if_t<2 <= sizeof...(Args), std::tuple<Args...>> init_discretization(
        std::tuple<D, Args...>&& a)
{
    if (detail::Discretization<std::remove_cv_t<std::remove_reference_t<D>>>::s_disc) {
        throw std::runtime_error("Discretization function already initialized.");
    }
    detail::Discretization<std::remove_cv_t<std::remove_reference_t<D>>>::s_disc.emplace(
            std::move(std::get<0>(a)));
    return detail::extract_after(std::move(a), std::index_sequence_for<Args...>());
}

template <class D, class... Args>
void init_discretization(Args&&... a)
{
    if (detail::Discretization<std::remove_cv_t<std::remove_reference_t<D>>>::s_disc) {
        throw std::runtime_error("Discretization function already initialized.");
    }
    detail::Discretization<std::remove_cv_t<std::remove_reference_t<D>>>::s_disc.emplace(
            std::forward<Args>(a)...);
}

template <class T>
T const& discretization()
{
    return *detail::Discretization<T>::s_disc;
}
