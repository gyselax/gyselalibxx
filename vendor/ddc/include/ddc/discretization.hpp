#pragma once

#include <optional>

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

} // namespace detail


template <class D, class... Args>
void init_discretization(Args&&... a)
{
    detail::Discretization<D>::s_disc.emplace(std::forward<Args>(a)...);
}

template <class T>
T const& discretization()
{
    return *detail::Discretization<T>::s_disc;
}
