// SPDX-License-Identifier: MIT

#pragma once

#include <utility>

#include "ddc/chunk_span.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_vector.hpp"

namespace detail {

template <class RetType, class Element, std::size_t N, class Functor, class... Is>
inline void for_each_serial(
        std::array<Element, N> const& start,
        std::array<Element, N> const& end,
        Functor const& f,
        Is const&... is) noexcept
{
    static constexpr std::size_t I = sizeof...(Is);
    if constexpr (I == N) {
        f(RetType(is...));
    } else {
        for (Element ii = start[I]; ii <= end[I]; ++ii) {
            for_each_serial<RetType>(start, end, f, is..., ii);
        }
    }
}

template <class RetType, class Element, std::size_t N, class Functor>
inline void for_each_omp(
        std::array<Element, N> const& start,
        std::array<Element, N> const& end,
        Functor&& f) noexcept
{
    Element const ib = start[0];
    Element const ie = end[0];
#pragma omp parallel for default(none) shared(ib, ie, start, end, f)
    for (Element ii = ib; ii <= ie; ++ii) {
        if constexpr (N == 1) {
            f(RetType(ii));
        } else {
            detail::for_each_serial<RetType>(start, end, f, ii);
        }
    }
}

} // namespace detail

/// Serial execution
struct serial_policy
{
};

/** iterates over a nD domain using the serial execution policy
 * @param[in] domain the domain over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each(serial_policy, DiscreteDomain<DDims...> const& domain, Functor&& f) noexcept
{
    detail::for_each_serial<DiscreteCoordinate<DDims...>>(
            detail::array(domain.front()),
            detail::array(domain.back()),
            std::forward<Functor>(f));
}

/** iterates over a nD extent using the serial execution policy
 * @param[in] extent the extent over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each_n(serial_policy, DiscreteVector<DDims...> const& extent, Functor&& f) noexcept
{
    detail::for_each_serial<DiscreteVector<DDims...>>(
            std::array<DiscreteVectorElement, sizeof...(DDims)> {},
            std::array<DiscreteVectorElement, sizeof...(DDims)> {get<DDims>(extent) - 1 ...},
            std::forward<Functor>(f));
}

/// OpenMP parallel execution on the outer loop with default scheduling
struct omp_policy
{
};

/** iterates over a nD domain using the OpenMP execution policy
 * @param[in] domain the domain over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each(omp_policy, DiscreteDomain<DDims...> const& domain, Functor&& f) noexcept
{
    detail::for_each_omp<DiscreteCoordinate<DDims...>>(
            detail::array(domain.front()),
            detail::array(domain.back()),
            std::forward<Functor>(f));
}

/** iterates over a nD extent using the OpenMP execution policy
 * @param[in] extent the extent over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each_n(omp_policy, DiscreteVector<DDims...> const& extent, Functor&& f) noexcept
{
    detail::for_each_omp<DiscreteVector<DDims...>>(
            std::array<DiscreteVectorElement, sizeof...(DDims)> {},
            std::array<DiscreteVectorElement, sizeof...(DDims)> {get<DDims>(extent) - 1 ...},
            std::forward<Functor>(f));
}

using default_policy = serial_policy;

namespace policies {

inline constexpr omp_policy omp;
inline constexpr serial_policy serial;

}; // namespace policies

/** iterates over a nD domain using the default execution policy
 * @param[in] domain the domain over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each(DiscreteDomain<DDims...> const& domain, Functor&& f) noexcept
{
    for_each(default_policy(), domain, std::forward<Functor>(f));
}

/** iterates over a nD extent using the default execution policy
 * @param[in] extent the extent over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each_n(DiscreteVector<DDims...> const& extent, Functor&& f) noexcept
{
    for_each_n(default_policy(), extent, std::forward<Functor>(f));
}

template <
        class ExecutionPolicy,
        class ElementType,
        class... DDims,
        class LayoutPolicy,
        class Functor>
inline void for_each_elem(
        ExecutionPolicy&& policy,
        ChunkSpan<ElementType, DiscreteDomain<DDims...>, LayoutPolicy> chunk_span,
        Functor&& f) noexcept
{
    for_each(std::forward<ExecutionPolicy>(policy), chunk_span.domain(), std::forward<Functor>(f));
}

template <class ElementType, class... DDims, class LayoutPolicy, class Functor>
inline void for_each_elem(
        ChunkSpan<ElementType, DiscreteDomain<DDims...>, LayoutPolicy> chunk_span,
        Functor&& f) noexcept
{
    for_each(default_policy(), chunk_span.domain(), std::forward<Functor>(f));
}
