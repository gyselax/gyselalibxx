// SPDX-License-Identifier: MIT

#pragma once

#include <utility>

#include <Kokkos_Core.hpp>

#include "ddc/chunk_span.hpp"
#include "ddc/detail/macros.hpp"
#include "ddc/detail/tagged_vector.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_vector.hpp"


namespace detail {

template <class F, class... DDims>
class ForEachKokkosLambdaAdapter
{
    template <class T>
    using index_type = std::size_t;

    F m_f;

public:
    ForEachKokkosLambdaAdapter(F const& f) : m_f(f) {}

    DDC_FORCEINLINE_FUNCTION void operator()(index_type<DDims>... ids) const
    {
        m_f(DiscreteCoordinate<DDims...>(ids...));
    }
};

template <class ExecSpace, class Functor, class DDim0>
inline void for_each_kokkos(DiscreteDomain<DDim0> const& domain, Functor const& f) noexcept
{
    Kokkos::parallel_for(
            Kokkos::RangePolicy<ExecSpace>(
                    select<DDim0>(domain).front().uid(),
                    select<DDim0>(domain).back().uid() + 1),
            ForEachKokkosLambdaAdapter<Functor, DDim0>(f));
}

template <class ExecSpace, class Functor, class DDim0, class DDim1, class... DDims>
inline void for_each_kokkos(
        DiscreteDomain<DDim0, DDim1, DDims...> const& domain,
        Functor&& f) noexcept
{
    Kokkos::Array<std::size_t, 2 + sizeof...(DDims)> const
            begin {select<DDim0>(domain).front().uid(),
                   select<DDim1>(domain).front().uid(),
                   select<DDims>(domain).front().uid()...};
    Kokkos::Array<std::size_t, 2 + sizeof...(DDims)> const
            end {(select<DDim0>(domain).back().uid() + 1),
                 (select<DDim1>(domain).back().uid() + 1),
                 (select<DDims>(domain).back().uid() + 1)...};
    Kokkos::parallel_for(
            Kokkos::MDRangePolicy<
                    ExecSpace,
                    Kokkos::Rank<
                            2 + sizeof...(DDims),
                            Kokkos::Iterate::Right,
                            Kokkos::Iterate::Right>>(begin, end),
            ForEachKokkosLambdaAdapter<Functor, DDim0, DDim1, DDims...>(f));
}

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

} // namespace detail

/// Serial execution on the host
struct serial_host_policy
{
};

/** iterates over a nD domain using the serial execution policy
 * @param[in] domain the domain over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each(
        serial_host_policy,
        DiscreteDomain<DDims...> const& domain,
        Functor&& f) noexcept
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
inline void for_each_n(
        serial_host_policy,
        DiscreteVector<DDims...> const& extent,
        Functor&& f) noexcept
{
    detail::for_each_serial<DiscreteVector<DDims...>>(
            std::array<DiscreteVectorElement, sizeof...(DDims)> {},
            std::array<DiscreteVectorElement, sizeof...(DDims)> {get<DDims>(extent) - 1 ...},
            std::forward<Functor>(f));
}

/// Parallel execution on the default device
struct parallel_host_policy
{
};

/** iterates over a nD domain using the serial execution policy
 * @param[in] domain the domain over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each(
        parallel_host_policy,
        DiscreteDomain<DDims...> const& domain,
        Functor&& f) noexcept
{
    detail::for_each_kokkos<Kokkos::DefaultHostExecutionSpace>(domain, std::forward<Functor>(f));
}

/// Kokkos parallel execution uisng MDRange policy
struct parallel_device_policy
{
};

/** iterates over a nD domain using the parallel_device_policy execution policy
 * @param[in] domain the domain over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each(
        parallel_device_policy,
        DiscreteDomain<DDims...> const& domain,
        Functor&& f) noexcept
{
    detail::for_each_kokkos<Kokkos::DefaultExecutionSpace>(domain, std::forward<Functor>(f));
}

using default_policy = serial_host_policy;

namespace policies {

inline constexpr serial_host_policy serial_host;
inline constexpr parallel_host_policy parallel_host;
inline constexpr parallel_device_policy parallel_device;

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
    for_each(chunk_span.domain(), std::forward<Functor>(f));
}
