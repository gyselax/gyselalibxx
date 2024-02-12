// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

#include <Kokkos_Core.hpp>

#include "ddc/chunk_span.hpp"
#include "ddc/detail/kokkos.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_vector.hpp"

namespace ddc {

namespace detail {

template <class F, class... DDims>
class ForEachKokkosLambdaAdapter
{
    template <class T>
    using index_type = std::size_t;

    F m_f;

public:
    ForEachKokkosLambdaAdapter(F const& f) : m_f(f) {}

    template <std::size_t N = sizeof...(DDims), std::enable_if_t<(N == 0), bool> = true>
    KOKKOS_IMPL_FORCEINLINE void operator()([[maybe_unused]] index_type<void> unused_id) const
    {
        m_f(DiscreteElement<>());
    }

    template <std::size_t N = sizeof...(DDims), std::enable_if_t<(N == 0), bool> = true>
    KOKKOS_FORCEINLINE_FUNCTION void operator()(
            use_annotated_operator,
            [[maybe_unused]] index_type<void> unused_id) const
    {
        m_f(DiscreteElement<>());
    }

    template <std::size_t N = sizeof...(DDims), std::enable_if_t<(N > 0), bool> = true>
    KOKKOS_IMPL_FORCEINLINE void operator()(index_type<DDims>... ids) const
    {
        m_f(DiscreteElement<DDims...>(ids...));
    }

    template <std::size_t N = sizeof...(DDims), std::enable_if_t<(N > 0), bool> = true>
    KOKKOS_FORCEINLINE_FUNCTION void operator()(use_annotated_operator, index_type<DDims>... ids)
            const
    {
        m_f(DiscreteElement<DDims...>(ids...));
    }
};

template <class ExecSpace, class Functor>
inline void for_each_kokkos(
        [[maybe_unused]] DiscreteDomain<> const& domain,
        Functor const& f) noexcept
{
    if constexpr (need_annotated_operator<ExecSpace>()) {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace, use_annotated_operator>(0, 1),
                ForEachKokkosLambdaAdapter<Functor>(f));
    } else {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace>(0, 1),
                ForEachKokkosLambdaAdapter<Functor>(f));
    }
}

template <class ExecSpace, class Functor, class DDim0>
inline void for_each_kokkos(DiscreteDomain<DDim0> const& domain, Functor const& f) noexcept
{
    DiscreteElement<DDim0> const ddc_begin = domain.front();
    DiscreteElement<DDim0> const ddc_end = domain.front() + domain.extents();
    std::size_t const begin = ddc::uid<DDim0>(ddc_begin);
    std::size_t const end = ddc::uid<DDim0>(ddc_end);
    if constexpr (need_annotated_operator<ExecSpace>()) {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace, use_annotated_operator>(begin, end),
                ForEachKokkosLambdaAdapter<Functor, DDim0>(f));
    } else {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace>(begin, end),
                ForEachKokkosLambdaAdapter<Functor, DDim0>(f));
    }
}

template <class ExecSpace, class Functor, class DDim0, class DDim1, class... DDims>
inline void for_each_kokkos(
        DiscreteDomain<DDim0, DDim1, DDims...> const& domain,
        Functor&& f) noexcept
{
    DiscreteElement<DDim0, DDim1, DDims...> const ddc_begin = domain.front();
    DiscreteElement<DDim0, DDim1, DDims...> const ddc_end = domain.front() + domain.extents();
    Kokkos::Array<std::size_t, 2 + sizeof...(DDims)> const
            begin {ddc::uid<DDim0>(ddc_begin),
                   ddc::uid<DDim1>(ddc_begin),
                   ddc::uid<DDims>(ddc_begin)...};
    Kokkos::Array<std::size_t, 2 + sizeof...(DDims)> const
            end {ddc::uid<DDim0>(ddc_end), ddc::uid<DDim1>(ddc_end), ddc::uid<DDims>(ddc_end)...};
    if constexpr (need_annotated_operator<ExecSpace>()) {
        Kokkos::parallel_for(
                Kokkos::MDRangePolicy<
                        ExecSpace,
                        Kokkos::Rank<
                                2 + sizeof...(DDims),
                                Kokkos::Iterate::Right,
                                Kokkos::Iterate::Right>,
                        use_annotated_operator>(begin, end),
                ForEachKokkosLambdaAdapter<Functor, DDim0, DDim1, DDims...>(f));
    } else {
        Kokkos::parallel_for(
                Kokkos::MDRangePolicy<
                        ExecSpace,
                        Kokkos::Rank<
                                2 + sizeof...(DDims),
                                Kokkos::Iterate::Right,
                                Kokkos::Iterate::Right>>(begin, end),
                ForEachKokkosLambdaAdapter<Functor, DDim0, DDim1, DDims...>(f));
    }
}

template <class RetType, class Element, std::size_t N, class Functor, class... Is>
inline void for_each_serial(
        std::array<Element, N> const& begin,
        std::array<Element, N> const& end,
        Functor const& f,
        Is const&... is) noexcept
{
    static constexpr std::size_t I = sizeof...(Is);
    if constexpr (I == N) {
        f(RetType(is...));
    } else {
        for (Element ii = begin[I]; ii < end[I]; ++ii) {
            for_each_serial<RetType>(begin, end, f, is..., ii);
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
    DiscreteElement<DDims...> const ddc_begin = domain.front();
    DiscreteElement<DDims...> const ddc_end = domain.front() + domain.extents();
    std::array const begin = detail::array(ddc_begin);
    std::array const end = detail::array(ddc_end);
    detail::for_each_serial<DiscreteElement<DDims...>>(begin, end, std::forward<Functor>(f));
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
    DiscreteVector<DDims...> const ddc_begin {};
    DiscreteVector<DDims...> const ddc_end = extent;
    std::array const begin = detail::array(ddc_begin);
    std::array const end = detail::array(ddc_end);
    detail::for_each_serial<DiscreteVector<DDims...>>(begin, end, std::forward<Functor>(f));
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

template <typename ExecSpace>
constexpr auto policy([[maybe_unused]] ExecSpace exec_space)
{
    if constexpr (std::is_same_v<ExecSpace, Kokkos::Serial>) {
        return ddc::policies::serial_host;
#ifdef KOKKOS_ENABLE_OPENMP
    } else if constexpr (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
        return ddc::policies::parallel_host;
#endif
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    } else {
        return ddc::policies::parallel_device;
#endif
    }
}

} // namespace policies

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

} // namespace ddc
