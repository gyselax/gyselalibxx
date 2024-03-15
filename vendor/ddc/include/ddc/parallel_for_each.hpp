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
    explicit ForEachKokkosLambdaAdapter(F const& f) : m_f(f) {}

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
        ExecSpace const& execution_space,
        [[maybe_unused]] DiscreteDomain<> const& domain,
        Functor const& f) noexcept
{
    if constexpr (need_annotated_operator<ExecSpace>()) {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace, use_annotated_operator>(execution_space, 0, 1),
                ForEachKokkosLambdaAdapter<Functor>(f));
    } else {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace>(execution_space, 0, 1),
                ForEachKokkosLambdaAdapter<Functor>(f));
    }
}

template <class ExecSpace, class Functor, class DDim0>
inline void for_each_kokkos(
        ExecSpace const& execution_space,
        DiscreteDomain<DDim0> const& domain,
        Functor const& f) noexcept
{
    DiscreteElement<DDim0> const ddc_begin = domain.front();
    DiscreteElement<DDim0> const ddc_end = domain.front() + domain.extents();
    std::size_t const begin = ddc::uid<DDim0>(ddc_begin);
    std::size_t const end = ddc::uid<DDim0>(ddc_end);
    if constexpr (need_annotated_operator<ExecSpace>()) {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace, use_annotated_operator>(execution_space, begin, end),
                ForEachKokkosLambdaAdapter<Functor, DDim0>(f));
    } else {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace>(execution_space, begin, end),
                ForEachKokkosLambdaAdapter<Functor, DDim0>(f));
    }
}

template <class ExecSpace, class Functor, class DDim0, class DDim1, class... DDims>
inline void for_each_kokkos(
        ExecSpace const& execution_space,
        DiscreteDomain<DDim0, DDim1, DDims...> const& domain,
        Functor const& f) noexcept
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
                        use_annotated_operator>(execution_space, begin, end),
                ForEachKokkosLambdaAdapter<Functor, DDim0, DDim1, DDims...>(f));
    } else {
        Kokkos::parallel_for(
                Kokkos::MDRangePolicy<
                        ExecSpace,
                        Kokkos::Rank<
                                2 + sizeof...(DDims),
                                Kokkos::Iterate::Right,
                                Kokkos::Iterate::Right>>(execution_space, begin, end),
                ForEachKokkosLambdaAdapter<Functor, DDim0, DDim1, DDims...>(f));
    }
}

} // namespace detail

/** iterates over a nD domain using a given `Kokkos` execution space
 * @param[in] execution_space a Kokkos execution space where the loop will be executed on
 * @param[in] domain the domain over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class ExecSpace, class... DDims, class Functor>
inline void parallel_for_each(
        ExecSpace const& execution_space,
        DiscreteDomain<DDims...> const& domain,
        Functor&& f) noexcept
{
    detail::for_each_kokkos(execution_space, domain, std::forward<Functor>(f));
}

/** iterates over a nD domain using the `Kokkos` default execution space
 * @param[in] domain the domain over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void parallel_for_each(DiscreteDomain<DDims...> const& domain, Functor&& f) noexcept
{
    parallel_for_each(Kokkos::DefaultExecutionSpace(), domain, std::forward<Functor>(f));
}

} // namespace ddc
