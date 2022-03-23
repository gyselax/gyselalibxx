// SPDX-License-Identifier: MIT

#pragma once

#include <utility>

#include "ddc/chunk_span.hpp"
#include "ddc/detail/tagged_vector.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_domain.hpp"

namespace detail {

template <class Element, class... DDims, class Functor, class... DCoords>
inline void for_each_serial(
        detail::TaggedVector<Element, DDims...> const& start,
        detail::TaggedVector<Element, DDims...> const& end,
        Functor const& f,
        DCoords const&... dcs) noexcept
{
    if constexpr (sizeof...(DCoords) == sizeof...(DDims)) {
        f(detail::TaggedVector<Element, DDims...> {dcs...});
    } else {
        using CurrentDDim = type_seq_element_t<sizeof...(DCoords), detail::TypeSeq<DDims...>>;
        for (Element ii = select<CurrentDDim>(start); ii <= select<CurrentDDim>(end); ++ii) {
            for_each_serial(start, end, f, dcs..., ii);
        }
    }
}

template <class Element, class... DDims, class Functor>
inline void for_each_omp(
        detail::TaggedVector<Element, DDims...> const& start,
        detail::TaggedVector<Element, DDims...> const& end,
        Functor&& f) noexcept
{
    using FirstDDim = type_seq_element_t<0, detail::TypeSeq<DDims...>>;
    Element const ib = select<FirstDDim>(start);
    Element const ie = select<FirstDDim>(end);
#pragma omp parallel for default(none) shared(ib, ie, start, end, f)
    for (Element ii = ib; ii <= ie; ++ii) {
        if constexpr (sizeof...(DDims) == 1) {
            f(detail::TaggedVector<Element, FirstDDim> {ii});
        } else {
            detail::for_each_serial(start, end, f, ii);
        }
    }
}

template <class X, class T, T v>
constexpr T type_constant_v = std::integral_constant<T, v>::value;

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
    detail::for_each_serial(domain.front(), domain.back(), std::forward<Functor>(f));
}

/** iterates over a nD extent using the serial execution policy
 * @param[in] extent the extent over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each_n(serial_policy, DiscreteVector<DDims...> const& extent, Functor&& f) noexcept
{
    detail::for_each_serial(
            DiscreteVector<DDims...> {detail::type_constant_v<DDims, std::ptrdiff_t, 0>...},
            DiscreteVector<DDims...> {get<DDims>(extent) - 1 ...},
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
    detail::for_each_omp(domain.front(), domain.back(), std::forward<Functor>(f));
}

/** iterates over a nD extent using the OpenMP execution policy
 * @param[in] extent the extent over which to iterate
 * @param[in] f      a functor taking an index as parameter
 */
template <class... DDims, class Functor>
inline void for_each_n(omp_policy, DiscreteVector<DDims...> const& extent, Functor&& f) noexcept
{
    detail::for_each_omp(
            DiscreteVector<DDims...> {detail::type_constant_v<DDims, std::ptrdiff_t, 0>...},
            DiscreteVector<DDims...> {get<DDims>(extent) - 1 ...},
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
