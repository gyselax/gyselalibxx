// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <type_traits>
#include <utility>

#include "ddc/coordinate.hpp"
#include "ddc/detail/macros.hpp"
#include "ddc/detail/type_seq.hpp"
#include "ddc/discrete_vector.hpp"

namespace ddc {

template <class...>
class DiscreteElement;

template <class T>
struct IsDiscreteElement : std::false_type
{
};

template <class... Tags>
struct IsDiscreteElement<DiscreteElement<Tags...>> : std::true_type
{
};

template <class T>
inline constexpr bool is_discrete_element_v = IsDiscreteElement<T>::value;


namespace detail {

template <class... Tags>
struct ToTypeSeq<DiscreteElement<Tags...>>
{
    using type = TypeSeq<Tags...>;
};

} // namespace detail

/** A DiscreteCoordElement is a scalar that identifies an element of the discrete dimension
 */
using DiscreteElementType = std::size_t;

template <class Tag>
KOKKOS_FUNCTION constexpr DiscreteElementType const& uid(DiscreteElement<Tag> const& tuple) noexcept
{
    return tuple.uid();
}

template <class Tag>
KOKKOS_FUNCTION constexpr DiscreteElementType& uid(DiscreteElement<Tag>& tuple) noexcept
{
    return tuple.uid();
}

template <class QueryTag, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteElementType const& uid(
        DiscreteElement<Tags...> const& tuple) noexcept
{
    return tuple.template uid<QueryTag>();
}

template <class QueryTag, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteElementType& uid(DiscreteElement<Tags...>& tuple) noexcept
{
    return tuple.template uid<QueryTag>();
}

template <class QueryTag, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteElementType const& uid_or(
        DiscreteElement<Tags...> const& tuple,
        DiscreteElementType const& default_value) noexcept
{
    return tuple.template uid_or<QueryTag>(default_value);
}

template <class... QueryTags, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteElement<QueryTags...> select(
        DiscreteElement<Tags...> const& arr) noexcept
{
    return DiscreteElement<QueryTags...>(arr);
}

template <class... QueryTags, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteElement<QueryTags...> select(
        DiscreteElement<Tags...>&& arr) noexcept
{
    return DiscreteElement<QueryTags...>(std::move(arr));
}

/// Returns a reference towards the DiscreteElement that contains the QueryTag
template <
        class QueryTag,
        class HeadDElem,
        class... TailDElems,
        std::enable_if_t<
                is_discrete_element_v<HeadDElem> && (is_discrete_element_v<TailDElems> && ...),
                int> = 1>
KOKKOS_FUNCTION constexpr auto const& take(HeadDElem const& head, TailDElems const&... tail)
{
    DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
    if constexpr (type_seq_contains_v<detail::TypeSeq<QueryTag>, to_type_seq_t<HeadDElem>>) {
        static_assert(
                (!type_seq_contains_v<detail::TypeSeq<QueryTag>, to_type_seq_t<TailDElems>> && ...),
                "ERROR: tag redundant");
        return head;
    } else {
        static_assert(sizeof...(TailDElems) > 0, "ERROR: tag not found");
        return take<QueryTag>(tail...);
    }
    DDC_IF_NVCC_THEN_POP
}

namespace detail {

/// Returns a reference to the underlying `std::array`
template <class... Tags>
KOKKOS_FUNCTION constexpr std::array<DiscreteElementType, sizeof...(Tags)>& array(
        DiscreteElement<Tags...>& v) noexcept
{
    return v.m_values;
}

/// Returns a reference to the underlying `std::array`
template <class... Tags>
KOKKOS_FUNCTION constexpr std::array<DiscreteElementType, sizeof...(Tags)> const& array(
        DiscreteElement<Tags...> const& v) noexcept
{
    return v.m_values;
}

} // namespace detail

/** A DiscreteElement identifies an element of the discrete dimension
 *
 * Each one is tagged by its associated dimensions.
 */
template <class... Tags>
class DiscreteElement
{
    using tags_seq = detail::TypeSeq<Tags...>;

    friend KOKKOS_FUNCTION constexpr std::array<DiscreteElementType, sizeof...(Tags)>& detail::
            array<Tags...>(DiscreteElement<Tags...>& v) noexcept;

    friend KOKKOS_FUNCTION constexpr std::array<DiscreteElementType, sizeof...(Tags)> const&
    detail::array<Tags...>(DiscreteElement<Tags...> const& v) noexcept;

private:
    std::array<DiscreteElementType, sizeof...(Tags)> m_values;

public:
    using value_type = DiscreteElementType;

    static KOKKOS_FUNCTION constexpr std::size_t size() noexcept
    {
        return sizeof...(Tags);
    }

public:
    KOKKOS_DEFAULTED_FUNCTION constexpr DiscreteElement() = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DiscreteElement(DiscreteElement const&) = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DiscreteElement(DiscreteElement&&) = default;

    template <class... DElems, class = std::enable_if_t<(is_discrete_element_v<DElems> && ...)>>
    explicit KOKKOS_FUNCTION constexpr DiscreteElement(DElems const&... delems) noexcept
        : m_values {take<Tags>(delems...).template uid<Tags>()...}
    {
    }

    template <
            class... Params,
            class = std::enable_if_t<(!is_discrete_element_v<Params> && ...)>,
            class = std::enable_if_t<(std::is_integral_v<Params> && ...)>,
            class = std::enable_if_t<sizeof...(Params) == sizeof...(Tags)>>
    explicit KOKKOS_FUNCTION constexpr DiscreteElement(Params const&... params) noexcept
        : m_values {static_cast<value_type>(params)...}
    {
    }

    KOKKOS_DEFAULTED_FUNCTION ~DiscreteElement() = default;

    KOKKOS_DEFAULTED_FUNCTION DiscreteElement& operator=(DiscreteElement const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION DiscreteElement& operator=(DiscreteElement&& other) = default;

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteElement& operator=(
            DiscreteElement<OTags...> const& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteElement& operator=(DiscreteElement<OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr value_type const& uid_or(value_type const& default_value) const&
    {
        DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
        if constexpr (in_tags_v<QueryTag, tags_seq>) {
            return m_values[type_seq_rank_v<QueryTag, tags_seq>];
        } else {
            return default_value;
        }
        DDC_IF_NVCC_THEN_POP
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr value_type& uid() noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from DiscreteElement");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr value_type const& uid() const noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from DiscreteElement");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <std::size_t N = sizeof...(Tags)>
    KOKKOS_FUNCTION constexpr std::enable_if_t<N == 1, value_type&> uid() noexcept
    {
        return m_values[0];
    }

    template <std::size_t N = sizeof...(Tags)>
    KOKKOS_FUNCTION constexpr std::enable_if_t<N == 1, value_type const&> uid() const noexcept
    {
        return m_values[0];
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    KOKKOS_FUNCTION constexpr DiscreteElement& operator++()
    {
        ++m_values[0];
        return *this;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    KOKKOS_FUNCTION constexpr DiscreteElement operator++(int)
    {
        DiscreteElement const tmp = *this;
        ++m_values[0];
        return tmp;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    KOKKOS_FUNCTION constexpr DiscreteElement& operator--()
    {
        ++m_values[0];
        return *this;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    KOKKOS_FUNCTION constexpr DiscreteElement operator--(int)
    {
        DiscreteElement const tmp = *this;
        ++m_values[0];
        return tmp;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteElement& operator+=(DiscreteVector<OTags...> const& rhs)
    {
        static_assert(((type_seq_contains_v<detail::TypeSeq<OTags>, tags_seq>)&&...));
        ((m_values[type_seq_rank_v<OTags, tags_seq>] += rhs.template get<OTags>()), ...);
        return *this;
    }

    template <
            class IntegralType,
            std::size_t N = sizeof...(Tags),
            class = std::enable_if_t<N == 1>,
            class = std::enable_if_t<std::is_integral_v<IntegralType>>>
    KOKKOS_FUNCTION constexpr DiscreteElement& operator+=(IntegralType const& rhs)
    {
        m_values[0] += rhs;
        return *this;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteElement& operator-=(DiscreteVector<OTags...> const& rhs)
    {
        static_assert(((type_seq_contains_v<detail::TypeSeq<OTags>, tags_seq>)&&...));
        ((m_values[type_seq_rank_v<OTags, tags_seq>] -= rhs.template get<OTags>()), ...);
        return *this;
    }

    template <
            class IntegralType,
            std::size_t N = sizeof...(Tags),
            class = std::enable_if_t<N == 1>,
            class = std::enable_if_t<std::is_integral_v<IntegralType>>>
    KOKKOS_FUNCTION constexpr DiscreteElement& operator-=(IntegralType const& rhs)
    {
        m_values[0] -= rhs;
        return *this;
    }
};

inline std::ostream& operator<<(std::ostream& out, DiscreteElement<> const&)
{
    out << "()";
    return out;
}

template <class Head, class... Tags>
std::ostream& operator<<(std::ostream& out, DiscreteElement<Head, Tags...> const& arr)
{
    out << "(";
    out << uid<Head>(arr);
    ((out << ", " << uid<Tags>(arr)), ...);
    out << ")";
    return out;
}


template <class... Tags, class... OTags>
KOKKOS_FUNCTION constexpr bool operator==(
        DiscreteElement<Tags...> const& lhs,
        DiscreteElement<OTags...> const& rhs) noexcept
{
    return ((lhs.template uid<Tags>() == rhs.template uid<Tags>()) && ...);
}

template <class... Tags, class... OTags>
KOKKOS_FUNCTION constexpr bool operator!=(
        DiscreteElement<Tags...> const& lhs,
        DiscreteElement<OTags...> const& rhs) noexcept
{
    return !(lhs == rhs);
}

template <class Tag>
KOKKOS_FUNCTION constexpr bool operator<(
        DiscreteElement<Tag> const& lhs,
        DiscreteElement<Tag> const& rhs)
{
    return lhs.uid() < rhs.uid();
}

template <class Tag>
KOKKOS_FUNCTION constexpr bool operator<=(
        DiscreteElement<Tag> const& lhs,
        DiscreteElement<Tag> const& rhs)
{
    return lhs.uid() <= rhs.uid();
}

template <class Tag>
KOKKOS_FUNCTION constexpr bool operator>(
        DiscreteElement<Tag> const& lhs,
        DiscreteElement<Tag> const& rhs)
{
    return lhs.uid() > rhs.uid();
}

template <class Tag>
KOKKOS_FUNCTION constexpr bool operator>=(
        DiscreteElement<Tag> const& lhs,
        DiscreteElement<Tag> const& rhs)
{
    return lhs.uid() >= rhs.uid();
}

/// right external binary operators: +, -

template <class... Tags, class... OTags>
KOKKOS_FUNCTION constexpr DiscreteElement<Tags...> operator+(
        DiscreteElement<Tags...> const& lhs,
        DiscreteVector<OTags...> const& rhs)
{
    using detail::TypeSeq;
    static_assert(((type_seq_contains_v<TypeSeq<OTags>, TypeSeq<Tags...>>)&&...));
    DiscreteElement<Tags...> result(lhs);
    result += rhs;
    return result;
}

template <
        class Tag,
        class IntegralType,
        class = std::enable_if_t<std::is_integral_v<IntegralType>>,
        class = std::enable_if_t<!is_discrete_vector_v<IntegralType>>>
KOKKOS_FUNCTION constexpr DiscreteElement<Tag> operator+(
        DiscreteElement<Tag> const& lhs,
        IntegralType const& rhs)
{
    return DiscreteElement<Tag>(uid<Tag>(lhs) + rhs);
}

template <class... Tags, class... OTags>
KOKKOS_FUNCTION constexpr DiscreteElement<Tags...> operator-(
        DiscreteElement<Tags...> const& lhs,
        DiscreteVector<OTags...> const& rhs)
{
    using detail::TypeSeq;
    static_assert(((type_seq_contains_v<TypeSeq<OTags>, TypeSeq<Tags...>>)&&...));
    DiscreteElement<Tags...> result(lhs);
    result -= rhs;
    return result;
}

template <
        class Tag,
        class IntegralType,
        class = std::enable_if_t<std::is_integral_v<IntegralType>>,
        class = std::enable_if_t<!is_discrete_vector_v<IntegralType>>>
KOKKOS_FUNCTION constexpr DiscreteElement<Tag> operator-(
        DiscreteElement<Tag> const& lhs,
        IntegralType const& rhs)
{
    return DiscreteElement<Tag>(uid<Tag>(lhs) - rhs);
}

/// binary operator: -

template <class... Tags, class... OTags>
KOKKOS_FUNCTION constexpr DiscreteVector<Tags...> operator-(
        DiscreteElement<Tags...> const& lhs,
        DiscreteElement<OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    return DiscreteVector<Tags...>((uid<Tags>(lhs) - uid<Tags>(rhs))...);
}

// Gives access to the type of the coordinates of a discrete element
// Example usage : "using Coords = coordinate_of_t<DElem>;"
template <class... DDims>
struct coordinate_of<ddc::DiscreteElement<DDims...>>
{
    // maybe a static_assert on DDims ?
    using type = Coordinate<typename DDims::continuous_dimension_type...>;
};

} // namespace ddc
