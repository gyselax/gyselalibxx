#pragma once

#include <limits>
#include <utility>

namespace ddc::reducer {

template <class T>
struct sum
{
    using value_type = T;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs + rhs;
    }
};

template <class T>
struct prod
{
    using value_type = T;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs * rhs;
    }
};

struct land
{
    using value_type = bool;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs && rhs;
    }
};

struct lor
{
    using value_type = bool;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs || rhs;
    }
};

template <class T>
struct band
{
    using value_type = T;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs & rhs;
    }
};

template <class T>
struct bor
{
    using value_type = T;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs | rhs;
    }
};

template <class T>
struct bxor
{
    using value_type = T;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs ^ rhs;
    }
};

template <class T>
struct min
{
    using value_type = T;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs < rhs ? lhs : rhs;
    }
};

template <class T>
struct max
{
    using value_type = T;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return lhs > rhs ? lhs : rhs;
    }
};

template <class T>
struct minmax
{
    using value_type = std::pair<T, T>;

    constexpr value_type operator()(value_type const& lhs, value_type const& rhs) const noexcept
    {
        return value_type(
                lhs.first < rhs.first ? lhs.first : rhs.first,
                lhs.second > rhs.second ? lhs.second : rhs.second);
    }
};

} // namespace ddc::reducer
