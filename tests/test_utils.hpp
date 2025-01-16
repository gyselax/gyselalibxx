// SPDX-License-Identifier: MIT
#pragma once

#include <tuple>
#include <type_traits>
#include <utility>

#include <gtest/gtest.h>

namespace detail {
template <class S>
struct to_tuple
{
    using type = S;
};

template <class T, T... Ints>
struct to_tuple<std::integer_sequence<T, Ints...>>
{
    using type = std::tuple<std::integral_constant<T, Ints>...>;
};

template <class T, class U>
struct to_tuple<std::pair<T, U>>
{
    using type = std::tuple<T, U>;
};
} // namespace detail

/// Transform a sequence S to a tuple:
/// - a std::integer_sequence<T, Ints...> to a std::tuple<std::integral_constant<T, Ints>...>
/// - a std::pair<T, U> to a std::tuple<T, U>
/// - identity otherwise (std::tuple)
template <class S>
using to_tuple_t = typename detail::to_tuple<S>::type;

namespace detail {
template <class TupleOfTuples, class Tuple>
struct for_each_tuple_cat;

template <class... Tuples, class Tuple>
struct for_each_tuple_cat<std::tuple<Tuples...>, Tuple>
{
    using type = std::tuple<
            decltype(std::tuple_cat(std::declval<Tuples>(), std::declval<Tuple>()))...>;
};

} // namespace detail

/// Construct a tuple of tuples that is the result of the concatenation of the tuples in TupleOfTuples with Tuple.
template <class TupleOfTuples, class Tuple>
using for_each_tuple_cat_t = typename detail::for_each_tuple_cat<TupleOfTuples, Tuple>::type;

static_assert(std::is_same_v<
              for_each_tuple_cat_t<
                      std::tuple<std::tuple<double, double>, std::tuple<int, double>>,
                      std::tuple<int>>,
              std::tuple<std::tuple<double, double, int>, std::tuple<int, double, int>>>);

static_assert(std::is_same_v<
              for_each_tuple_cat_t<std::tuple<std::tuple<double, double>>, std::tuple<int>>,
              std::tuple<std::tuple<double, double, int>>>);

namespace detail {
template <class InTupleOfTuples, class OutTupleOfTuples>
struct cartesian_product_impl;

template <class... HeadArgs, class... TailTuples, class OutTupleOfTuples>
struct cartesian_product_impl<std::tuple<std::tuple<HeadArgs...>, TailTuples...>, OutTupleOfTuples>
    : cartesian_product_impl<
              std::tuple<TailTuples...>,
              decltype(std::tuple_cat(
                      std::declval<
                              for_each_tuple_cat_t<OutTupleOfTuples, std::tuple<HeadArgs>>>()...))>
{
};

template <class OutTupleOfTuples>
struct cartesian_product_impl<std::tuple<>, OutTupleOfTuples>
{
    using type = OutTupleOfTuples;
};
} // namespace detail

/// Generate a std::tuple cartesian product from multiple tuple-like structures (std::tuple, std::integer_sequence and std::pair)
/// Do not rely on the ordering result.
template <class... InTuplesLike>
using cartesian_product_t = typename detail::cartesian_product_impl<
        std::tuple<to_tuple_t<InTuplesLike>...>,
        std::tuple<std::tuple<>>>::type;

static_assert(std::is_same_v<
              cartesian_product_t<std::tuple<int, float>, std::tuple<double>>,
              std::tuple<std::tuple<int, double>, std::tuple<float, double>>>);

namespace detail {
template <class T>
struct tuple_to_types
{
    using type = T;
};

template <class... Args>
struct tuple_to_types<std::tuple<Args...>>
{
    using type = testing::Types<Args...>;
};
} // namespace detail

/// Transform a std::tuple<Args...> to a testing::Types<Args...>, identity otherwise
template <class T>
using tuple_to_types_t = typename detail::tuple_to_types<T>::type;
