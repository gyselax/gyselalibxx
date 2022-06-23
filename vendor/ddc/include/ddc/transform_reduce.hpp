// SPDX-License-Identifier: MIT

#pragma once

#include <utility>

#include <Kokkos_Core.hpp>

#include "ddc/detail/macros.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/for_each.hpp"
#include "ddc/reducer.hpp"

namespace detail {

template <class Reducer>
struct KokkosReducer;

template <class Reducer>
using KokkosReducer_t = typename KokkosReducer<Reducer>::type;

template <class T>
struct KokkosReducer<reducer::sum<T>>
{
    using type = Kokkos::Sum<T>;
};

template <class T>
struct KokkosReducer<reducer::prod<T>>
{
    using type = Kokkos::Prod<T>;
};

template <>
struct KokkosReducer<reducer::land>
{
    using type = Kokkos::LAnd<void>;
};

template <>
struct KokkosReducer<reducer::lor>
{
    using type = Kokkos::LOr<void>;
};

template <class T>
struct KokkosReducer<reducer::band<T>>
{
    using type = Kokkos::BAnd<T>;
};

template <class T>
struct KokkosReducer<reducer::bor<T>>
{
    using type = Kokkos::BOr<T>;
};

template <class T>
struct KokkosReducer<reducer::bxor<T>>
{
    static_assert(std::is_same_v<T, T>, "This reducer is not yet implemented");
};

template <class T>
struct KokkosReducer<reducer::min<T>>
{
    using type = Kokkos::Min<T>;
};

template <class T>
struct KokkosReducer<reducer::max<T>>
{
    using type = Kokkos::Max<T>;
};

template <class T>
struct KokkosReducer<reducer::minmax<T>>
{
    using type = Kokkos::MinMax<T>;
};

/** A serial reduction over a nD domain
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] neutral the neutral element of the reduction operation
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and neutral.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <
        class... DDims,
        class T,
        class BinaryReductionOp,
        class UnaryTransformOp,
        class... DCoords>
inline T transform_reduce_serial(
        DiscreteDomain<DDims...> const& domain,
        T const neutral,
        BinaryReductionOp const& reduce,
        UnaryTransformOp const& transform,
        DCoords const&... dcoords) noexcept
{
    DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
    if constexpr (sizeof...(DCoords) == sizeof...(DDims)) {
        return transform(DiscreteElement<DDims...>(dcoords...));
    } else {
        using CurrentDDim = type_seq_element_t<sizeof...(DCoords), detail::TypeSeq<DDims...>>;
        T result = neutral;
        for (DiscreteElement<CurrentDDim> const ii : select<CurrentDDim>(domain)) {
            result = reduce(
                    result,
                    transform_reduce_serial(domain, neutral, reduce, transform, dcoords..., ii));
        }
        return result;
    }
    DDC_IF_NVCC_THEN_POP
}

template <class Reducer, class Functor, class... DDims>
class TransformReducerKokkosLambdaAdapter
{
    template <class T>
    using index_type = std::size_t;

    Reducer reducer;

    Functor functor;

public:
    DDC_FORCEINLINE_FUNCTION
    TransformReducerKokkosLambdaAdapter(Reducer const& r, Functor const& f) : reducer(r), functor(f)
    {
    }

    DDC_FORCEINLINE_FUNCTION void operator()(
            index_type<DDims>... ids,
            typename Reducer::value_type& a) const
    {
        a = reducer(a, functor(DiscreteElement<DDims...>(ids...)));
    }
};

/** A parallel reduction over a nD domain using the default Kokkos execution space
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] neutral the neutral element of the reduction operation
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and neutral.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class ExecSpace, class DDim0, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce_kokkos(
        DiscreteDomain<DDim0> const& domain,
        T neutral,
        BinaryReductionOp const& reduce,
        UnaryTransformOp const& transform) noexcept
{
    T result = neutral;
    Kokkos::parallel_reduce(
            Kokkos::RangePolicy<ExecSpace>(
                    select<DDim0>(domain).front().uid(),
                    select<DDim0>(domain).back().uid() + 1),
            TransformReducerKokkosLambdaAdapter<
                    BinaryReductionOp,
                    UnaryTransformOp,
                    DDim0>(reduce, transform),
            KokkosReducer_t<BinaryReductionOp>(result));
    return result;
}

/** A parallel reduction over a nD domain using the default Kokkos execution space
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] neutral the neutral element of the reduction operation
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and neutral.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <
        class ExecSpace,
        class DDim0,
        class DDim1,
        class... DDims,
        class T,
        class BinaryReductionOp,
        class UnaryTransformOp>
inline T transform_reduce_kokkos(
        DiscreteDomain<DDim0, DDim1, DDims...> const& domain,
        T neutral,
        BinaryReductionOp const& reduce,
        UnaryTransformOp const& transform) noexcept
{
    T result = neutral;
    Kokkos::Array<std::size_t, 2 + sizeof...(DDims)> const
            begin {select<DDim0>(domain).front().uid(),
                   select<DDim1>(domain).front().uid(),
                   select<DDims>(domain).front().uid()...};
    Kokkos::Array<std::size_t, 2 + sizeof...(DDims)> const
            end {select<DDim0>(domain).back().uid() + 1,
                 select<DDim1>(domain).back().uid() + 1,
                 (select<DDims>(domain).back().uid() + 1)...};
    Kokkos::parallel_reduce(
            Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2 + sizeof...(DDims)>>(begin, end),
            TransformReducerKokkosLambdaAdapter<
                    BinaryReductionOp,
                    UnaryTransformOp,
                    DDim0,
                    DDim1,
                    DDims...>(reduce, transform),
            KokkosReducer_t<BinaryReductionOp>(result));
    return result;
}

} // namespace detail

/** A reduction over a nD domain using the Serial execution policy
 * @param[in] policy the execution policy to use
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] neutral the neutral element of the reduction operation
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and neutral.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class... DDims, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce(
        [[maybe_unused]] serial_host_policy policy,
        DiscreteDomain<DDims...> const& domain,
        T neutral,
        BinaryReductionOp&& reduce,
        UnaryTransformOp&& transform) noexcept
{
    return detail::transform_reduce_serial(
            domain,
            neutral,
            std::forward<BinaryReductionOp>(reduce),
            std::forward<UnaryTransformOp>(transform));
}

/** A reduction over a nD domain using the Kokkos execution policy
 * @param[in] policy the execution policy to use
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] neutral the neutral element of the reduction operation
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and neutral.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class... DDims, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce(
        [[maybe_unused]] parallel_host_policy policy,
        DiscreteDomain<DDims...> const& domain,
        T neutral,
        BinaryReductionOp&& reduce,
        UnaryTransformOp&& transform) noexcept
{
    return detail::transform_reduce_kokkos<Kokkos::DefaultHostExecutionSpace>(
            domain,
            neutral,
            std::forward<BinaryReductionOp>(reduce),
            std::forward<UnaryTransformOp>(transform));
}

/** A reduction over a nD domain using the Kokkos execution policy
 * @param[in] policy the execution policy to use
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] neutral the neutral element of the reduction operation
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and neutral.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class... DDims, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce(
        [[maybe_unused]] parallel_device_policy policy,
        DiscreteDomain<DDims...> const& domain,
        T neutral,
        BinaryReductionOp&& reduce,
        UnaryTransformOp&& transform) noexcept
{
    return detail::transform_reduce_kokkos<Kokkos::DefaultExecutionSpace>(
            domain,
            neutral,
            std::forward<BinaryReductionOp>(reduce),
            std::forward<UnaryTransformOp>(transform));
}

/** A reduction over a nD domain using the default execution policy
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] neutral the neutral element of the reduction operation
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and neutral.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class... DDims, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce(
        DiscreteDomain<DDims...> const& domain,
        T neutral,
        BinaryReductionOp&& reduce,
        UnaryTransformOp&& transform) noexcept
{
    return transform_reduce(
            default_policy(),
            domain,
            neutral,
            std::forward<BinaryReductionOp>(reduce),
            std::forward<UnaryTransformOp>(transform));
}
