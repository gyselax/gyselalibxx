// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

/**
 * @brief A type alias to describe a set of dimensions that can be used
 * to index a vector (e.g. x to get E_x).
 */
template <class... Dims>
using VectorIndexSet = ddc::detail::TypeSeq<Dims...>;

/**
 * @brief A helper structure to recognise a VectorIndexSet type.
 */
template <class Type>
struct is_vector_index_set : std::false_type
{
};

template <class... Dims>
struct is_vector_index_set<VectorIndexSet<Dims...>> : std::true_type
{
};

template <class VectorIndexSet>
struct is_covariant_vector_index_set;
template <class VectorIndexSet>
struct is_contravariant_vector_index_set;

/**
 * @brief A helper structure to check if all the dimensions in a
 * VectorIndexSet can represent covariant indices.
 */
template <class... Dims>
struct is_covariant_vector_index_set<VectorIndexSet<Dims...>>
{
    /// Compile-time boolean indicating if the vector indices are covariant.
    static constexpr bool value
            = std::conjunction_v<std::integral_constant<bool, Dims::IS_COVARIANT>...>;
};

/**
 * @brief A helper structure to check if all the dimensions in a
 * VectorIndexSet can represent contravariant indices.
 */
template <class... Dims>
struct is_contravariant_vector_index_set<VectorIndexSet<Dims...>>
{
    /// Compile-time boolean indicating if the vector indices are contravariant.
    static constexpr bool value
            = std::conjunction_v<std::integral_constant<bool, Dims::IS_CONTRAVARIANT>...>;
};

/**
 * @brief A helper structure to find a VectorIndexSet describing
 * the covariant indices from a VectorIndexSet describing contravariant
 * indices or to find a VectorIndexSet describing the contravariant
 * indices from a VectorIndexSet describing covariant indices.
 */
template <class VectorIndexSet>
struct vector_index_set_dual;

/**
 * @brief The implementation of vector_index_set_dual for a VectorIndexSet.
 */
template <class... Dims>
struct vector_index_set_dual<VectorIndexSet<Dims...>>
{
    /// A type alias describing the matching VectorIndexSet.
    using type = VectorIndexSet<typename Dims::Dual...>;
};

/**
 * @brief A type alias to find a VectorIndexSet describing
 * the covariant indices from a VectorIndexSet describing contravariant
 * indices or to find a VectorIndexSet describing the contravariant
 * indices from a VectorIndexSet describing covariant indices.
 */
template <class VectorIndexSet>
using vector_index_set_dual_t = typename vector_index_set_dual<VectorIndexSet>::type;

/**
 * @brief A compile-time boolean to check if all the dimensions in a
 * VectorIndexSet can represent covariant indices.
 */
template <class VectorIndexSet>
static constexpr bool is_covariant_vector_index_set_v
        = is_covariant_vector_index_set<VectorIndexSet>::value;

/**
 * @brief A compile-time boolean to check if all the dimensions in a
 * VectorIndexSet can represent contravariant indices.
 */
template <class VectorIndexSet>
static constexpr bool is_contravariant_vector_index_set_v
        = is_contravariant_vector_index_set<VectorIndexSet>::value;
