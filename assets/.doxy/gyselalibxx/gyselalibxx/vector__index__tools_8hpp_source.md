

# File vector\_index\_tools.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**vector\_index\_tools.hpp**](vector__index__tools_8hpp.md)

[Go to the documentation of this file](vector__index__tools_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

template <class... Dims>
using VectorIndexSet = ddc::detail::TypeSeq<Dims...>;

namespace tensor_tools {

template <class AnyVectorIndexSet>
struct GetContravariantDims;

template <class... Dims>
struct GetContravariantDims<VectorIndexSet<Dims...>>
{
    using type = VectorIndexSet<
            std::conditional_t<Dims::IS_CONTRAVARIANT, Dims, typename Dims::Dual>...>;
};

template <class AnyVectorIndexSet>
struct GetCovariantDims;

template <class... Dims>
struct GetCovariantDims<VectorIndexSet<Dims...>>
{
    using type
            = VectorIndexSet<std::conditional_t<Dims::IS_COVARIANT, Dims, typename Dims::Dual>...>;
};

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

template <class... Dims>
struct is_covariant_vector_index_set<VectorIndexSet<Dims...>>
{
    static constexpr bool value = (Dims::IS_COVARIANT && ...);
};

template <class... Dims>
struct is_contravariant_vector_index_set<VectorIndexSet<Dims...>>
{
    static constexpr bool value = (Dims::IS_CONTRAVARIANT && ...);
};

template <class VectorIndexSet>
struct vector_index_set_dual;

template <class... Dims>
struct vector_index_set_dual<VectorIndexSet<Dims...>>
{
    using type = VectorIndexSet<typename Dims::Dual...>;
};

template <char Id, class AssociatedVectorIndexSet>
struct VectorIndexIdMap
{
    static_assert(
            tensor_tools::is_vector_index_set<AssociatedVectorIndexSet>::value,
            "The possible index values must be described by a VectorIndexSet.");
    static constexpr char id = Id;
    using possible_idx_values = AssociatedVectorIndexSet;
};

} // namespace tensor_tools

template <class Type>
static constexpr bool is_vector_index_set_v = tensor_tools::is_vector_index_set<Type>::value;

template <class AnyVectorIndexSet>
using get_contravariant_dims_t =
        typename tensor_tools::GetContravariantDims<AnyVectorIndexSet>::type;

template <class AnyVectorIndexSet>
using get_covariant_dims_t = typename tensor_tools::GetCovariantDims<AnyVectorIndexSet>::type;

template <class VectorIndexSet>
using vector_index_set_dual_t = typename tensor_tools::vector_index_set_dual<VectorIndexSet>::type;

template <class VectorIndexSet>
static constexpr bool is_covariant_vector_index_set_v
        = tensor_tools::is_covariant_vector_index_set<VectorIndexSet>::value;

template <class VectorIndexSet>
static constexpr bool is_contravariant_vector_index_set_v
        = tensor_tools::is_contravariant_vector_index_set<VectorIndexSet>::value;
```


