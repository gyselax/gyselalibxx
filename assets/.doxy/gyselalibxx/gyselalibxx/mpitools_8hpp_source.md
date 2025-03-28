

# File mpitools.hpp

[**File List**](files.md) **>** [**mpi\_parallelisation**](dir_a35b8fd75f8fad0c2619b083ab571e51.md) **>** [**mpitools.hpp**](mpitools_8hpp.md)

[Go to the documentation of this file](mpitools_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <mpi.h>

#include "ddc_aliases.hpp"

template <class DistributedDim>
struct MPIDim
{
    using distributed_dim = DistributedDim;
};

namespace detail {

template <class ElementType>
struct MPITypeDescriptor
{
    static_assert(!std::is_same_v<ElementType, ElementType>, "MPI Datatype not recognised.");
};

template <>
struct MPITypeDescriptor<double>
{
    static MPI_Datatype get_type() noexcept
    {
        return MPI_DOUBLE;
    }
};

template <>
struct MPITypeDescriptor<unsigned long>
{
    static MPI_Datatype get_type() noexcept
    {
        return MPI_UNSIGNED_LONG;
    }
};
} // namespace detail

template <class ElementType>
static const MPI_Datatype MPI_type_descriptor_t
        = detail::MPITypeDescriptor<ElementType>::get_type();

namespace detail {


template <class QueryTypeSeq, class ContainerTypeSeq>
struct DimensionsAreAdjacent;

template <class... Query, class ContainerTypeSeq>
struct DimensionsAreAdjacent<ddc::detail::TypeSeq<Query...>, ContainerTypeSeq>
{
    static constexpr bool check_adjacent()
    {
        static_assert((ddc::in_tags_v<Query, ContainerTypeSeq> && ...));
        std::array<std::size_t, sizeof...(Query)> ranks(
                {ddc::type_seq_rank_v<Query, ContainerTypeSeq>...});
        bool ordered = true;
        for (int i(1); i < sizeof...(Query); ++i) {
            ordered &= ((ranks[i] - ranks[i - 1]) == 1);
        }
        return ordered;
    }
    static constexpr bool value = check_adjacent();
};

template <class TagToInsert, std::size_t POS, class TypeSeq>
struct InsertIntoTypeSeq
{
    using first_type_seq_elem_seq = ddc::detail::TypeSeq<ddc::type_seq_element_t<0, TypeSeq>>;
    using type = ddc::type_seq_merge_t<
            first_type_seq_elem_seq,
            typename InsertIntoTypeSeq<
                    TagToInsert,
                    POS - 1,
                    ddc::type_seq_remove_t<TypeSeq, first_type_seq_elem_seq>>::type>;
};

template <class TagToInsert, class... Tags>
struct InsertIntoTypeSeq<TagToInsert, 0, ddc::detail::TypeSeq<Tags...>>
{
    using type = ddc::detail::TypeSeq<TagToInsert, Tags...>;
};

} // namespace detail

template <class Query, class ContainerTypeSeq>
constexpr bool dimensions_are_adjacent_v
        = detail::DimensionsAreAdjacent<Query, ContainerTypeSeq>::value;

template <class TagToInsert, std::size_t POS, class TypeSeq>
using insert_into_type_seq_t = typename detail::InsertIntoTypeSeq<TagToInsert, POS, TypeSeq>::type;

template <class TagToInsert, class SubSeq, class TypeSeq>
using insert_into_seq_before_t = typename detail::InsertIntoTypeSeq<
        TagToInsert,
        ddc::type_seq_rank_v<ddc::type_seq_element_t<0, SubSeq>, TypeSeq>,
        TypeSeq>::type;

namespace detail {

template <class MPISeq, class TypeSeq>
struct InsertMPITags;

template <class HeadDistribDim, class... TailDistribDims, class TypeSeq>
struct InsertMPITags<
        ddc::detail::TypeSeq<MPIDim<HeadDistribDim>, MPIDim<TailDistribDims>...>,
        TypeSeq>
{
    using type = typename InsertMPITags<
            ddc::detail::TypeSeq<MPIDim<TailDistribDims>...>,
            insert_into_type_seq_t<
                    MPIDim<HeadDistribDim>,
                    ddc::type_seq_rank_v<HeadDistribDim, TypeSeq>,
                    TypeSeq>>::type;
};

template <class TypeSeq>
struct InsertMPITags<ddc::detail::TypeSeq<>, TypeSeq>
{
    using type = TypeSeq;
};
} // namespace detail

template <class MPISeq, class TypeSeq>
using insert_mpi_tags_into_seq_t = typename detail::InsertMPITags<MPISeq, TypeSeq>::type;
```


