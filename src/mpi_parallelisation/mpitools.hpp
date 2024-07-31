// SPDX-License-Identifier: MIT
#pragma once
#include <mpi.h>

#include "ddc_aliases.hpp"

/**
 * @brief An internal tag used to dsecribe an artificial dimension describing the MPI rank
 * where the scattered information will be sent to or where the gathered information will
 * be collected from.
 */
template <class DistributedDim>
struct MPIDim
{
    /// The dimension which is distributed over the MPI ranks.
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

/// A helper to get the MPI type descriptor from an element type
template <class ElementType>
static const MPI_Datatype MPI_type_descriptor_t
        = detail::MPITypeDescriptor<ElementType>::get_type();

namespace detail {


/**
 * @brief A class to determine if a set of types found within a TypeSeq are adjacent to one another.
 * @tparam Query The TypeSeq describing the possibly adjacent subset.
 * @tparam ContainerTypeSeq The TypeSeq which contains the subset.
 */
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

/**
 * @brief A tool to insert a tag into an existing TypeSeq at a specified position.
 *
 * @tparam TagToInsert The tag to be inserted.
 * @tparam POS The position at which the tag should be inserted.
 * @tparam TypeSeq The TypeSeq into which the tag should be inserted.
 */
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

/**
 * @brief A helper constant to determine if a set of types found within a TypeSeq are adjacent to
 * one another.
 * @tparam Query The TypeSeq describing the possibly adjacent subset.
 * @tparam ContainerTypeSeq The TypeSeq which contains the subset.
 */
template <class Query, class ContainerTypeSeq>
constexpr bool dimensions_are_adjacent_v
        = detail::DimensionsAreAdjacent<Query, ContainerTypeSeq>::value;

/**
 * @brief A tool to insert a tag into an existing TypeSeq at a specified position.
 *
 * @tparam TagToInsert The tag to be inserted.
 * @tparam POS The position at which the tag should be inserted.
 * @tparam TypeSeq The TypeSeq into which the tag should be inserted.
 */
template <class TagToInsert, std::size_t POS, class TypeSeq>
using insert_into_type_seq_t = typename detail::InsertIntoTypeSeq<TagToInsert, POS, TypeSeq>::type;

/**
 * @brief A tool to insert a tag into an existing TypeSeq immediately preceeding an existing
 * subset of tags.
 *
 * @tparam TagToInsert The tag to be inserted.
 * @tparam SubSeq The subset of tags found in the TypeSeq before which the tag should be inserted.
 * @tparam TypeSeq The TypeSeq into which the tag should be inserted.
 */
template <class TagToInsert, class SubSeq, class TypeSeq>
using insert_into_seq_before_t = typename detail::InsertIntoTypeSeq<
        TagToInsert,
        ddc::type_seq_rank_v<ddc::type_seq_element_t<0, SubSeq>, TypeSeq>,
        TypeSeq>::type;

namespace detail {

/**
 * @brief Insert MPI distribution tags into an existing TypeSeq.
 *
 * The MPI tags are each associated with an index range. The MPI tags are inserted into the TypeSeq
 * immediately preceeding the tag with which they are associated. This allows an index range to be
 * split along the axes on which it will be scattered.
 *
 * E.g. if MPI<Phi> is inserted into <R, Theta, Phi, VPar> we would obtain:
 * <R, Theta, MPI<Phi>, Phi, VPar>
 *
 * @tparam MPISeq A type sequence containing the MPI tags which should be inserted
 * @tparam TypeSeq The type sequence which should be inserted.
 */
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

/**
 * @brief Insert MPI distribution tags into an existing TypeSeq.
 *
 * The MPI tags are each associated with an index range. The MPI tags are inserted into the TypeSeq
 * immediately preceeding the tag with which they are associated. This allows an index range to be
 * split along the axes on which it will be scattered.
 *
 * E.g. if MPI<Phi> is inserted into <R, Theta, Phi, VPar> we would obtain:
 * <R, Theta, MPI<Phi>, Phi, VPar>
 *
 * @tparam MPISeq A type sequence containing the MPI tags which should be inserted
 * @tparam TypeSeq The type sequence which should be inserted.
 */
template <class MPISeq, class TypeSeq>
using insert_mpi_tags_into_seq_t = typename detail::InsertMPITags<MPISeq, TypeSeq>::type;
