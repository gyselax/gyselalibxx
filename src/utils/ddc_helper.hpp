// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


namespace ddcHelper {
//TODO: this should be directly handled by ddc::Discretisation really,
//      in the meantime, we do it ourselves
/**
 * Calculate the total length of a non-periodic domain.
 *
 * @param idx_range The domain on which the length should be calculated.
 *
 * @return The length of the domain.
 */
template <class IDim>
constexpr std::enable_if_t<!IDim::continuous_dimension_type::PERIODIC, double>
total_interval_length(IdxRange<IDim> const& idx_range)
{
    return std::fabs(ddc::rlength(idx_range));
}

//TODO: this should be directly handled by ddc::Discretisation really,
//      in the meantime, we do it ourselves
/**
 * Calculate the total length of a uniform periodic domain.
 *
 * @param idx_range The domain on which the length should be calculated.
 *
 * @return The length of the domain.
 */
template <class IDim>
constexpr std::enable_if_t<
        IDim::continuous_dimension_type::PERIODIC && ddc::is_uniform_point_sampling_v<IDim>,
        double>
total_interval_length(IdxRange<IDim> const& idx_range)
{
    return std::fabs(ddc::rlength(idx_range) + ddc::step<IDim>());
}

//TODO: this should be directly handled by ddc::Discretisation really,
//      in the meantime, we do it ourselves
/**
 * Calculate the total length of a non-uniform periodic domain.
 *
 * @param idx_range The domain on which the length should be calculated.
 *
 * @return The length of the domain.
 */
template <class IDim>
constexpr std::enable_if_t<
        IDim::continuous_dimension_type::PERIODIC && ddc::is_non_uniform_point_sampling_v<IDim>,
        double>
total_interval_length(IdxRange<IDim> const& idx_range)
{
    IdxRange<IDim> dom_periodic(idx_range.front(), idx_range.extents() + 1);
    return std::fabs(ddc::rlength(dom_periodic));
}

//TODO: this should be directly handled by ddc::Discretisation really,
//      in the meantime, we do it ourselves
/**
 * @brief Calculate the Coordinate inside the domain.
 *
 * In the case of a periodic domain, a Coordinate can sometimes be found
 * outside the domain. In this case it is useful to be able to find the
 * equivalent coordinate inside the domain. This function makes that
 * possible.
 *
 * @param[in] coord
 *      The 1D coordinate we want to compute inside the domain.
 * @param[in] idx_range
 *      The domain where the coordinate is defined.
 *
 * @return The equivalent coordinate inside the domain.
 */
template <class IDim>
constexpr std::enable_if_t<
        IDim::continuous_dimension_type::PERIODIC,
        Coord<typename IDim::continuous_dimension_type>>
restrict_to_idx_range(
        Coord<typename IDim::continuous_dimension_type> coord,
        IdxRange<IDim> const& idx_range)
{
    using Coord1D = Coord<typename IDim::continuous_dimension_type>;
    double const x_min = ddc::rmin(idx_range);
    double const length = total_interval_length(idx_range);
    double const x_max = x_min + length;

    assert(length > 0);
    coord -= x_min;
    if (fabs(coord) > 10 * length) {
        double periodic_factor = 2 * M_PI / length;
        double coord_2pi = double(coord) * periodic_factor;
        coord_2pi = std::copysign(std::acos(std::cos(coord_2pi)), std::sin(coord_2pi));
        coord = coord_2pi < 0 ? Coord1D((coord_2pi + 2 * M_PI) / periodic_factor)
                              : Coord1D((coord_2pi) / periodic_factor);
    }
    coord += x_min;
    while (coord < x_min)
        coord += length;
    while (coord >= x_max)
        coord -= length;
    return coord;
}

/**
 * @brief Calculate the Coordinate inside the domain.
 *
 * In the case of a periodic domain, a Coordinate can sometimes be found
 * outside the domain. In this case it is useful to be able to find the
 * equivalent coordinate inside the domain. This function makes that
 * possible.
 *
 * @param[inout] coord
 *      The 1D coordinate we want to compute inside the domain.
 *
 * @return The equivalent coordinate inside the domain.
 */
template <class BSpline>
KOKKOS_INLINE_FUNCTION void restrict_to_bspline_domain(
        Coord<typename BSpline::continuous_dimension_type>& coord)
{
    if constexpr (BSpline::is_periodic()) {
        if (coord < ddc::discrete_space<BSpline>().rmin()
            || coord > ddc::discrete_space<BSpline>().rmax()) {
            coord -= Kokkos::floor(
                             (coord - ddc::discrete_space<BSpline>().rmin())
                             / ddc::discrete_space<BSpline>().length())
                     * ddc::discrete_space<BSpline>().length();
        }
    }
}

/**
 * @brief Dump the coordinates found on a domain into a span.
 *
 * @param[out] dump_coord The span which will contain the coordinates.
 * @param[in] sampler The domain indicating the coordinates.
 */
template <class ExecSpace, class Grid1D, class Layout, class MemorySpace>
inline void dump_coordinates(
        ExecSpace exec_space,
        DField<IdxRange<Grid1D>, Layout, MemorySpace> dump_coord)
{
    ddc::parallel_for_each(
            exec_space,
            dump_coord.domain(),
            KOKKOS_LAMBDA(Idx<Grid1D> i) { dump_coord(i) = ddc::coordinate(i); });
}

/**
 * @brief Dump the coordinates found on a domain into a span.
 *
 * @param[out] dump_coord The span which will contain the coordinates.
 * @param[in] sampler The domain indicating the coordinates.
 */
template <class ExecSpace, class Grid1D, class Layout, class MemorySpace>
inline void dump_coordinates(
        ExecSpace exec_space,
        Field<Coord<typename Grid1D::continuous_dimension_type>,
              IdxRange<Grid1D>,
              Layout,
              MemorySpace> dump_coord)
{
    ddc::parallel_for_each(
            exec_space,
            dump_coord.domain(),
            KOKKOS_LAMBDA(Idx<Grid1D> i) { dump_coord(i) = ddc::coordinate(i); });
}

/**
 * @brief Calculate the maximum distance between two adjacent points 
 * within an IdxRange
 *
 * @param idx_range The domain on which the distance should be calculated.
 *
 * @return The maximum distance between two adjacent points.
 */
template <class GridDim>
double maximum_distance_between_two_points(IdxRange<GridDim> const& idx_range)
{
    using IdxStep = IdxStep<GridDim>;
    using IdxDim = Idx<GridDim>;

    IdxDim const ifront = idx_range.front();
    IdxStep const step(1);

    IdxRange<GridDim> idx_range_chopped = idx_range.remove_first(step);
    double const max_dist = ddc::parallel_transform_reduce(
            Kokkos::DefaultHostExecutionSpace(),
            idx_range_chopped,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxDim const ix) {
                return ddc::coordinate(ix) - ddc::coordinate(ix - step);
            });


    return max_dist;
}


} // namespace ddcHelper

//-----------------------------------------------------------------------------

namespace detail {
template <class, class>
struct OnMemorySpace
{
};

/**
 * @brief Set a `ddc::Chunk` on a given NewMemorySpace.
 * @tparam NewMemorySpace The new memory space. 
 * @tparam ElementType Type of the elements in the ddc::Chunk.
 * @tparam SupportType Type of the domain of the ddc::Chunk.
 * @tparam Allocator Allocator type (see ddc::KokkosAllocator).
 * @see ddc::Chunk
 */
template <class NewMemorySpace, class ElementType, class SupportType, class Allocator>
struct OnMemorySpace<NewMemorySpace, ddc::Chunk<ElementType, SupportType, Allocator>>
{
    using type = typename ddc::
            Chunk<ElementType, SupportType, ddc::KokkosAllocator<ElementType, NewMemorySpace>>;
};

/**
 * @brief Get a new `ddc::ChunkSpan` type with the same parametrisation
 * except in the memory space which is set to NewMemorySpace.
 * @tparam NewMemorySpace The new memory space. 
 * @tparam ElementType Type of the elements in the ddc::Chunk.
 * @tparam SupportType Type of the domain of the ddc::Chunk.
 * @tparam Layout Layout tag (see Kokkos).
 * @tparam MemorySpace The original memory space of the chunk.
 * @see ddc::ChunkSpan
 */
template <
        class NewMemorySpace,
        class ElementType,
        class SupportType,
        class Layout,
        class MemorySpace>
struct OnMemorySpace<NewMemorySpace, ddc::ChunkSpan<ElementType, SupportType, Layout, MemorySpace>>
{
    using type = typename ddc::ChunkSpan<ElementType, SupportType, Layout, NewMemorySpace>;
};

template <template <class Tag> class Templ, class TypeSeq>
struct ApplyTemplateToTypeSeq;

template <template <class Tag> class Templ, class... Tags>
struct ApplyTemplateToTypeSeq<Templ, ddc::detail::TypeSeq<Tags...>>
{
    using type = ddc::detail::TypeSeq<Templ<Tags>...>;
};

/// R contains all elements that are in A and B.
/// Remark 1: This operation preserves the order from B.
/// Remark 2: It is similar to the set intersection in the set theory.
/// Example: A = [a, b, c], B = [z, c, y], R = [c]
template <class TagSeqA, class TagSeqB, class TagSeqR>
struct TypeSeqIntersection;

template <class TypeSeqA, class TypeSeqR>
struct TypeSeqIntersection<TypeSeqA, ddc::detail::TypeSeq<>, TypeSeqR>
{
    using type = TypeSeqR;
};

template <class TypeSeqA, class HeadTagsB, class... TailTagsB, class... TagsR>
struct TypeSeqIntersection<
        TypeSeqA,
        ddc::detail::TypeSeq<HeadTagsB, TailTagsB...>,
        ddc::detail::TypeSeq<TagsR...>>
{
    using type = std::conditional_t<
            ddc::in_tags_v<HeadTagsB, TypeSeqA>,
            typename TypeSeqIntersection<
                    TypeSeqA,
                    ddc::detail::TypeSeq<TailTagsB...>,
                    ddc::detail::TypeSeq<TagsR..., HeadTagsB>>::type,
            typename TypeSeqIntersection<
                    TypeSeqA,
                    ddc::detail::TypeSeq<TailTagsB...>,
                    ddc::detail::TypeSeq<TagsR...>>::type>;
};


/// \endcond

} // namespace detail

/**
 * @brief Alias template helper returning the type of a `ddc::Chunk`, a `ddc::ChunkSpan`, a `VectorFieldMem`
 * or a `VectorField` on a MemorySpace.
 */
template <class MemorySpace, class C>
using on_memory_space_t = typename detail::OnMemorySpace<MemorySpace, C>::type;

/**
 * @brief Alias template helper returning the "host" version of a `ddc::Chunk`, a `ddc::ChunkSpan`,
 * a `VectorFieldMem` or a `VectorField`.
 */
template <class C>
using host_t = on_memory_space_t<Kokkos::HostSpace, C>;

/**
 * @brief Alias template helper returning the "device" version of a `ddc::Chunk`, a `ddc::ChunkSpan`,
 * a `VectorFieldMem` or a `VectorField`.
 */
template <class C>
using device_t = on_memory_space_t<Kokkos::DefaultExecutionSpace::memory_space, C>;



namespace ddcHelper {
/// A helper to get a type sequence by applying a template to a sequence of type tags.
template <template <class Tag> class Templ, class TypeSeq>
using apply_template_to_type_seq_t = typename detail::ApplyTemplateToTypeSeq<Templ, TypeSeq>::type;

/// A helper to find all types which are found in both TypeSeq1 and TypeSeq2
template <class TypeSeq1, class TypeSeq2>
using type_seq_intersection_t =
        typename detail::TypeSeqIntersection<TypeSeq1, TypeSeq2, ddc::detail::TypeSeq<>>::type;
} // namespace ddcHelper
