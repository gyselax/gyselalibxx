// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_aliases.hpp"
#include "indexed_tensor.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "l_norm_tools.hpp"
#include "multipatch_type.hpp"
#include "tensor.hpp"
#include "types.hpp"
#include "vector_index_tools.hpp"


/**
 * @brief Patch locator specialised for the Restelli mesh.
 *
 * This operator locates on which patch a given physical coordinate is.
 * 
 * @tparam MultipatchMapping A MultipatchType type containing the coordinate 
 *          transformations on each patch.
 */
template <class MultipatchMapping>
class PatchLocatorNewton;


/**
 * @brief Patch locator specialised for the Restelli mesh.
 *
 * @tparam MultipatchMapping A MultipatchType type containing the coordinate 
 *          transformations on each patch.
 */
template <class... Patches, class... Mappings>
class PatchLocatorNewton<
        detail::MultipatchType<ddc::TypeSeq<Patches...>, ddc::TypeSeq<Mappings...>>>
{
public:
    /// @brief Sequence ddc::TypeSeq of patch tags.
    using patch_seq = ddc::TypeSeq<Patches...>;

    /// @brief Default value to define outside domain (not a patch).
    static constexpr int outside_domain = -1;

private:
    using mapping_on_first_patch = ddc::type_seq_element_t<0, ddc::TypeSeq<Mappings...>>;

    using ArgBasis = ddc::to_type_seq_t<typename mapping_on_first_patch::CoordArg>;
    using ResBasis = ddc::to_type_seq_t<typename mapping_on_first_patch::CoordResult>;

    static_assert(
            (std::is_same_v<ArgBasis, ddc::to_type_seq_t<typename Mappings::CoordArg>> && ...),
            "The argument continuous dimensions of all the coordinate transformations should be "
            "the same type on all the patches.");

    static_assert(
            (std::is_same_v<ResBasis, ddc::to_type_seq_t<typename Mappings::CoordResult>> && ...),
            "The result continuous dimensions of all the coordinate transformations should be "
            "the same type on all the patches.");


    using DimArg1 = ddc::type_seq_element_t<0, ArgBasis>;
    using DimArg2 = ddc::type_seq_element_t<1, ArgBasis>;

    using DimRes1 = ddc::type_seq_element_t<0, ResBasis>;
    using DimRes2 = ddc::type_seq_element_t<1, ResBasis>;

    using MultipatchMapping
            = detail::MultipatchType<ddc::TypeSeq<Patches...>, ddc::TypeSeq<Mappings...>>;

private:
    static constexpr std::size_t n_patches = ddc::type_seq_size_v<patch_seq>;

    MultipatchMapping const m_coord_transformations;

    const std::size_t m_counter_max;
    const double m_TOL;

public:
    /**
     * @brief Instantiate the patch locator solving with a Newton algorithm.
     * 
     * @param coord_transformations A MultipatchType of coordinate transformations from the logical 
     *          domain to the physical domain of each patch
     * @param counter_max Maximum number of loops for the Newton algorithm.
     * @param TOL Tolerance on the error for the Newton algorithm.
    */
    PatchLocatorNewton(
            MultipatchMapping const& coord_transformations,
            std::size_t const counter_max = 10,
            double const TOL = 1e-12)
        : m_coord_transformations(coord_transformations)
        , m_counter_max(counter_max)
        , m_TOL(TOL)
    {
    }

    /**
     * @brief Get the patch where the given physical coordinate is.
     * 
     * We apply a Newton algorithm to determine the patch.
     *
     * @param[in] coord The given physical coordinate.
     * @return [int] The patch index where the physical coordinate. If the coordinate
     *              is outside of the domain, it returns a negative value.
     */
    KOKKOS_INLINE_FUNCTION int operator()(Coord<DimRes1, DimRes2> const coord) const
    {
        Coord<DimArg1, DimArg2> coord_log;
        return find_patch(coord_log, coord);
    }

    /**
     * @brief Get the patch where the given physical coordinate is.
     * 
     * We apply a Newton algorithm to determine the patch.
     *
     * @param[in, out] coord_log The logical coordinate corresponding to the given 
     *              physical coordinate. The given value is used as initial condition
     *              in the Newton algorithm.
     * @param[in] coord The given physical coordinate.
     * @return [int] The patch index where the physical coordinate. If the coordinate
     *              is outside of the domain, it returns a negative value.
     */
    KOKKOS_INLINE_FUNCTION int operator()(
            Coord<DimArg1, DimArg2>& coord_log,
            Coord<DimRes1, DimRes2> const coord) const
    {
        return find_patch(coord_log, coord);
    }

private:
    /**
     * @brief Find recursively on which patch the physical coordinate 
     * is located and update the logical coordinate corresponding to
     * the physical coordinate.
     */
    template <std::size_t I = 0>
    KOKKOS_INLINE_FUNCTION int find_patch(
            Coord<DimArg1, DimArg2>& coord_log,
            Coord<DimRes1, DimRes2> const& coord) const
    {
        if constexpr (I == n_patches) {
            return outside_domain;
        } else {
            using P = ddc::type_seq_element_t<I, patch_seq>;
            if (is_in_patch<P>(coord_log, coord)) {
                return static_cast<int>(I);
            }
            return find_patch<I + 1>(coord_log, coord);
        }
    }

    /**
     * @brief Check if the physical coordinate is located inside the patch P.
     * The coordinate transformation is inverted using the Newton algorithm.
     * The input logical cooordinate is modified into the computed associated 
     * logical coordinate of the given physical coordinate.   
     */
    template <class P>
    KOKKOS_INLINE_FUNCTION bool is_in_patch(
            Coord<DimArg1, DimArg2>& logical_coord,
            Coord<DimRes1, DimRes2> const& physical_coord) const
    {
        Coord<DimArg1> dim1_min = ddc::discrete_space<typename P::BSplines1>().rmin();
        Coord<DimArg1> dim1_max = ddc::discrete_space<typename P::BSplines1>().rmax();

        Coord<DimArg2> dim2_min = ddc::discrete_space<typename P::BSplines2>().rmin();
        Coord<DimArg2> dim2_max = ddc::discrete_space<typename P::BSplines2>().rmax();

        replace_inside<P>(logical_coord, dim1_min, dim1_max, dim2_min, dim2_max);
        Coord<DimRes1, DimRes2> computed_physical_coord
                = m_coord_transformations.template get<P>()(logical_coord);

        InverseJacobianMatrix get_inverse_jacobian(m_coord_transformations.template get<P>());

        std::size_t counter = 0;
        while (!(norm_inf(computed_physical_coord - physical_coord) <= m_TOL)
               && (counter < m_counter_max)) {
            counter++;
            DTensor<ArgBasis, VectorIndexSet<typename DimRes1::Dual, typename DimRes2::Dual>> J_inv
                    = get_inverse_jacobian(logical_coord);

            DVector<DimRes1, DimRes2> vec(computed_physical_coord - physical_coord);
            DVector<DimArg1, DimArg2> delta = tensor_mul(index<'i', 'j'>(J_inv), index<'j'>(vec));
            Coord<DimArg1, DimArg2>
                    delta_coord(ddcHelper::get<DimArg1>(delta), ddcHelper::get<DimArg2>(delta));
            logical_coord = logical_coord - delta_coord;

            // Replace the coordinate inside the domain of definition.
            replace_inside<P>(logical_coord, dim1_min, dim1_max, dim2_min, dim2_max);

            computed_physical_coord = m_coord_transformations.template get<P>()(logical_coord);
        }
        return (norm_inf(computed_physical_coord - physical_coord) <= m_TOL);
    }

    /// @brief Replace the logical coordinate into the logical domain.
    template <class P>
    KOKKOS_INLINE_FUNCTION void replace_inside(
            Coord<DimArg1, DimArg2>& coord_log,
            Coord<DimArg1> dim1_min,
            Coord<DimArg1> dim1_max,
            Coord<DimArg2> dim2_min,
            Coord<DimArg2> dim2_max) const
    {
        constexpr bool dim1_periodic = P::BSplines1::is_periodic();
        constexpr bool dim2_periodic = P::BSplines2::is_periodic();

        Coord<DimArg1> coord1_log(coord_log);
        Coord<DimArg2> coord2_log(coord_log);

        if (coord1_log > dim1_max) {
            coord1_log = dim1_periodic ? coord1_log - (dim1_max - dim1_min) : dim1_max;
        }
        if (coord1_log < dim1_min) {
            coord1_log = dim1_periodic ? coord1_log + (dim1_max - dim1_min) : dim1_min;
        }
        if (coord2_log > dim2_max) {
            coord2_log = dim2_periodic ? coord2_log - (dim2_max - dim2_min) : dim2_max;
        }
        if (coord2_log < dim2_min) {
            coord2_log = dim2_periodic ? coord2_log + (dim2_max - dim2_min) : dim2_min;
        }
        coord_log = Coord<DimArg1, DimArg2>(coord1_log, coord2_log);
    }
};


// To help the template deduction.
template <class MultipatchMapping>
PatchLocatorNewton(MultipatchMapping const& coord_transformations)
        -> PatchLocatorNewton<MultipatchMapping>;

template <class MultipatchMapping>
PatchLocatorNewton(
        MultipatchMapping const& coord_transformations,
        std::size_t const counter_max,
        double const TOL) -> PatchLocatorNewton<MultipatchMapping>;
