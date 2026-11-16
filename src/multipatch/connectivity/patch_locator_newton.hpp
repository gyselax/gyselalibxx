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
        double const nan = Kokkos::Experimental::quiet_NaN_v<double>;
        Coord<DimArg1, DimArg2> coord_log(nan, nan);
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
            // Each patch starts from the same initial guess (not from the end point of a
            // failed Newton on the previous patch).
            Coord<DimArg1, DimArg2> guess = coord_log;
            if (is_in_patch<P>(guess, coord)) {
                coord_log = guess;
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

        // If the coordinate is not initialised, initialise at the centre of the patch.
        if (Kokkos::isnan(ddc::get<DimArg1>(logical_coord))
            || Kokkos::isnan(ddc::get<DimArg2>(logical_coord))) {
            logical_coord
                    = initial_guess<P>(physical_coord, dim1_min, dim1_max, dim2_min, dim2_max);
        }

        // If the coordinate is outside the patch, replace inside.
        replace_inside<P>(logical_coord);
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
            replace_inside<P>(logical_coord);

            computed_physical_coord = m_coord_transformations.template get<P>()(logical_coord);
        }
        return (norm_inf(computed_physical_coord - physical_coord) <= m_TOL);
    }

    /// @brief Replace the logical coordinate into the logical domain.
    template <class P>
    KOKKOS_INLINE_FUNCTION static void replace_inside(Coord<DimArg1, DimArg2>& coord_log)
    {
        auto const& bsplines1 = ddc::discrete_space<typename P::BSplines1>();
        auto const& bsplines2 = ddc::discrete_space<typename P::BSplines2>();

        double const dim1 = put_inside<P::BSplines1::is_periodic()>(
                ddc::get<DimArg1>(coord_log),
                bsplines1.rmin(),
                bsplines1.rmax());
        double const dim2 = put_inside<P::BSplines2::is_periodic()>(
                ddc::get<DimArg2>(coord_log),
                bsplines2.rmin(),
                bsplines2.rmax());
        coord_log = Coord<DimArg1, DimArg2>(dim1, dim2);
    }

    /**
     * @brief Put a 1D logical coordinate inside [xmin, xmax].
     * Periodic dimension: the coordinate is wrapped (any number of periods).
     * Non-periodic dimension: the coordinate is clamped to the closest bound.
     */
    template <bool Periodic>
    KOKKOS_INLINE_FUNCTION static double put_inside(
            double const x,
            double const xmin,
            double const xmax)
    {
        if constexpr (Periodic) {
            double const L = xmax - xmin;
            return x - L * Kokkos::floor((x - xmin) / L);
        } else {
            return Kokkos::clamp(x, xmin, xmax);
        }
    }

    /**
     * @brief Compute an initial guess for the Newton algorithm by a coarse-to-fine
     * search of the closest point on the patch P.
     *
     * At each level, the mapping is evaluated on a 3x3 stencil of the current cell
     * and the cell is recentred on the closest point and halved.
     *
     * @param physical_coord Physical coordinate we want to locate. 
     * @param dim1_min Minimum value on the first logical dimension of the patch.
     * @param dim1_max Maximum value on the first logical dimension of the patch.
     * @param dim2_min Minimum value on the second logical dimension of the patch.
     * @param dim2_max Maximum value on the second logical dimension of the patch.
     * @param err_rel Relative size of the final cell with respect to the patch.
     */
    template <class P>
    KOKKOS_INLINE_FUNCTION Coord<DimArg1, DimArg2> initial_guess(
            Coord<DimRes1, DimRes2> const& physical_coord,
            double const dim1_min,
            double const dim1_max,
            double const dim2_min,
            double const dim2_max,
            double const err_rel = 1e-2) const
    {
        // Current cell: centre and half-widths.
        double dim1_mid = 0.5 * (dim1_min + dim1_max);
        double dim2_mid = 0.5 * (dim2_min + dim2_max);
        double delta_dim1 = 0.5 * (dim1_max - dim1_min);
        double delta_dim2 = 0.5 * (dim2_max - dim2_min);

        int const n_levels = static_cast<int>(Kokkos::ceil(Kokkos::log2(1. / err_rel)));
        for (int level = 0; level < n_levels; ++level) {
            double best_dist = Kokkos::Experimental::infinity_v<double>;
            double best1 = dim1_mid;
            double best2 = dim2_mid;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    double dim1 = Kokkos::clamp(dim1_mid + i * delta_dim1, dim1_min, dim1_max);
                    double dim2 = Kokkos::clamp(dim2_mid + j * delta_dim2, dim2_min, dim2_max);
                    Coord<DimRes1, DimRes2> computed_physical_coord
                            = m_coord_transformations.template get<P>()(
                                    Coord<DimArg1, DimArg2>(dim1, dim2));
                    double const dist = norm_inf(computed_physical_coord - physical_coord);
                    if (dist < best_dist) {
                        best_dist = dist;
                        best1 = dim1;
                        best2 = dim2;
                    }
                }
            }
            dim1_mid = best1;
            dim2_mid = best2;
            delta_dim1 *= 0.5;
            delta_dim2 *= 0.5;
        }
        return Coord<DimArg1, DimArg2>(dim1_mid, dim2_mid);
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
