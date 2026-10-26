// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

namespace detail {

/**
 * @brief Not all simulations expose theta,r,phi (tor1, tor2, tor3) so we "spoof" the presence of the dimensions. This
 * is the class from which fake dimensions inherit. These fake dimensions are inserted in order to match Koliop's
 * interface. These fake dimensions are then assigned to an: IdxRange<Grid1D>(Idx<Grid1D> {0}, IdxStep<Grid1D> {1}); See
 * create_scalar_index_range();
 */
struct InternalSpoofGrid
{
};

/**
 * @brief Fake toroidal dimension to be used if there is no toroidal dimension in the simulation.
 */
struct InternalSpoofGridTor3 : InternalSpoofGrid
{
};

/**
 * @brief Fake radial dimension to be used if there is no radial dimension in the simulation.
 */
struct InternalSpoofGridTor2 : InternalSpoofGrid
{
};

/**
 * @brief Fake poloidal dimension to be used if there is no poloidal dimension in the simulation.
 */
struct InternalSpoofGridTor1 : InternalSpoofGrid
{
};

/**
 * @brief A storage class for the fields the collision operator requires. For a tokamak geometry, tor1 is expected to be
 * theta, tor2 r and tor3 phi.
 */
template <
        class IdxRangeDistributionFunction,
        class Species,
        class GridTor3,
        class GridTor2,
        class GridTor1,
        class GridVpar,
        class GridMu>
class CollisionConfigurationData
{
public:
    static_assert(
            ddc::is_uniform_point_sampling_v<GridVpar>,
            "The grid must be uniform in the vpar direction.");
    static_assert(
            ddc::is_uniform_point_sampling_v<GridVpar>,
            "The grid must be uniform in the mu direction.");

public:
    /**
     * @brief Construct a new Collision Configuration Data object/
     *
     * @param[in] yaml_input_file simulation input file.
     * @param[in] index_range_fdistribution the index range that will represent the distribution function given to the
     * operator.
     * @param[in] coeff_intdmu quadrature coefficient.
     * @param[in] coeff_intdvpar quadrature coefficient.
     */
    CollisionConfigurationData(
            PC_tree_t const& yaml_input_file,
            IdxRangeDistributionFunction index_range_fdistribution,
            DConstField<IdxRange<GridMu>> coeff_intdmu,
            DConstField<IdxRange<GridVpar>> coeff_intdvpar)
        : m_hat_As_allocation {extract_md_index_range<Species>(index_range_fdistribution)}
        , m_hat_Zs_allocation {extract_md_index_range<Species>(index_range_fdistribution)}
        , m_B_norm_allocation {extract_md_index_range<GridTor3, GridTor1, GridTor2>(
                  index_range_fdistribution)}
        , m_Bstar_s_allocation {extract_md_index_range<
                  Species,
                  GridTor3,
                  GridTor1,
                  GridTor2,
                  GridVpar>(index_range_fdistribution)}
        , m_mug_allocation {extract_md_index_range<GridMu>(index_range_fdistribution)}
        , m_vparg_allocation {extract_md_index_range<GridVpar>(index_range_fdistribution)}
        , m_coeff_intdmu {coeff_intdmu}
        , m_coeff_intdvpar {coeff_intdvpar}
        , m_coeff_AD_allocation {extract_md_index_range<GridTor2>(index_range_fdistribution)}
        , m_nu_ref {::PCpp_double(yaml_input_file, ".Collisions.nu_ref")}
        , m_collisions_interspecies {::PCpp_bool(yaml_input_file, ".Collisions.interspecies")}
        , m_mu_extent {extract_md_index_range<GridMu>(index_range_fdistribution).size()}
        , m_vpar_extent {extract_md_index_range<GridVpar>(index_range_fdistribution).size()}
        , m_theta_extent {extract_md_index_range<GridTor1>(index_range_fdistribution).size()}
        , m_r_extent {extract_md_index_range<GridTor2>(index_range_fdistribution).size()}
        , m_phi_extent {extract_md_index_range<GridTor3>(index_range_fdistribution).size()}
        , m_species_extent {extract_md_index_range<Species>(index_range_fdistribution).size()}
    {
        /* Define default array values. The user of CollisionConfigurationData is free to modify these.
         */
        {
            IdxRange<Species> const index_range_sp
                    = extract_md_index_range<Species>(index_range_fdistribution);

            /* Initialise the mass species.
             */
            host_t<DConstFieldSp> hat_As_host
                    = ddc::host_discrete_space<Species>().masses()[index_range_sp];
            ddc::parallel_deepcopy(get_field(m_hat_As_allocation), hat_As_host);

            /* Initialise the species charge.
             */
            host_t<DConstFieldSp> hat_Zs_host
                    = ddc::host_discrete_space<Species>().charges()[index_range_sp];
            ddc::parallel_deepcopy(get_field(m_hat_Zs_allocation), hat_Zs_host);
        }
        {
            ddcHelper::
                    dump_coordinates(Kokkos::DefaultExecutionSpace(), get_field(m_mug_allocation));
            ddcHelper::dump_coordinates(
                    Kokkos::DefaultExecutionSpace(),
                    get_field(m_vparg_allocation));
        }

        Kokkos::fence();
    }

    /**
     * @brief We have some fields whose dimension vary depending on the simulation. Bstar_s in a typical 5D simulation
     * has dimension [Sp, theta, R, Vpar]. Now some simulation do not have all these dimension, say only [Sp, Vpar, Mu].
     * We have to satisfy both situations. TransposeDeepCopy takes two fields, expecting the source field to be copied
     * into the destination field. It may be that the source has less dimension but never more. If the source has less
     * dimensions, the data is replicated. For instance from [Vpar, Sp] into [Sp, R, Vpar], the [Vpar, Sp] plane with be
     * replicated along Rand the transpose done appropriately.
     *
     * @param[out] destination a field.
     * @param[in] source a field.
     */
    template <class IdxRange0, class IdxRange1, class IdxRangeSlice = IdxRange1>
    static void transpose_replicate_deep_copy(
            DField<IdxRange0> destination,
            DConstField<IdxRange1> source)
    {
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(destination),
                KOKKOS_LAMBDA(typename IdxRange0::discrete_element_type destination_index) {
                    const typename IdxRangeSlice::discrete_element_type source_index(
                            destination_index);
                    destination(destination_index) = source(source_index);
                });
    }

    /// @brief No-op conversion. a_view^t^t=a_view, so if we keep the underlying
    /// view pointer and simply permute the indicies, we get a layout left view
    /// without copying.
    /// @tparam View
    /// @param a_view A layout left or right, view.
    /// @return A layout left view.
    template <class View>
    static auto as_layout_left(const View& a_view)
    {
        if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutRight>) {
            return [&]<std::size_t... I>(std::index_sequence<I...>) {
                return Kokkos::View<
                        typename View::data_type,
                        Kokkos::LayoutLeft,
                        typename View::device_type,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>> {
                        a_view.data(),
                        a_view.extent(View::rank - 1 - I)...};
            }(std::make_index_sequence<View::rank> {});
        } else {
            return a_view;
        }
    }

    /// @brief Cast from layout right with potentially no all dimension (missing r, theta or phi) into a full blown 6D
    /// array with the proper extents.
    /// @tparam View
    /// @param a_view A distribution function view.
    /// @return A views in 6D, layout left and proper extents.
    template <class View>
    auto cast_fdistrib(const View& a_view) const
    {
        return Kokkos::View<
                typename View::value_type******,
                Kokkos::LayoutLeft,
                typename View::device_type,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                a_view.data(),
                m_mu_extent,
                m_vpar_extent,
                m_theta_extent,
                m_r_extent,
                m_phi_extent,
                m_species_extent);
    }

protected:
    /**
     * @brief Check if a dimension is spoofed but is not present in the actual simulation.
     */
    template <class Grid1D>
    static constexpr bool is_spoofed_dimension_v = std::is_base_of_v<InternalSpoofGrid, Grid1D>;

    /**
     * @brief Create a scalar index range object.
     *
     * @tparam Grid1D
     * @return IdxRange<Grid1D>
     */
    template <class Grid1D>
    static IdxRange<Grid1D> create_scalar_index_range()
    {
        return IdxRange<Grid1D>(Idx<Grid1D> {0}, IdxStep<Grid1D> {1});
    }

    /**
     * @brief Get the index range for a specific grid from a multi-D index range. If the dimension is spoofed and does
     * not appear in the multi-D index range then an index range which only iterates over the index(0) is returned.
     *
     * @tparam Grid The tag for the specific grid.
     * @param[in] index_range The multi-D index range.
     * @returns The index range for the specific grid.
     */
    template <class Grid1D, class IdxRange0>
    static IdxRange<Grid1D> extract_1d_index_range(IdxRange0 index_range)
    {
        if constexpr (is_spoofed_dimension_v<Grid1D>) {
            return create_scalar_index_range<Grid1D>();
        } else {
            return ddc::select<Grid1D>(index_range);
        }
    }

    /**
     * @brief Get the index range for specific grid dimensions from a multi-D index range.
     *
     * @tparam Grid1D The tags for the specific grid dimensions.
     * @param[in] index_range The multi-D index range.
     * @returns The index range for the specific grid.
     */
    template <class... Grid1D, class IdxRange0>
    static IdxRange<Grid1D...> extract_md_index_range(IdxRange0 index_range)
    {
        return IdxRange<Grid1D...>(extract_1d_index_range<Grid1D>(index_range)...);
    }

public:
    /**
     * @brief Normalised masses for all species.
     */
    DFieldMem<IdxRange<Species>> m_hat_As_allocation;

    /**
     * @brief Normalised charges for all species.
     */
    DFieldMem<IdxRange<Species>> m_hat_Zs_allocation;

    /**
     * @brief B norm in (phi,theta,r), ordered to favor the stellarator/tokamak geometry distribution.
     */
    DFieldMem<IdxRange<GridTor3, GridTor1, GridTor2>> m_B_norm_allocation;

    /**
     * @brief Bstar_s(species,phi,theta,r,vpar). We should avoid copying it, it can get quite large. This should match what
     * koliop expects but can also be cheaply changed in koliop for whatever fits better the geometry of the simulation.
     */
    DFieldMem<IdxRange<Species, GridTor3, GridTor1, GridTor2, GridVpar>> m_Bstar_s_allocation;

    /**
     * @brief Grid in mu direction.
     */
    DFieldMem<IdxRange<GridMu>> m_mug_allocation;

    /**
     * @brief Grid in vpar direction.
     */
    DFieldMem<IdxRange<GridVpar>> m_vparg_allocation;

    /**
     * @brief Quadrature coefficients in mu direction.
     */
    DConstField<IdxRange<GridMu>> m_coeff_intdmu;

    /**
     * @brief Quadrature coefficients in vpar direction.
     */
    DConstField<IdxRange<GridVpar>> m_coeff_intdvpar;

    /**
     * @brief Radial coefficients of AD, almost equivalent to the collision frequency. We do not reduce this to
     * nu_ref*sqrt(2) because masking on the r profile can be useful, notably because collision can behave badly on the
     * colder edges.
     */
    DFieldMem<IdxRange<GridTor2>> m_coeff_AD_allocation;

    /**
     * @brief Collisionality, not used as is, instead, used to build a radial profile (coeff_AD).
     */
    double m_nu_ref;

    /**
     * @brief Boolean that is equal to true if inter-species collisions are taken into account. Read in the YAML input
     * file.
     */
    bool m_collisions_interspecies;

    /**
     * @brief mu dimension extent.
     */
    std::size_t m_mu_extent;

    /**
     * @brief vpar dimension extent.
     */
    std::size_t m_vpar_extent;

    /**
     * @brief Theta dimension extent.
     */
    std::size_t m_theta_extent;

    /**
     * @brief r dimension extent.
     */
    std::size_t m_r_extent;

    /**
     * @brief Phi dimension extent.
     */
    std::size_t m_phi_extent;

    /**
     * @brief Species count.
     */
    std::size_t m_species_extent;
};
} // namespace detail
