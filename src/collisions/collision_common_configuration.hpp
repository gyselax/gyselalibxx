// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

namespace detail {

/**
 * @brief Not all simulations expose r,theta,phi (tor1, tor2, tor3) so we
 * "spoof" the presence of the dimensions. This is the class from which fake
 * dimensions inherit. These fake dimensions are inserted in order to match
 * Koliop's interface.
 * These fake dimensions are then assigned to an:
 *      IdxRange<Grid1D>(Idx<Grid1D> {0}, IdxStep<Grid1D> {1});
 *  See create_scalar_index_range();
 */
struct InternalSpoofGrid
{
};

/**
 * @brief Fake radial dimension to be used if there is no radial dimension in
 * the simulation.
 */
struct InternalSpoofGridR : InternalSpoofGrid
{
};

/**
 * @brief Fake poloidal dimension to be used if there is no poloidal dimension
 * in the simulation.
 */
struct InternalSpoofGridTheta : InternalSpoofGrid
{
};

/**
 * @brief Fake toroidal dimension to be used if there is no toroidal dimension
 * in the simulation.
 */
struct InternalSpoofGridPhi : InternalSpoofGrid
{
};

/**
 * @brief A storage class for the fields the collision operator requires.
 */
template <
        class IdxRangeDistributionFunction,
        class GridSp,
        class GridPhi,
        class GridR,
        class GridTheta,
        class GridVpar,
        class GridMu>
class CollisionConfigurationData
{
public:
    /**
     * @brief Distribution function index range.
     */
    using IdxRangeDistributionFunctionType = IdxRangeDistributionFunction;

    /**
     * @brief Mu grid.
     */
    using GridMuType = GridMu;
    /**
     * @brief Vpar grid.
     */
    using GridVparType = GridVpar;
    /**
     * @brief R grid.
     */
    using GridRType = GridR;
    /**
     * @brief Theta grid.
     */
    using GridThetaType = GridTheta;
    /**
     * @brief Phi grid.
     */
    using GridPhiType = GridPhi;
    /**
     * @brief Sp grid.
     */
    using GridSpType = GridSp;

    /**
     * @brief Sp,Theta,R,Vpar index range.
     */
    using IdxRangeSpThetaRVparType = IdxRange<GridSpType, GridThetaType, GridRType, GridVparType>;
    /**
     * @brief Theta,R index range.
     */
    using IdxRangeThetaRType = IdxRange<GridThetaType, GridRType>;

    /**
     * @brief Mu index range.
     */
    using IdxRangeMuType = IdxRange<GridMuType>;
    /**
     * @brief Vpar index range.
     */
    using IdxRangeVparType = IdxRange<GridVparType>;
    /**
     * @brief R index range.
     */
    using IdxRangeRType = IdxRange<GridRType>;
    /**
     * @brief Species index range.
     */
    using IdxRangeSpType = IdxRange<GridSpType>;

    static_assert(
            ddc::is_uniform_point_sampling_v<GridVparType>,
            "The grid must be uniform in the vpar direction.");
    static_assert(
            ddc::is_uniform_point_sampling_v<GridVparType>,
            "The grid must be uniform in the mu direction.");

public:
    /**
     * @brief Construct a new Collision Configuration Data object/
     *
     * @param[in] yaml_input_file simulation input file.
     * @param[in] index_range_fdistribution the index range that will represent the
     * distribution function given to the operator.
     * @param[in] coeff_intdmu quadrature coefficient.
     * @param[in] coeff_intdvpar quadrature coefficient.
     */
    CollisionConfigurationData(
            PC_tree_t const& yaml_input_file,
            IdxRangeDistributionFunctionType index_range_fdistribution,
            DConstField<IdxRangeMuType> coeff_intdmu,
            DConstField<IdxRangeVparType> coeff_intdvpar)
        : m_hat_As_allocation {extract_md_index_range<GridSpType>(index_range_fdistribution)}
        , m_hat_Zs_allocation {extract_md_index_range<GridSpType>(index_range_fdistribution)}
        , m_mask_buffer_r_allocation {extract_md_index_range<GridRType>(index_range_fdistribution)}
        , m_mask_LIM_allocation {extract_md_index_range<GridThetaType, GridRType>(
                  index_range_fdistribution)}
        , m_B_norm_allocation {extract_md_index_range<GridThetaType, GridRType>(
                  index_range_fdistribution)}
        , m_Bstar_s_allocation {extract_md_index_range<
                  GridSpType,
                  GridThetaType,
                  GridRType,
                  GridVparType>(index_range_fdistribution)}
        , m_mug_allocation {extract_md_index_range<GridMuType>(index_range_fdistribution)}
        , m_vparg_allocation {extract_md_index_range<GridVparType>(index_range_fdistribution)}
        , m_coeff_AD_allocation {extract_md_index_range<GridRType>(index_range_fdistribution)}
        , m_coeff_intdmu {coeff_intdmu}
        , m_coeff_intdvpar {coeff_intdvpar}
        , m_collisions_interspecies {::PCpp_bool(yaml_input_file, ".Collisions.interspecies")}
        , m_mu_extent {extract_md_index_range<GridMuType>(index_range_fdistribution).size()}
        , m_vpar_extent {extract_md_index_range<GridVparType>(index_range_fdistribution).size()}
        , m_r_extent {extract_md_index_range<GridRType>(index_range_fdistribution).size()}
        , m_theta_extent {extract_md_index_range<GridThetaType>(index_range_fdistribution).size()}
        , m_phi_extent {extract_md_index_range<GridPhiType>(index_range_fdistribution).size()}
        , m_sp_extent {extract_md_index_range<GridSpType>(index_range_fdistribution).size()}
    {
        if (!((m_mu_extent > 0) && ((m_mu_extent & (m_mu_extent - 1)) == 0) && (m_vpar_extent > 0)
              && ((m_vpar_extent & (m_vpar_extent - 1)) == 0) && (m_r_extent > 0)
              && ((m_r_extent & (m_r_extent - 1)) == 0) && (m_theta_extent > 0)
              && ((m_theta_extent & (m_theta_extent - 1)) == 0) && (m_phi_extent > 0)
              && ((m_phi_extent & (m_phi_extent - 1)) == 0) && (m_sp_extent > 0)
              && ((m_sp_extent & (m_sp_extent - 1)) == 0))) {
            throw std::runtime_error("For performance reason, the number of points in "
                                     "each dimension must be a power of 2.");
        }

        /* Define default array values. The user of CollisionConfigurationData
         * is free to modify these.
         */

        {
            IdxRangeSpType const index_range_sp
                    = extract_md_index_range<GridSpType>(index_range_fdistribution);

            /* Initialise the mass species.
             */
            host_t<DConstFieldSp> hat_As_host
                    = ddc::host_discrete_space<GridSpType>().masses()[index_range_sp];
            ddc::parallel_deepcopy(get_field(m_hat_As_allocation), hat_As_host);

            /* Initialise the species charge.
             */
            host_t<DConstFieldSp> hat_Zs_host
                    = ddc::host_discrete_space<GridSpType>().charges()[index_range_sp];
            ddc::parallel_deepcopy(get_field(m_hat_Zs_allocation), hat_Zs_host);
        }

        {
            /* Disabling Initialisation of the masks that have no sense here to
             * 0.
             * Masked if >= 0.99
             */
            ddc::parallel_fill(get_field(m_mask_buffer_r_allocation), 0.0);
            ddc::parallel_fill(get_field(m_mask_LIM_allocation), 0.0);
        }

        {
            /* Initialise grid coordinates.
             */
            ddcHelper::
                    dump_coordinates(Kokkos::DefaultExecutionSpace(), get_field(m_mug_allocation));
            ddcHelper::dump_coordinates(
                    Kokkos::DefaultExecutionSpace(),
                    get_field(m_vparg_allocation));
        }

        Kokkos::fence();
    }

    /**
     * @brief We have some fields whose dimension vary depending on the
     * simulation. Bstar_s in a typical 5D simulation has dimension [Sp, theta,
     * R, Vpar]. Now some simulation do not have all these dimension, say only
     * [Sp, Vpar, Mu]. We have to satisfy both situations. TransposeDeepCopy
     * takes two fields, expecting the source field to be copied into the
     * destination field. It may be that the source has less dimension but never
     * more. If the source has less dimensions, the data is replicated. For
     * instance from [Vpar, Sp] into [Sp, R, Vpar], the [Vpar, Sp] plane with be
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
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(destination),
                KOKKOS_LAMBDA(typename IdxRange0::discrete_element_type destination_index) {
                    const typename IdxRangeSlice::discrete_element_type source_index(
                            destination_index);
                    destination(destination_index) = source(source_index);
                });
    }

protected:
    /**
     * @brief Check if a dimension is spoofed but is not present in the actual
     * simulation.
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
     * @brief Get the index range for a specific grid from a multi-D index range.
     * If the dimension is spoofed and does not appear in the multi-D index range
     * then an index range which only iterates over the index(0) is returned.
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
     * @brief Get the index range for specific grid dimensions from a multi-D index
     * range.
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

    /* Not protected nor private because there is no internal state to
     * protect/condition to satisfy.
     */
public:
    /**
     * @brief Normalised masses for all species.
     */
    DFieldMem<IdxRangeSpType> m_hat_As_allocation;

    /**
     * @brief Normalised charges for all species.
     */
    DFieldMem<IdxRangeSpType> m_hat_Zs_allocation;

    /**
     * @brief Mask used to avoid applying collision in certain region.
     * [TODO]: This mask is not exposed to the user in the C++ version.
     */
    DFieldMem<IdxRangeRType> m_mask_buffer_r_allocation;

    /**
     * @brief Limiter mask in (theta,r).
     */
    DFieldMem<IdxRangeThetaRType> m_mask_LIM_allocation;

    /**
     * @brief B norm in (theta,r).
     * [TODO] Attention this must be 3D for generalisation to 3D geometry.
     */
    DFieldMem<IdxRangeThetaRType> m_B_norm_allocation;

    /**
     * @brief Bstar_s(species,theta,r,vpar).
     * [TODO] Must be 5D for full 3D geometry.
     */
    DFieldMem<IdxRangeSpThetaRVparType> m_Bstar_s_allocation;

    /**
     * @brief Grid in mu direction.
     */
    DFieldMem<IdxRangeMuType> m_mug_allocation;

    /**
     * @brief Grid in vpar direction.
     */
    DFieldMem<IdxRangeVparType> m_vparg_allocation;

    /**
     * @brief Radial coefficients of AD, almost equivalent to the collision
     * frequency.
     */
    DFieldMem<IdxRangeRType> m_coeff_AD_allocation;

    /**
     * @brief Quadrature coefficients in mu direction.
     */
    DConstField<IdxRangeMuType> m_coeff_intdmu;

    /**
     * @brief Quadrature coefficients in vpar direction.
     */
    DConstField<IdxRangeVparType> m_coeff_intdvpar;

    /**
     * @brief Boolean that is equal to true if inter-species collisions are
     * taken into account.
     * Read in the YAML input file.
     */
    bool m_collisions_interspecies;

    /**
     * @brief Mu dimension extent.
     */
    std::size_t m_mu_extent;

    /**
     * @brief Vpar dimension extent.
     */
    std::size_t m_vpar_extent;

    /**
     * @brief R dimension extent.
     */
    std::size_t m_r_extent;

    /**
     * @brief Theta dimension extent.
     */
    std::size_t m_theta_extent;

    /**
     * @brief Phi dimension extent.
     */
    std::size_t m_phi_extent;

    /**
     * @brief Species dimension extent.
     */
    std::size_t m_sp_extent;
};
} // namespace detail
