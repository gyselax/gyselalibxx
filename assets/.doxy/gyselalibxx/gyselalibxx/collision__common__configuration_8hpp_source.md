

# File collision\_common\_configuration.hpp

[**File List**](files.md) **>** [**collisions**](dir_64163437c27c8707f17f92558da22106.md) **>** [**collision\_common\_configuration.hpp**](collision__common__configuration_8hpp.md)

[Go to the documentation of this file](collision__common__configuration_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

namespace detail {

struct InternalSpoofGrid
{
};

struct InternalSpoofGridR : InternalSpoofGrid
{
};

struct InternalSpoofGridTheta : InternalSpoofGrid
{
};

struct InternalSpoofGridPhi : InternalSpoofGrid
{
};

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
    using IdxRangeDistributionFunctionType = IdxRangeDistributionFunction;

    using GridMuType = GridMu;
    using GridVparType = GridVpar;
    using GridRType = GridR;
    using GridThetaType = GridTheta;
    using GridPhiType = GridPhi;
    using GridSpType = GridSp;

    using IdxRangeSpThetaRVparType = IdxRange<GridSpType, GridThetaType, GridRType, GridVparType>;
    using IdxRangeThetaRType = IdxRange<GridThetaType, GridRType>;

    using IdxRangeMuType = IdxRange<GridMuType>;
    using IdxRangeVparType = IdxRange<GridVparType>;
    using IdxRangeRType = IdxRange<GridRType>;
    using IdxRangeSpType = IdxRange<GridSpType>;

    static_assert(
            ddc::is_uniform_point_sampling_v<GridVparType>,
            "The grid must be uniform in the vpar direction.");
    static_assert(
            ddc::is_uniform_point_sampling_v<GridVparType>,
            "The grid must be uniform in the mu direction.");

public:
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
        if (!((m_r_extent > 0) && ((m_r_extent & (m_r_extent - 1)) == 0) && (m_theta_extent > 0)
              && ((m_theta_extent & (m_theta_extent - 1)) == 0) && (m_phi_extent > 0)
              && ((m_phi_extent & (m_phi_extent - 1)) == 0))) {
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
    template <class Grid1D>
    static constexpr bool is_spoofed_dimension_v = std::is_base_of_v<InternalSpoofGrid, Grid1D>;

    template <class Grid1D>
    static IdxRange<Grid1D> create_scalar_index_range()
    {
        return IdxRange<Grid1D>(Idx<Grid1D> {0}, IdxStep<Grid1D> {1});
    }

    template <class Grid1D, class IdxRange0>
    static IdxRange<Grid1D> extract_1d_index_range(IdxRange0 index_range)
    {
        if constexpr (is_spoofed_dimension_v<Grid1D>) {
            return create_scalar_index_range<Grid1D>();
        } else {
            return ddc::select<Grid1D>(index_range);
        }
    }

    template <class... Grid1D, class IdxRange0>
    static IdxRange<Grid1D...> extract_md_index_range(IdxRange0 index_range)
    {
        return IdxRange<Grid1D...>(extract_1d_index_range<Grid1D>(index_range)...);
    }

    /* Not protected nor private because there is no internal state to
     * protect/condition to satisfy.
     */
public:
    DFieldMem<IdxRangeSpType> m_hat_As_allocation;

    DFieldMem<IdxRangeSpType> m_hat_Zs_allocation;

    DFieldMem<IdxRangeRType> m_mask_buffer_r_allocation;

    DFieldMem<IdxRangeThetaRType> m_mask_LIM_allocation;

    DFieldMem<IdxRangeThetaRType> m_B_norm_allocation;

    DFieldMem<IdxRangeSpThetaRVparType> m_Bstar_s_allocation;

    DFieldMem<IdxRangeMuType> m_mug_allocation;

    DFieldMem<IdxRangeVparType> m_vparg_allocation;

    DFieldMem<IdxRangeRType> m_coeff_AD_allocation;

    DConstField<IdxRangeMuType> m_coeff_intdmu;

    DConstField<IdxRangeVparType> m_coeff_intdvpar;

    bool m_collisions_interspecies;

    std::size_t m_mu_extent;

    std::size_t m_vpar_extent;

    std::size_t m_r_extent;

    std::size_t m_theta_extent;

    std::size_t m_phi_extent;

    std::size_t m_sp_extent;
};
} // namespace detail
```


