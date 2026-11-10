

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

struct InternalSpoofGrid : ddc::DiscreteDimension
{
};

struct InternalSpoofGridTor3 : InternalSpoofGrid
{
};

struct InternalSpoofGridTor2 : InternalSpoofGrid
{
};

struct InternalSpoofGridTor1 : InternalSpoofGrid
{
};

template <
        class IdxRangeDistributionFunction,
        class GridTor3Type,
        class GridTor2Type,
        class GridTor1Type,
        class GridVparType,
        class GridMuType>
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
    CollisionConfigurationData(
            PC_tree_t const& yaml_input_file,
            IdxRangeDistributionFunction index_range_fdistribution,
            DConstField<IdxRange<GridMu>> coeff_intdmu,
            DConstField<IdxRange<GridVpar>> coeff_intdvpar)
        : m_hat_As_allocation {extract_md_index_range<Species>(index_range_fdistribution)}
        , m_hat_Zs_allocation {extract_md_index_range<Species>(index_range_fdistribution)}
        , m_B_norm_allocation {extract_md_index_range<GridTor3Type, GridTor1Type, GridTor2Type>(
                  index_range_fdistribution)}
        , m_Bstar_s_allocation {extract_md_index_range<
                  Species,
                  GridTor3Type,
                  GridTor1Type,
                  GridTor2Type,
                  GridVpar>(index_range_fdistribution)}
        , m_mug_allocation {extract_md_index_range<GridMu>(index_range_fdistribution)}
        , m_vparg_allocation {extract_md_index_range<GridVpar>(index_range_fdistribution)}
        , m_coeff_intdmu {coeff_intdmu}
        , m_coeff_intdvpar {coeff_intdvpar}
        , m_coeff_AD_allocation {extract_md_index_range<GridTor1Type>(index_range_fdistribution)}
        , m_nu_ref {::PCpp_double(yaml_input_file, ".Collisions.nu_ref")}
        , m_collisions_interspecies {::PCpp_bool(yaml_input_file, ".Collisions.interspecies")}
        , m_mu_extent {extract_md_index_range<GridMu>(index_range_fdistribution).size()}
        , m_vpar_extent {extract_md_index_range<GridVpar>(index_range_fdistribution).size()}
        , m_theta_extent {extract_md_index_range<GridTor2Type>(index_range_fdistribution).size()}
        , m_r_extent {extract_md_index_range<GridTor1Type>(index_range_fdistribution).size()}
        , m_phi_extent {extract_md_index_range<GridTor3Type>(index_range_fdistribution).size()}
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

    template <class View>
    static auto as_layout_left(const View& a_view)
    {
        if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutRight>) {
            return [&]<std::size_t... I>(std::index_sequence<I...>)
            {
                return Kokkos::View<
                        typename View::data_type,
                        Kokkos::LayoutLeft,
                        typename View::device_type,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>> {
                        a_view.data(),
                        a_view.extent(View::rank - 1 - I)...};
            }
            (std::make_index_sequence<View::rank> {});
        } else {
            return a_view;
        }
    }

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

public:
    DFieldMem<IdxRange<Species>> m_hat_As_allocation;

    DFieldMem<IdxRange<Species>> m_hat_Zs_allocation;

    DFieldMem<IdxRange<GridTor3Type, GridTor1Type, GridTor2Type>> m_B_norm_allocation;

    DFieldMem<IdxRange<Species, GridTor3Type, GridTor1Type, GridTor2Type, GridVpar>>
            m_Bstar_s_allocation;

    DFieldMem<IdxRange<GridMu>> m_mug_allocation;

    DFieldMem<IdxRange<GridVpar>> m_vparg_allocation;

    DConstField<IdxRange<GridMu>> m_coeff_intdmu;

    DConstField<IdxRange<GridVpar>> m_coeff_intdvpar;

    DFieldMem<IdxRange<GridTor1Type>> m_coeff_AD_allocation;

    double m_nu_ref;

    bool m_collisions_interspecies;

    std::size_t m_mu_extent;

    std::size_t m_vpar_extent;

    std::size_t m_theta_extent;

    std::size_t m_r_extent;

    std::size_t m_phi_extent;

    std::size_t m_species_extent;
};
} // namespace detail
```


