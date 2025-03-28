

# File onion\_patch\_locator.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**onion\_patch\_locator.hpp**](onion__patch__locator_8hpp.md)

[Go to the documentation of this file](onion__patch__locator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <stdexcept>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "multipatch_type.hpp"
#include "types.hpp"



template <
        class MultipatchIdxRanges,
        class LogicalToPhysicalMapping,
        class PhysicalToLogicalMapping,
        class ExecSpace = Kokkos::DefaultExecutionSpace>
class OnionPatchLocator;


template <
        class... Patches,
        class LogicalToPhysicalMapping,
        class PhysicalToLogicalMapping,
        class ExecSpace>
class OnionPatchLocator<
        MultipatchType<IdxRangeOnPatch, Patches...>,
        LogicalToPhysicalMapping,
        PhysicalToLogicalMapping,
        ExecSpace>
{
    static_assert(is_curvilinear_2d_mapping_v<LogicalToPhysicalMapping>);

    using X = typename LogicalToPhysicalMapping::cartesian_tag_x;
    using Y = typename LogicalToPhysicalMapping::cartesian_tag_y;
    using R = typename LogicalToPhysicalMapping::curvilinear_tag_r;
    using Theta = typename LogicalToPhysicalMapping::curvilinear_tag_theta;

    static_assert(Theta::PERIODIC, "Theta dimension must be periodic.");

public:
    using MultipatchIdxRanges = MultipatchType<IdxRangeOnPatch, Patches...>;

    using PatchOrdering = ddc::detail::TypeSeq<Patches...>;

    using exec_space = ExecSpace;

    static constexpr int outside_rmax_domain = -1;

    static constexpr int outside_rmin_domain = -2;

private:
    static constexpr std::size_t n_patches = ddc::type_seq_size_v<PatchOrdering>;

    static_assert(
            std::is_invocable_r_v<Coord<R, Theta>, PhysicalToLogicalMapping, Coord<X, Y>>,
            "The mapping has to contain an operator from the physical domain to the logical "
            "domain.");
    static_assert(
            (((std::is_same_v<typename Patches::Dim1, R>)&&(
                     std::is_same_v<typename Patches::Dim2, Theta>))
             && ...),
            "The mappings and the patches have to be defined on the same dimensions.");


    LogicalToPhysicalMapping const m_to_cartesian_mapping;
    PhysicalToLogicalMapping const m_to_curvilinear_mapping;
    MultipatchIdxRanges const m_all_idx_ranges;

    Kokkos::View<Coord<R>*, typename ExecSpace::memory_space> m_radii;

public:
    OnionPatchLocator(
            MultipatchIdxRanges const& all_idx_ranges,
            LogicalToPhysicalMapping const& to_physical_mapping,
            PhysicalToLogicalMapping const& to_logical_mapping)
        : m_to_cartesian_mapping(to_physical_mapping)
        , m_to_curvilinear_mapping(to_logical_mapping)
        , m_all_idx_ranges(all_idx_ranges)
        , m_radii("m_radii", n_patches + 1)
    {
        set_and_check_radii();
    }


    ~OnionPatchLocator() = default;

    KOKKOS_INLINE_FUNCTION int operator()(Coord<X, Y> const coord) const
    {
        int patch_index_min = 0;
        int patch_index_max = n_patches - 1;

        Coord<R> r_min = m_radii(patch_index_min);
        Coord<R> r_max = m_radii(patch_index_max + 1);
        Coord<R> r(m_to_curvilinear_mapping(coord));
        KOKKOS_ASSERT(!Kokkos::isnan(double(r)));

        if (r > r_max) {
            return outside_rmax_domain;
        }
        if (r < r_min) {
            return outside_rmin_domain;
        }

        while (patch_index_max >= patch_index_min) {
            int patch_index_mid = (patch_index_min + patch_index_max) / 2;
            r_min = m_radii(patch_index_mid);
            r_max = m_radii(patch_index_mid + 1);

            if (r_min <= r && r < r_max) {
                return patch_index_mid;
            } else if (r < r_min) {
                patch_index_max = patch_index_mid - 1;
            } else if (r_max <= r) {
                patch_index_min = patch_index_mid + 1;
            }
        };

        if (r != m_radii(n_patches)) {
            Kokkos::abort("Dichotomy method failed to find the right patch.");
        }
        return n_patches - 1;
    }


    template <class Patch>
    KOKKOS_FUNCTION LogicalToPhysicalMapping get_mapping_on_patch() const
    {
        return m_to_cartesian_mapping;
    }

    template <class Dim1, class Dim2>
    KOKKOS_FUNCTION LogicalToPhysicalMapping get_mapping_on_logical_dim() const
    {
        static_assert(
                (std::is_same_v<Dim1, R>)&&(std::is_same_v<Dim2, Theta>),
                "Wrong continuous dimensions.");
        return m_to_cartesian_mapping;
    }

    template <class Patch>
    using get_mapping_on_patch_t = LogicalToPhysicalMapping;

    template <class Dim1, class Dim2>
    using get_mapping_on_logical_dim_t = LogicalToPhysicalMapping;

private:
    void set_and_check_radii()
    {
        Kokkos::View<Coord<R>*, Kokkos::HostSpace> radii_host("radii_host", n_patches + 1);

        std::array<Coord<R>, n_patches> r_min {
                (Coord<R>(ddc::coordinate(m_all_idx_ranges.template get<Patches>().front())))...};
        std::array<Coord<R>, n_patches> r_max {
                (Coord<R>(ddc::coordinate(m_all_idx_ranges.template get<Patches>().back())))...};

        radii_host(0) = r_min[0];
        for (std::size_t i(0); i < n_patches - 1; i++) {
            if (abs(double(r_min[i + 1] - r_max[i])) > 1e-14) {
                throw std::invalid_argument("The patches listed in PatchOrdering must be ordered.");
            }
            radii_host(i + 1) = r_max[i];
        }
        radii_host(n_patches) = r_max[n_patches - 1];

        Kokkos::deep_copy(m_radii, radii_host);
    }
};


// To help the template deduction.
template <
        class MultipatchIdxRanges,
        class LogicalToPhysicalMapping,
        class PhysicalToLogicalMapping,
        class ExecSpace = Kokkos::DefaultExecutionSpace>
OnionPatchLocator(
        MultipatchIdxRanges const& all_idx_ranges,
        LogicalToPhysicalMapping const& to_physical_mapping,
        PhysicalToLogicalMapping const& to_logical_mapping)
        -> OnionPatchLocator<
                MultipatchIdxRanges,
                LogicalToPhysicalMapping,
                PhysicalToLogicalMapping,
                ExecSpace>;
```


