// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "onion_patch_locator.hpp"
#include "utils_patch_locators.hpp"


/**
 * @brief Define constant extrapolation rule for onion shape geometries. 
 * Struct useful for the MultipatchSplineEvaluator types.  
 * @tparam PatchLocator A patch locator specialised for onion shape geometries. 
 */
template <class PatchLocator>
struct ConstantExtrapolationRuleOnion
{
private:
    static_assert(
            is_onion_patch_locator_v<PatchLocator>,
            "Extrapolation rule only defined for OnionPatchLocator.");

    using PatchOrdering = typename PatchLocator::PatchOrdering;

    static constexpr std::size_t n_patches = ddc::type_seq_size_v<PatchOrdering>;

    using MinRadiusPatch = ddc::type_seq_element_t<0, PatchOrdering>;
    using MaxRadiusPatch = ddc::type_seq_element_t<n_patches - 1, PatchOrdering>;

    using R_min = typename MinRadiusPatch::Dim1;
    using Theta_min = typename MinRadiusPatch::Dim2;

    using R_max = typename MaxRadiusPatch::Dim1;
    using Theta_max = typename MaxRadiusPatch::Dim2;

    ddc::ConstantExtrapolationRule<R_min, Theta_min> const bc_r_min;
    ddc::ConstantExtrapolationRule<R_max, Theta_max> const bc_r_max;

public:
    /**
     * @brief Instantiate a ConstantExtrapolationRuleOnion.
     * The R1 and R2 templates are needed for GPU.
     * @tparam R1 Continuous dimension of the minimum radial coordinate of the domain.
     * @tparam R2 Continuous dimension of the maximum radial coordinate of the domain.
     * @param[in] r_min Minimum radial coordinate of the domain.
     * @param[in] r_max Maximum radial coordinate of the domain.
     */
    explicit ConstantExtrapolationRuleOnion(Coord<R_min> const& r_min, Coord<R_max> const& r_max)
        : bc_r_min(r_min)
        , bc_r_max(r_max)
    {
    }

    /**
     * @brief Evaluate at a given outside coordinate. 
     * @tparam Dim Continuous dimensions where the given coordinate is defined. 
     * @tparam SplinesOnPatch Field of spline coefficients template on the Patch.
     * @tparam Patches Patch types.
     * @param[in] coord_extrap Coordinate where we want to evaluate. 
     * @param[in] patches_splines Splines stored in a MultipatchType. 
     * @param[in] out_of_bounds_idx Index of the localisation of the coordinate. 
     *      It is supposed to be negative to be considered as outside of the domain.
     * @return A double with the evaluation outside of the domain. 
     */
    template <class... Dim, template <typename P> typename SplinesOnPatch, class... Patches>
    KOKKOS_FUNCTION double operator()(
            Coord<Dim...> const& coord_extrap,
            MultipatchField<SplinesOnPatch, Patches...> const& patches_splines,
            int const out_of_bounds_idx) const
    {
        assert((out_of_bounds_idx == PatchLocator::outside_rmin_domain)
               || (out_of_bounds_idx == PatchLocator::outside_rmax_domain));

        if (out_of_bounds_idx == PatchLocator::outside_rmin_domain) {
            SplinesOnPatch<MinRadiusPatch> const min_spline
                    = patches_splines.template get<MinRadiusPatch>();
            return bc_r_min(coord_extrap, min_spline);
        } else {
            SplinesOnPatch<MaxRadiusPatch> const max_spline
                    = patches_splines.template get<MaxRadiusPatch>();
            return bc_r_max(coord_extrap, max_spline);
        }
    }
};
