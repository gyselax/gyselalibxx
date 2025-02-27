// SPDX-License-Identifier: MIT

#pragma once
#include <stdexcept>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

/**
 * @brief Define null extrapolation rule common to all geometries. 
 * @see MultipatchSplineEvaluator2D.
 */
struct NullExtrapolationRule
{
    /**
     * @brief Evaluate to zero the splines at a given coordinate outside of the splines domain. 
     * @param[in] coord_extrap Coordinate where we want to evaluate. 
     * @param[in] patches_splines Splines stored in a MultipatchType. 
     * @param[in] out_of_bounds_idx Index of the localisation of the coordinate. 
     *      It is supposed to be negative to be considered as outside of the domain.
     * @return A null double. 
     */
    template <class... Dim, template <typename P> typename SplinesOnPatch, class... Patches>
    KOKKOS_FUNCTION double operator()(
            Coord<Dim...> const& coord_extrap,
            MultipatchType<SplinesOnPatch, Patches...> const& patches_splines,
            int const out_of_bounds_idx) const
    {
        assert(out_of_bounds_idx < 0);
        return 0.0;
    }
};
