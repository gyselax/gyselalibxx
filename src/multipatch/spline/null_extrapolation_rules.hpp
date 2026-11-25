// SPDX-License-Identifier: MIT

#pragma once
#include <stdexcept>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

/**
 * @brief Define null extrapolation rule common to all geometries. 
 * @see InterpolationEvaluator2DCollection.
 */
struct NullExtrapolationRule
{
    /**
     * @brief Evaluate to zero the splines at a given coordinate outside of the splines domain. 
     * @param[in] coord_extrap Coordinate where we want to evaluate. 
     * @param[in] patches_splines Splines stored in a MultipatchType. The coefficients on each
     *      patch are always defined on the unbatched coefficient index range of the patch.
     *      When the evaluated values are batched, the InterpolationEvaluator2DCollection
     *      provides the batch slice of the coefficients corresponding to the evaluated point.
     * @param[in] out_of_bounds_idx Index of the localisation of the coordinate. 
     *      It is supposed to be negative to be considered as outside of the domain.
     * @return A null double. 
     */
    template <class... Dim, class TypeSeqPatches, class TypeSeqInternalTypes>
    KOKKOS_FUNCTION double operator()(
            Coord<Dim...> const& coord_extrap,
            detail::MultipatchType<TypeSeqPatches, TypeSeqInternalTypes> const& patches_splines,
            int const out_of_bounds_idx) const
    {
        assert(out_of_bounds_idx < 0);
        return 0.0;
    }
};
