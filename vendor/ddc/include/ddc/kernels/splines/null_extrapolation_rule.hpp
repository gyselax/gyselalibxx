// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Kokkos_Core.hpp>

namespace ddc {

/**
 * @brief A functor describing a null extrapolation boundary value for 1D spline evaluator.
 */
struct NullExtrapolationRule
{
    /**
     * @brief Evaluates the spline at a coordinate outside of the domain.
     *
     * @return A double with the value of the function outside the domain (here, 0.).
     */
    template <class CoordType, class ChunkSpan>
    KOKKOS_FUNCTION double operator()(CoordType, ChunkSpan) const
    {
        return 0.0;
    }
};
} // namespace ddc
