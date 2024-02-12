#pragma once

namespace ddc {

struct NullExtrapolationRule
{
    template <class CoordType, class ChunkSpan>
    KOKKOS_FUNCTION double operator()(CoordType, ChunkSpan) const
    {
        return 0.0;
    }
};
} // namespace ddc
