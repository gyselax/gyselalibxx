#pragma once

#include <array>

namespace Dim {

struct X {
    static constexpr bool PERIODIC = true;
};

struct VX {
    static constexpr bool PERIODIC = false;
};

struct Time {
    static constexpr bool PERIODIC = false;
};

} // namespace Dim

struct RSpace;

template <int NDIMS>
using RCoordND = std::array<double, NDIMS>;

using RCoord1D = RCoordND<1>;

using RCoord2D = RCoordND<2>;

class RDimension
{
    /** Periodicity of the dimension
     *
     * - in case the dimension is periodic, m_periodicity != 0.
     *   && normalize(x+n.m_periodicity)== x,
     * - otherwise, m_periodicity==0. and is unused
     */
    RCoord1D m_periodicity;

public:
    RCoord1D normalize(RCoord1D) const;
};

template <int NDIMS>
struct RDomainND {
    const RCoordND<NDIMS> start;

    const RCoordND<NDIMS> end;

    constexpr RDomainND(RCoordND<NDIMS> start, RCoordND<NDIMS> end) noexcept
        : start(start)
        , end(end)
    {
    }
};

using RDomain1D = RDomainND<1>;

using RDomain2D = RDomainND<2>;
