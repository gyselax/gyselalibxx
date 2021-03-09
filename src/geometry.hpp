#pragma once

#include <array>
#include <memory>
#include <vector>

#include <bsplines.h>

#include "realdomain.hpp"

template <int N>
class GeometryND;

template <>
class GeometryND<1> {
public:
    std::unique_ptr<BSplines> m_spline;

    std::array<std::reference_wrapper<BSplines>, 1> splines() const;
};

using Geometry1D = GeometryND<1>;

using Geometry2D = GeometryND<2>;

template <int N>
class GeometryND {
public:
    std::array<Geometry1D, N> m_subgeom;

    std::array<std::reference_wrapper<BSplines>, N> splines() const;
};
