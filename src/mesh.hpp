#pragma once

#include <array>
#include <memory>

#include <experimental/mdspan>

#include "realdomain.hpp"

template <int N>
class MeshND;

template <>
class MeshND<1> {
private:
    std::vector<double> m_interp_pts;

public:
    RealDomain1D& domain;

    std::experimental::mdspan<double, std::experimental::dynamic_extent> interp_pts();
};

using Mesh1D = MeshND<1>;

using Mesh2D = MeshND<2>;

template <int N>
class MeshND {

    std::array<Mesh1D, N> m_submeshs;

    std::array<std::experimental::mdspan<double, std::experimental::dynamic_extent>, N> interp_pts();
};
