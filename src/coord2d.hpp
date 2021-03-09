#pragma once

// standard C++ library headers
#include <array>

/// The dimension identifiers
enum Dimension : int {
    DVX = 0,
    DX = 1
};

/// a 2D coordinate
using Coord2D = std::array<int, 2>;
template <int N>
using Coord = std::array<int, N>;

/// a 2D shape a.k.a. extent, i.e. a ND generalization of size
using Shape2D = std::array<int, 2>;
template <int N>
using Shape = std::array<int, N>;
