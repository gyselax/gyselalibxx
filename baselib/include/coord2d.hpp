#pragma once

// standard C++ library headers
#include <array>

/// The dimension identifiers
enum Dimension2D : int { DY = 0, DX = 1 };

/// a 2D coordinate
using Coord2D = std::array<int, 2>;

/// a 2D shape a.k.a. extent, i.e. a ND generalization of size
using Shape2D = std::array<int, 2>;

/// A 2D direction
enum Direction2D { LEFT, RIGHT, UP, DOWN };
