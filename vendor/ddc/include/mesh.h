#pragma once

#include <type_traits>

/// @class Mesh
/// @brief `Mesh` should be used as a base class of classes actually modelling what we consider to be a "mesh". The derived class should implement the method `to_real` and `operator==`.
struct Mesh
{
};

template <class T>
inline constexpr bool is_mesh_v = std::is_base_of_v<Mesh, T>;
