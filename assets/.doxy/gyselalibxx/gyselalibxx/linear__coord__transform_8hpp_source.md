

# File linear\_coord\_transform.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**linear\_coord\_transform.hpp**](linear__coord__transform_8hpp.md)

[Go to the documentation of this file](linear__coord__transform_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "tensor.hpp"

template <class InputDim, class OutputDim, class CoordJacob = Coord<InputDim>>
class LinearCoordTransform
{
public:
    using CoordArg = Coord<InputDim>;
    using CoordResult = Coord<OutputDim>;
    using CoordJacobian = CoordJacob;

private:
    Coord<InputDim> m_reference_point_on_input_dim;
    Coord<OutputDim> m_reference_point_on_output_dim;
    double m_scaling_factor;

public:
    explicit KOKKOS_FUNCTION LinearCoordTransform(
            Coord<InputDim> reference_point_on_input_dim,
            Coord<OutputDim> reference_point_on_output_dim,
            double scaling_factor)
        : m_reference_point_on_input_dim(reference_point_on_input_dim)
        , m_reference_point_on_output_dim(reference_point_on_output_dim)
        , m_scaling_factor(scaling_factor)
    {
    }

    KOKKOS_DEFAULTED_FUNCTION LinearCoordTransform(LinearCoordTransform const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION ~LinearCoordTransform() = default;

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        double arg_dist = coord - m_reference_point_on_input_dim;
        double res_dist = arg_dist * m_scaling_factor;
        return m_reference_point_on_output_dim + res_dist;
    }

    KOKKOS_FUNCTION double jacobian(CoordArg const& coord) const
    {
        return m_scaling_factor;
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<OutputDim>, VectorIndexSet<InputDim>> jacobian_matrix(
            CoordArg const& coord) const
    {
        DTensor<VectorIndexSet<OutputDim>, VectorIndexSet<InputDim>> jacobian_matrix;
        ddcHelper::get<OutputDim, InputDim>(jacobian_matrix) = m_scaling_factor;
        return jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(std::is_same_v<IndexTag1, OutputDim>);
        static_assert(std::is_same_v<IndexTag2, InputDim>);
        return m_scaling_factor;
    }


    KOKKOS_FUNCTION DTensor<VectorIndexSet<InputDim>, VectorIndexSet<OutputDim>>
    inv_jacobian_matrix(CoordArg const& coord) const
    {
        DTensor<VectorIndexSet<InputDim>, VectorIndexSet<OutputDim>> matrix;
        ddcHelper::get<InputDim, OutputDim>(matrix) = 1.0 / m_scaling_factor;
        return matrix;
    }


    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double inv_jacobian_component(CoordArg const& coord) const
    {
        static_assert(std::is_same_v<IndexTag1, InputDim>);
        static_assert(std::is_same_v<IndexTag2, OutputDim>);
        return 1.0 / m_scaling_factor;
    }

    KOKKOS_INLINE_FUNCTION LinearCoordTransform<OutputDim, InputDim> get_inverse_mapping() const
    {
        return LinearCoordTransform<OutputDim, InputDim>(
                m_reference_point_on_output_dim,
                m_reference_point_on_input_dim,
                1.0 / m_scaling_factor);
    }
};
```


