// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "tensor.hpp"

/**
 * @brief A class describing a linear coordinate transformation.
 *
 * A class describing a linear coordinate transformation of the form:
 * @f$ x_1 = \alpha x_2 + \beta @f$
 * where @f$ x_1 @f$ and @f$ x_2 @f$ are coordinates in the input and
 * output dimensions, and @f$ \alpha @f$ and @f$ \beta @f$ are coefficients.
 */
template <class InputDim, class OutputDim, class CoordJacob = Coord<InputDim>>
class LinearCoordTransform
{
public:
    /// The type of the argument of the function described by this mapping
    using CoordArg = Coord<InputDim>;
    /// The type of the result of the function described by this mapping
    using CoordResult = Coord<OutputDim>;
    /// The type of the coordinate that can be used to evaluate the Jacobian of this mapping
    using CoordJacobian = CoordJacob;

private:
    Coord<InputDim> m_reference_point_on_input_dim;
    Coord<OutputDim> m_reference_point_on_output_dim;
    double m_scaling_factor;

public:
    /**
     * @brief A constructor for the linear coordinate transformation.
     *
     * The constructor infers the linear coordinate transformation from a
     * reference point provided in both coordinate systems and a scaling factor.
     *
     * The transformation is then defined as:
     * @f$ x_{out} = x_{out}^* + s * (x_{in} - x_{in}^*) @f$
     * where @f$ x_{in} @f$ and @f$ x_{out} @f$ denote the point expressed in the
     * input and output coordinate system, @f$ \cdot^* @f$ denotes the reference
     * point and @f$ s @f$ denotes the scaling factor.
     *
     * @param reference_point_on_input_dim
     *              The reference point expressed in the input coordinate system.
     * @param reference_point_on_output_dim
     *              The reference point expressed in the output coordinate system.
     * @param scaling_factor
     *              The scaling factor describing how distances in the output
     *              coordinate system scale compared to distances in the input
     *              coordinate system.
     */
    explicit KOKKOS_FUNCTION LinearCoordTransform(
            Coord<InputDim> reference_point_on_input_dim,
            Coord<OutputDim> reference_point_on_output_dim,
            double scaling_factor)
        : m_reference_point_on_input_dim(reference_point_on_input_dim)
        , m_reference_point_on_output_dim(reference_point_on_output_dim)
        , m_scaling_factor(scaling_factor)
    {
    }

    /**
     * Copy constructor for LinearCoordTransform.
     * @param[in] other The operator to be copied.
     */
    KOKKOS_DEFAULTED_FUNCTION LinearCoordTransform(LinearCoordTransform const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION ~LinearCoordTransform() = default;

    /**
     * @brief Convert the coordinate on the input dimension to a coordinate on the output dimension.
     *
     * @param[in] coord The coordinate to be converted expressed on the input coordinate system.
     *
     * @return The coordinate expressed on the output coordinate system.
     */
    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        double arg_dist = coord - m_reference_point_on_input_dim;
        double res_dist = arg_dist * m_scaling_factor;
        return m_reference_point_on_output_dim + res_dist;
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(CoordArg const& coord) const
    {
        return m_scaling_factor;
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficient can be given as a scalar with the function jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<OutputDim>, VectorIndexSet<InputDim>> jacobian_matrix(
            CoordArg const& coord) const
    {
        DTensor<VectorIndexSet<OutputDim>, VectorIndexSet<InputDim>> jacobian_matrix;
        ddcHelper::get<OutputDim, InputDim>(jacobian_matrix) = m_scaling_factor;
        return jacobian_matrix;
    }

    /**
     * @brief Compute the (i,j) coefficient of the Jacobian matrix.
     * 
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(std::is_same_v<IndexTag1, OutputDim>);
        static_assert(std::is_same_v<IndexTag2, InputDim>);
        return m_scaling_factor;
    }


    /**
     * @brief Compute full inverse Jacobian matrix.
     *
     * For some computations, we need the complete inverse Jacobian matrix or just the
     * coefficients.
     * The coefficient can be given as a scalar with the function inv_jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The inverse Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<InputDim>, VectorIndexSet<OutputDim>>
    inv_jacobian_matrix(CoordArg const& coord) const
    {
        DTensor<VectorIndexSet<InputDim>, VectorIndexSet<OutputDim>> matrix;
        ddcHelper::get<InputDim, OutputDim>(matrix) = 1.0 / m_scaling_factor;
        return matrix;
    }


    /**
     * @brief Compute the (i,j) coefficient of the inverse Jacobian matrix.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the inverse Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double inv_jacobian_component(CoordArg const& coord) const
    {
        static_assert(std::is_same_v<IndexTag1, InputDim>);
        static_assert(std::is_same_v<IndexTag2, OutputDim>);
        return 1.0 / m_scaling_factor;
    }

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    KOKKOS_INLINE_FUNCTION LinearCoordTransform<OutputDim, InputDim> get_inverse_mapping() const
    {
        return LinearCoordTransform<OutputDim, InputDim>(
                m_reference_point_on_output_dim,
                m_reference_point_on_input_dim,
                1.0 / m_scaling_factor);
    }
};
