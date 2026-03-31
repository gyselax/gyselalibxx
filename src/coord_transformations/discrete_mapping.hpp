// SPDX-License-Identifier: MIT
#pragma once

#include "i_interpolation_evaluator.hpp"

namespace details {

template <class RowDim, class DataType, class... RDim, class... ADim, class Mapping>
void fill_jacobian_matrix_row(
        Tensor<DataType, VectorIndexSet<RDim...>, VectorIndexSet<ADim...>> jacobian_matrix,
        Mapping mapping)
{
    ((ddc::get<RowDim, ADim>(jacobian_matrix) = mapping.jacobian_component<RowDim, ADim>), ...);
}

template <class DataType, class... RDim, class... ADim, class Mapping>
void fill_jacobian_matrix(
        Tensor<DataType, VectorIndexSet<RDim...>, VectorIndexSet<ADim...>> jacobian_matrix,
        Mapping mapping)
{
    fill_jacobian_matrix_row<RDim>(jacobian_matrix, mapping);
}

} // namespace details

template <class StartCoord, class EndCoord, class NDEvaluator>
class DiscreteMapping
{
    /// The type of the argument of the function described by this mapping
    using CoordArg = StartCoord;
    /// The type of the result of the function described by this mapping
    using CoordResult = EndCoord;
    /// The type of the coordinate that can be used to evaluate the Jacobian of this mapping
    using CoordJacobian = StartCoord;

    using DataType = typename NDEvaluator::data_type;
    using ArgBasis = ddc::to_type_seq_t<StartCoord>;
    using ResultBasis = ddc::to_type_seq_t<EndCoord>;

private:
    using CoeffField = VectorField<
            typename NDEvaluator::DataType,
            typename NDEvaluator::coeff_idx_range_type,
            ddc::to_type_seq_t<CoordResult>>;

private:
    CoeffField m_coeff_representation;
    NDEvaluator m_evaluator;

public:
    /**
     * @brief Instantiate a DiscreteMapping from the coefficients of an interpolating function which
     * approximates the mapping.
     *
     * A discrete mapping is a mapping whose values are known only at the mesh points of the grid.
     * To interpolate the mapping, we use an interpolation on a basis. The DiscreteMapping is initialised
     * from the coefficients in front of the basis functions which arise when we approximate the
     * functions @f$ E(S) @f$ (with @f$ E @f$ the coordinates in the End vector space and S the coordinates
     * in the Start vector space). Then to interpolate the mapping, we will evaluate the decomposed
     * functions on the chosen interpolating function (see DiscreteMapping::operator()).
     *
     * Here, the evaluator is given as input.
     *
     * @param[in] coeff_representation
     * 		The coefficients of the interpolating function representing the coords in the
     *      target domain.
     *
     * @param[in] evaluator
     * 		The evaluator used to evaluate the mapping.
     */
    DiscreteMapping(CoeffField coeff_representation, SplineEvaluator const& evaluator)
        : m_coeff_representation(coeff_representation)
        , m_spline_evaluator(evaluator)
    {
    }

    /**
     * @brief Compute the target coordinates from the start coordinates.
     *
     * @param[in] coord
     * 			The coordinates in the start domain.
     *
     * @return The coordinates of the mapping in the end domain.
     */
    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        return CoordResult(to_coord(evaluate(m_evaluator, coord, m_coeff_representation)));
    }

    /**
     * @brief Compute the physical coordinates at the points on the logical grid.
     *
     * It evaluates the decomposed mapping on B-splines at the coordinate point
     * with a SplineEvaluator2D.
     *
     * @param[in] exec_space
     *          The execution space where the calculation should be carried out.
     * @param[out] coords
     *          The coordinates of the mapping in the physical domain at the grid
     *          points.
     *
     * @see SplineEvaluator2D
     */
    template <class ExecSpace, class... GridType>
    void operator()(ExecSpace exec_space, Field<CoordResult, IdxRange<GridType...>> coords)
    {
        static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible);
        static_assert(((ddc::in_tags_v<
                        typename GridType::continuous_dimension_type,
                        ddc::to_type_seq_t<StartCoord>>)&&...));
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space,
                get_idx_range(coords),
                KOKKOS_CLASS_LAMBDA(Idx<GridType...> idx) {
                    coords(idx) = (*this)(ddc::coordinate(idx));
                });
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the function jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     */
    KOKKOS_FUNCTION Tensor<DataType, ResultBasis, get_covariant_dims_t<ArgBasis>> jacobian_matrix(
            CoordJacobian const& coord) const
    {
        Tensor<DataType, ResultBasis, get_covariant_dims_t<ArgBasis>> jacobian_matrix;
        details::fill_jacobian_matrix(jacobian_matrix, *this);
        return jacobian_matrix;
    }

    /**
     * @brief Compute the (i,j) coefficient of the Jacobian matrix.
     * 
     * For a mapping given by @f$ \mathcal{F} : {q_i}\mapsto {x_i} @f$,
     *  with @f${q_i}@f$ the curvilinear coordinates and @f${x_i}@f$ the Cartesian coordinates, the
     * (i,j) coefficient of the Jacobian matrix is given by @f$ J^i_j\frac{\partial x_i}{\partial q_j} @f$.
     * 
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the Jacobian matrix.
     *
     * @see SplineEvaluator2D
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(Coord<R, Theta> coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ddc::to_type_seq_t<EndCoord>>);
        static_assert(
                ddc::in_tags_v<IndexTag2, get_covariant_dims_t<ddc::to_type_seq_t<StartCoord>>>);

        return m_evaluator
                .deriv(Idx < ddc::Deriv<typename IndexTag2::Dual>(1),
                       coord,
                       get_const_field(ddc::get<IndexTag1>(m_coeff_representation)));
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(CoordJacobian const& coord) const
    {
        Tensor J = jacobian_matrix(coord);
        return determinant(J);
    }
};


namespace mapping_detail {
template <class StartCoord, class EndCoord, class NDEvaluator, class ExecSpace>
struct MappingAccessibility<ExecSpace, DiscreteMapping<StartCoord, EndCoord, NDEvaluator>>
{
    static constexpr bool value
            = Kokkos::SpaceAccessibility<ExecSpace, typename NDEvaluator::memory_space>::accessible;
};

} // namespace mapping_detail
