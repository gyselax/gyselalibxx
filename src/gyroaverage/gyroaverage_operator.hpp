// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "i_interpolation.hpp"


/**
 * @brief Operator to compute the gyroaverage of a field as a weighted sum of the its value at pre-
 * determined sample points along the gyro-motion.
 *
 * Interpolation is used to evaluate the field at points which are not on the grid. The choice of
 * treatment for points outside the domain (Dirichlet, constant outside, etc) is determined by the
 * interpolation class itself.
 *
 * We distinguish 3 types of dimensions:
 *  - "Params" refers to the dimensions which parametrise the gyroaverage of a given field; this
 *    usually corresponds to species and mu.
 *  - "InterpSpace" refers to the space dimensions along which the interpolation is performed, e.g.
 *    poloidal cross-sections when the approximation B // e_phi is assumed.
 *  - "BatchSpace" refers to the remaining dimensions, if any. The same coefficients are used
 *    for every batch index.
 *
 * @tparam InterpolatorType
 *     Type of the interpolator.
 * @tparam IdxRangeParams
 *     The additional dimensions parametrising the gyroaverage (usually species and mu).
 * @tparam IdxRangeBatch
 *     The batch dimensions on which the gyroaverage does not depend (default: none).
 */
template <
        concepts::Interpolation InterpolatorType,
        class IdxRangeParams,
        class IdxRangeBatch = IdxRange<>>
class GyroAverageOperator
{
    using InterpSpaceBuilderType = typename InterpolatorType::BuilderType;
    using InterpSpaceEvaluatorType = typename InterpolatorType::EvaluatorType;
    using ExecutionSpace = typename InterpSpaceBuilderType::exec_space;

    using IdxRangeInterpSpace
            = InterpolationBuilderTraits<InterpSpaceBuilderType>::interpolation_idx_range_type;
    using IdxInterpSpace = typename IdxRangeInterpSpace::discrete_element_type;
    using CoordsInterpSpace = ddc::coordinate_of_t<IdxInterpSpace>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxRangeInterpBatchedSpace = ddc::combine_t<IdxRangeBatch, IdxRangeInterpSpace>;
    using IdxInterpBatchedSpace = typename IdxRangeInterpBatchedSpace::discrete_element_type;
    using IdxParams = typename IdxRangeParams::discrete_element_type;
    using IdxRangeGyroSpace = ddc::combine_t<IdxRangeParams, IdxRangeInterpSpace>;
    using IdxGyroSpace = typename IdxRangeGyroSpace::discrete_element_type;
    using IdxRangeGyroBatchedSpace = ddc::combine_t<IdxRangeParams, IdxRangeInterpBatchedSpace>;

    /// Interpolation coefficients
    using IdxRangeInterpCoeffs
            = InterpolationBuilderTraits<InterpSpaceBuilderType>::coeff_idx_range_type;
    using IdxRangeBatchedInterpCoeffs =
            typename InterpolationBuilderTraits<InterpSpaceBuilderType>::
                    template batched_basis_idx_range_type<IdxRangeInterpBatchedSpace>;
    using DFieldMemBatchedInterpCoeffs = DFieldMem<IdxRangeBatchedInterpCoeffs>;
    using DFieldBatchedInterpCoeffs = DField<IdxRangeBatchedInterpCoeffs>;

public:
    /// @brief Dimension indexing the sample points (must be public to comply with nvcc rules).
    struct CoeffList;

private:
    /// Types for the arrays storing the coefficients
    using IdxRangeList = IdxRange<CoeffList>;
    using IdxList = typename IdxRangeList::discrete_element_type;
    using IdxRangeGyroCoeffs = ddc::combine_t<IdxRangeGyroSpace, IdxRangeList>;

public:
    /// @brief A FieldMem of point coordinates for the sample points.
    using GyroCoeffsCoordsMem = FieldMem<CoordsInterpSpace, IdxRangeGyroCoeffs>;
    /// @brief A FieldMem of doubles containing the sample points' weights.
    using GyroCoeffsWeightsMem = DFieldMem<IdxRangeGyroCoeffs>;

    /// @brief A Field of point coordinates for the sample points.
    using GyroCoeffsCoords = Field<CoordsInterpSpace, IdxRangeGyroCoeffs>;
    /// @brief A Field of doubles containing the sample points' weights.
    using GyroCoeffsWeights = DField<IdxRangeGyroCoeffs>;

private:
    /**
     * @brief The builder for the interpolation
     */
    InterpSpaceBuilderType const& m_builder;

    /**
     * @brief The evaluator for the interpolation
     */
    InterpSpaceEvaluatorType const& m_evaluator;

    /**
     * @brief Coordinates and weights of the sample points
     */
    GyroCoeffsCoords const m_gyrocoeffs_coords;
    GyroCoeffsWeights const m_gyrocoeffs_weights;

    /**
     * @brief Number of points (size of the coefficients in the `CoeffList` dimension)
     */
    size_t const m_nb_gyro_points;

public:
    /**
     * @brief Helper function allocating empty coefficients of the right size and type.
     *
     * @param[in] idx_range_params The index range for the parameter sets.
     * @param[in] idx_range_space The index range for the interpolation space.
     * @param[in] nb_gyro_points The number of sample points on the gyro-trajectory (default: 8).
     *
     * @return A tuple containing fields for the coordinates and their weights.
     */
    static std::tuple<GyroCoeffsCoordsMem, GyroCoeffsWeightsMem> build_empty_coefficients(
            IdxRangeParams idx_range_params,
            IdxRangeInterpSpace idx_range_space,
            size_t const nb_gyro_points = 8)
    {
        IdxRangeGyroSpace idx_range_gyrospace(idx_range_params, idx_range_space);

        IdxRange<CoeffList> idx_range_gyrolist(IdxList(0), IdxStep<CoeffList>(nb_gyro_points));
        IdxRangeGyroCoeffs idx_range_gyrocoeffs(idx_range_gyrospace, idx_range_gyrolist);
        GyroCoeffsCoordsMem gyrocoeffs_coords_alloc(
                "gyrocoeffs_coords (GyroAverageOperator::build_coefficients_circular)",
                idx_range_gyrocoeffs);
        GyroCoeffsWeightsMem gyrocoeffs_weights_alloc(
                "gyrocoeffs_weights (GyroAverageOperator::build_coefficients_circular)",
                idx_range_gyrocoeffs);

        return {std::move(gyrocoeffs_coords_alloc), std::move(gyrocoeffs_weights_alloc)};
    }

    /**
     * @brief Helper function building coefficients for a given set of Larmor radii, in geometries
     * where an inverse mapping to Cartesian coordinates is available.
     *
     * The Larmor radii are of the form `larmor_factor_params * larmor_factor_space`.
     *
     * @param[in] transform_from_cartesian Transformation from Cartesian to mesh coordinates.
     * @param[in] larmor_factor_params Factor of the Larmor radius depending on the parameter sets.
     * @param[in] larmor_factor_space Factor of the Larmor radius depending on space.
     * @param[in] nb_gyro_points Number of sample points on the gyro-trajectory (default: 8).
     *
     * @return A tuple containing the point coordinates and their weights that were computed.
     */
    template <concepts::AnalyticalMapping TransformFromCartesian>
    static std::tuple<GyroCoeffsCoordsMem, GyroCoeffsWeightsMem> build_circular_coefficients(
            TransformFromCartesian transform_from_cartesian,
            DField<IdxRangeParams> const larmor_factor_params,
            DField<IdxRangeInterpSpace> const larmor_factor_space,
            size_t const nb_gyro_points = 8)
    {
        using CoordCartesian = typename TransformFromCartesian::CoordArg;

        IdxRangeParams idx_range_params = get_idx_range(larmor_factor_params);
        IdxRangeInterpSpace idx_range_space = get_idx_range(larmor_factor_space);

        auto [gyrocoeffs_coords_alloc, gyrocoeffs_weights_alloc]
                = build_empty_coefficients(idx_range_params, idx_range_space, nb_gyro_points);
        GyroCoeffsCoords gyrocoeffs_coords = get_field(gyrocoeffs_coords_alloc);
        GyroCoeffsWeights gyrocoeffs_weights = get_field(gyrocoeffs_weights_alloc);

        inverse_mapping_t<TransformFromCartesian> transform_to_cartesian
                = transform_from_cartesian.get_inverse_mapping();

        IdxRangeGyroSpace idx_range_gyrospace(idx_range_params, idx_range_space);

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                Kokkos::DefaultExecutionSpace(),
                idx_range_gyrospace,
                KOKKOS_LAMBDA(IdxGyroSpace const idx_gyrospace) {
                    IdxInterpSpace const idx_space(idx_gyrospace);
                    IdxParams const idx_params(idx_gyrospace);

                    CoordCartesian const gyrocentre
                            = transform_to_cartesian(ddc::coordinate(idx_space));
                    double const rho
                            = larmor_factor_params(idx_params) * larmor_factor_space(idx_space);

                    for (std::size_t igyro = 0; igyro < nb_gyro_points; igyro++) {
                        double const alpha = 2 * M_PI * igyro / static_cast<double>(nb_gyro_points);
                        CoordCartesian const sample_point = gyrocentre
                                                            + CoordCartesian(
                                                                    rho * Kokkos::cos(alpha),
                                                                    rho * Kokkos::sin(alpha));

                        IdxList idx_gyro(igyro);
                        gyrocoeffs_coords(idx_gyrospace, idx_gyro)
                                = transform_from_cartesian(sample_point);
                        gyrocoeffs_weights(idx_gyrospace, idx_gyro) = 1. / nb_gyro_points;
                    }
                });

        return {std::move(gyrocoeffs_coords_alloc), std::move(gyrocoeffs_weights_alloc)};
    }

    /**
     * @brief Constructor for the general case where coefficients are computed by the caller.
     *
     * @param[in] interpolator Interpolator for the evaluation at the sample points.
     * @param[in] gyrocoeffs_coords Coordinates of the sample points.
     * @param[in] gyrocoeffs_weights Weights associated with the sample points.
     */
    explicit GyroAverageOperator(
            InterpolatorType const& interpolator,
            GyroCoeffsCoordsMem& gyrocoeffs_coords,
            GyroCoeffsWeightsMem& gyrocoeffs_weights)
        : m_builder(interpolator.get_builder())
        , m_evaluator(interpolator.get_evaluator())
        , m_gyrocoeffs_coords(gyrocoeffs_coords)
        , m_gyrocoeffs_weights(gyrocoeffs_weights)
        , m_nb_gyro_points(get_idx_range(gyrocoeffs_coords).template extent<CoeffList>())
    {
    }

    /**
     * @brief Compute the gyroaverage of a field which depends both on space and parameter sets.
     *
     * The computation is done in-place for each parameter set.
     *
     * @param[in, out] func Field to process on input and its gyroaverage on output.
     */
    void operator()(DField<IdxRangeGyroBatchedSpace> const& func) const
    {
        IdxRangeParams idx_range_params(get_idx_range(func));

        ddc::host_for_each(idx_range_params, [&](IdxParams const idx_params) {
            (*this)(func[idx_params], idx_params);
        });
    }

    /**
     * @brief For a given parameter set, compute the gyroaverage of a field which only depends on
     * space (and on the batch dimensions).
     *
     * @param[in, out] func Field to process on input and its gyroaverage on output.
     * @param[in] idx_param Parameter set whose coefficients will be used.
     */
    void operator()(DField<IdxRangeInterpBatchedSpace> const& func, IdxParams idx_param) const
    {
        IdxRangeInterpSpace const idx_range_space(get_idx_range(func));
        // Ensure that there are coefficients for the points in the input data
        assert(IdxRangeInterpSpace(get_idx_range(m_gyrocoeffs_weights))
                       .contains(idx_range_space.front()));
        assert(IdxRangeInterpSpace(get_idx_range(m_gyrocoeffs_weights))
                       .contains(idx_range_space.back()));
        // Ensure that there are coefficients for the requested parameter set
        assert(IdxRangeParams(get_idx_range(m_gyrocoeffs_weights)).contains(idx_param));

        // Initialise the interpolator
        DFieldMemBatchedInterpCoeffs scoef_alloc(
                "scoef (GyroAvrageOperator::operator())",
                batched_basis_idx_range(m_builder, get_idx_range(func)));
        DFieldBatchedInterpCoeffs scoef = get_field(scoef_alloc);
        m_builder(scoef, get_const_field(func));
        InterpSpaceEvaluatorType evaluator = m_evaluator;

        std::size_t const nb_gyro_points_proxy = m_nb_gyro_points;
        GyroCoeffsCoords const gyrocoeffs_coords_proxy = m_gyrocoeffs_coords;
        GyroCoeffsWeights const gyrocoeffs_weights_proxy = m_gyrocoeffs_weights;
        IdxList const idx_gyro0(get_idx_range(gyrocoeffs_coords_proxy).front());

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(func),
                KOKKOS_LAMBDA(IdxInterpBatchedSpace const idx_interpbatchedspace) {
                    IdxInterpSpace const idx_interpspace(idx_interpbatchedspace);
                    IdxBatch const idx_batch(idx_interpbatchedspace);
                    DConstField<IdxRangeInterpCoeffs> const scoef_slice
                            = get_const_field(scoef[idx_batch]);
                    // Average over gyro points
                    double gavg = 0.0;
                    for (std::size_t igyro = 0; igyro < nb_gyro_points_proxy; igyro++) {
                        IdxList const idx_gyro = idx_gyro0 + igyro;
                        double const value = evaluator(
                                gyrocoeffs_coords_proxy(idx_param, idx_interpspace, idx_gyro),
                                scoef_slice);
                        gavg += gyrocoeffs_weights_proxy(idx_param, idx_interpspace, idx_gyro)
                                * value;
                    }
                    func(idx_interpbatchedspace) = gavg;
                });
    }
};
