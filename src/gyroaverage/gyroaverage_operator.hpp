// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "deriv_details.hpp"
#include "i_interpolation.hpp"


/**
 * @brief Operator to compute the gyroaverage of a field as a weighted sum of the its value at pre-
 * determined sample points along the gyro-motion.
 *
 * Spline interpolation is used to evaluate the field at points which are not on the grid. The
 * choice of treatment for points outside the domain (Dirichlet, constant outside, etc) is determined
 * by the interpolation class itself.
 *
 * We distinguish 2 types of dimensions:
 *  - "Params" refers to the dimensions which parametrise the gyroaverage of a given field; this is
 * usually corresponds to species and mu.
 *  - "InterpSpace" refers to the space dimensions along which the interpolations is performed, e.g.
 * poloidal cross-sections when the approximation B // e_phi is assumed.
 *  - "BatchSpace" refers to the remaining dimensions, if any (not yet implemented).
 *
 * @tparam SplineInterpSpaceBuilder
 *     Type of the spline builder for the interpolations.
 * @tparam SplineInterpSpaceEvaluator
 *     Type of the spline evaluator for the interpolations.
 * @tparam IdxRangeParams
 *     The additional dimensions parametrising the gyroaverage (usually species and mu).
 */
template <class SplineInterpSpaceBuilder, class SplineInterpSpaceEvaluator, class IdxRangeParams>
class GyroAverageOperator
{
    //static_assert(
    //ddc::is_evaluator_admissible_v<SplineInterpSpaceBuilder, SplineInterpSpaceEvaluator>,
    //"SplineInterpSpaceEvaluator must be admissible to SplineInterpSpaceBuilder");

public:
    using IdxRangeInterpSpace
            = InterpolationBuilderTraits<SplineInterpSpaceBuilder>::interpolation_idx_range_type;
    using IdxInterpSpace = typename IdxRangeInterpSpace::discrete_element_type;
    using CoordsInterpSpace = ddc::coordinate_of_t<IdxInterpSpace>;
    using IdxParams = typename IdxRangeParams::discrete_element_type;
    using IdxRangeGyroSpace = detail::combine_t<IdxRangeParams, IdxRangeInterpSpace>;
    using IdxGyroSpace = typename IdxRangeGyroSpace::discrete_element_type;

    // Splines
    using IdxRangeSplines
            = InterpolationBuilderTraits<SplineInterpSpaceBuilder>::coeff_idx_range_type;
    using DFieldMemSplines = DFieldMem<IdxRangeSplines>;
    using DFieldSplines = DField<IdxRangeSplines>;

    // Types for the arrays storing the coefficients
    struct DimList;
    using IdxRangeList = IdxRange<DimList>;
    using IdxList = typename IdxRangeList::discrete_element_type;
    using IdxRangeGyroCoeffs = detail::combine_t<IdxRangeGyroSpace, IdxRangeList>;

    using GyroCoeffsCoordsMem = FieldMem<CoordsInterpSpace, IdxRangeGyroCoeffs>;
    using GyroCoeffsWeightsMem = FieldMem<double, IdxRangeGyroCoeffs>;
    using GyroCoeffsCoords = Field<CoordsInterpSpace, IdxRangeGyroCoeffs>;
    using GyroCoeffsWeights = Field<double, IdxRangeGyroCoeffs>;


private:
    /**
     * @brief The spline builder for the interpolation
     */
    SplineInterpSpaceBuilder const& m_spline_builder;

    /**
     * @brief The spline evaluator for the interpolation
     */
    SplineInterpSpaceEvaluator const& m_spline_evaluator;

    /**
     * @brief Coordinates and weights of the sample points
     */
    GyroCoeffsCoords const m_gyrocoeffs_coords;
    GyroCoeffsWeights const m_gyrocoeffs_weights;

    /**
     * @brief Number of points (size of the coefficients in the `DimList` dimension)
     */
    size_t const m_nb_gyro_points;

public:
    /**
     * @brief Helper function allocating empty coefficients of the right size and type.
     *
     * @param[in] idx_range_params The index range for the parameter sets.
     * @param[in] idx_range_space The index range for the interpolation space.
     * @param[in] nb_gyro_points The number of sample points in the gyro-trajectory.
     *
     * @return A tuple containing fields for the coordinates and their weights.
     */
    static std::tuple<GyroCoeffsCoordsMem, GyroCoeffsWeightsMem> build_empty_coefficients(
            IdxRangeParams idx_range_params,
            IdxRangeInterpSpace idx_range_space,
            size_t const nb_gyro_points)
    {
        IdxRangeGyroSpace idx_range_gyrospace(idx_range_params, idx_range_space);

        IdxRange<DimList> idx_range_gyrolist(IdxList(0), IdxStep<DimList>(nb_gyro_points));
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
     * @param[in] transform_from_cartesian Transformation from Cartesian to mesh cordinates.
     * @param[in] larmor_factor_params Factor of the Larmor radius depending on the parameter sets.
     * @param[in] larmor_factor_space Factor of the Larmor radius depending on space.
     * @param[in] nb_gyro_points Number of sample points in the gyro-trajectory.
     *
     * @return A tuple containing the point coordinates and their weights that were computed.
     */
    template <class TransformToCartesian>
    static std::tuple<GyroCoeffsCoordsMem, GyroCoeffsWeightsMem> build_circular_coefficients(
            TransformToCartesian transform_from_cartesian,
            DField<IdxRangeParams> const larmor_factor_params,
            DField<IdxRangeInterpSpace> const larmor_factor_space,
            size_t const nb_gyro_points)
    {
        using CoordCartesian = typename TransformToCartesian::CoordArg;

        IdxRangeParams idx_range_params = get_idx_range(larmor_factor_params);
        IdxRangeInterpSpace idx_range_space = get_idx_range(larmor_factor_space);

        auto [gyrocoeffs_coords_alloc, gyrocoeffs_weights_alloc]
                = build_empty_coefficients(idx_range_params, idx_range_space, 8);
        GyroCoeffsCoords gyrocoeffs_coords = get_field(gyrocoeffs_coords_alloc);
        GyroCoeffsWeights gyrocoeffs_weights = get_field(gyrocoeffs_weights_alloc);

        IdxRangeGyroSpace idx_range_gyrospace(idx_range_params, idx_range_space);

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                Kokkos::DefaultExecutionSpace(),
                idx_range_gyrospace,
                KOKKOS_LAMBDA(IdxGyroSpace const idx_gyrospace) {
                    IdxInterpSpace const idx_space(idx_gyrospace);
                    IdxParams const idx_params(idx_gyrospace);

                    inverse_mapping_t<TransformToCartesian> transform_to_cartesian
                            = transform_from_cartesian.get_inverse_mapping();

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
     * @param[in] spline_builder Spline builder for the interpolation.
     * @param[in] spline_evaluator Spline evaluator for the interpolation.
     * @param[in] gyrocoeffs_coords Coordinates of the sample points.
     * @param[in] gyrocoeffs_weights Weights associated with the sample points.
     */
    explicit GyroAverageOperator(
            SplineInterpSpaceBuilder const& spline_builder,
            SplineInterpSpaceEvaluator const& spline_evaluator,
            GyroCoeffsCoordsMem& gyrocoeffs_coords,
            GyroCoeffsWeightsMem& gyrocoeffs_weights)
        : m_spline_builder(spline_builder)
        , m_spline_evaluator(spline_evaluator)
        , m_gyrocoeffs_coords(gyrocoeffs_coords)
        , m_gyrocoeffs_weights(gyrocoeffs_weights)
        , m_nb_gyro_points(get_idx_range(gyrocoeffs_coords).template extent<DimList>())
    {
    }

    /**
     * @brief Compute the gyroaverage of a field which depends both on space and parameter sets.
     *
     * The computation is done in-place for each parameter set.
     *
     * @param[in, out] func Field to process on input and its gyroaverage on output.
     */
    void operator()(DField<IdxRangeGyroSpace> const& func) const
    {
        IdxRangeParams idx_range_params(get_idx_range(func));

        ddc::host_for_each(idx_range_params, [&](IdxParams const idx_params) {
            (*this)(func[idx_params], get_const_field(func[idx_params]), idx_params);
        });
    }

    /**
     * @brief For a given parameter set, compute the gyroaverage of a field which only depends on
     * space.
     *
     * The same field can be passed as both input and output to work in-place, as the spline
     * coefficients represent the input field during the computation.
     *
     * @param[out] func_gavg Gyroaverage of the input field.
     * @param[in] func The field to process.
     * @param[in] idx_param Parameter set whose coefficients will be used.
     */
    void operator()(
            DField<IdxRangeInterpSpace> const& func_gavg,
            DConstField<IdxRangeInterpSpace> const& func,
            IdxParams idx_param) const
    {
        // Ensure that the input and output have the same index ranges
        assert(get_idx_range(func).extents() == get_idx_range(func_gavg).extents());
        // Ensure that there are coefficients for the points in the input data
        assert(IdxRangeInterpSpace(get_idx_range(m_gyrocoeffs_weights))
                       .contains(get_idx_range(func).front()));
        assert(IdxRangeInterpSpace(get_idx_range(m_gyrocoeffs_weights))
                       .contains(get_idx_range(func).back()));
        // Ensure that there are coefficients for the requested parameter set
        assert(IdxRangeParams(get_idx_range(m_gyrocoeffs_weights)).contains(idx_param));

        // Initialise the splines
        DFieldMemSplines scoef_alloc(
                "scoef (GyroAvrageOperator::operator())",
                get_spline_idx_range(m_spline_builder));
        DFieldSplines scoef = get_field(scoef_alloc);
        m_spline_builder(scoef, get_const_field(func));
        SplineInterpSpaceEvaluator spline_evaluator = m_spline_evaluator;

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(func_gavg),
                KOKKOS_CLASS_LAMBDA(IdxInterpSpace const idx_interpspace) {
                    IdxList const idx_gyro0(get_idx_range(m_gyrocoeffs_coords).front());

                    // Average over gyro points
                    double gavg = 0.0;
                    for (std::size_t igyro = 0; igyro < m_nb_gyro_points; igyro++) {
                        IdxList const idx_gyro = idx_gyro0 + igyro;
                        double const value = spline_evaluator(
                                m_gyrocoeffs_coords(idx_param, idx_interpspace, idx_gyro),
                                get_const_field(scoef));
                        gavg += m_gyrocoeffs_weights(idx_param, idx_interpspace, idx_gyro) * value;
                    }
                    func_gavg(idx_interpspace) = gavg;
                });
    }
};
