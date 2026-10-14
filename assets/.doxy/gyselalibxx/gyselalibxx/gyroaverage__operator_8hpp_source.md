

# File gyroaverage\_operator.hpp

[**File List**](files.md) **>** [**gyroaverage**](dir_0488de83d5a156211532c59817863a53.md) **>** [**gyroaverage\_operator.hpp**](gyroaverage__operator_8hpp.md)

[Go to the documentation of this file](gyroaverage__operator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "deriv_details.hpp"
#include "i_interpolation.hpp"


template <concepts::Interpolation InterpolatorType, class IdxRangeParams>
class GyroAverageOperator
{
    using InterpSpaceBuilderType = typename InterpolatorType::BuilderType;
    using InterpSpaceEvaluatorType = typename InterpolatorType::EvaluatorType;
    using ExecutionSpace = typename InterpSpaceBuilderType::exec_space;

    using IdxRangeInterpSpace
            = InterpolationBuilderTraits<InterpSpaceBuilderType>::interpolation_idx_range_type;
    using IdxInterpSpace = typename IdxRangeInterpSpace::discrete_element_type;
    using CoordsInterpSpace = ddc::coordinate_of_t<IdxInterpSpace>;
    using IdxParams = typename IdxRangeParams::discrete_element_type;
    using IdxRangeGyroSpace = detail::combine_t<IdxRangeParams, IdxRangeInterpSpace>;
    using IdxGyroSpace = typename IdxRangeGyroSpace::discrete_element_type;

    using IdxRangeInterpCoeffs
            = InterpolationBuilderTraits<InterpSpaceBuilderType>::coeff_idx_range_type;
    using DFieldMemInterpCoeffs = DFieldMem<IdxRangeInterpCoeffs>;
    using DFieldInterpCoeffs = DField<IdxRangeInterpCoeffs>;

public:
    struct CoeffList;

private:
    using IdxRangeList = IdxRange<CoeffList>;
    using IdxList = typename IdxRangeList::discrete_element_type;
    using IdxRangeGyroCoeffs = detail::combine_t<IdxRangeGyroSpace, IdxRangeList>;

public:
    using GyroCoeffsCoordsMem = FieldMem<CoordsInterpSpace, IdxRangeGyroCoeffs>;
    using GyroCoeffsWeightsMem = DFieldMem<IdxRangeGyroCoeffs>;

    using GyroCoeffsCoords = Field<CoordsInterpSpace, IdxRangeGyroCoeffs>;
    using GyroCoeffsWeights = DField<IdxRangeGyroCoeffs>;

private:
    InterpSpaceBuilderType const& m_builder;

    InterpSpaceEvaluatorType const& m_evaluator;

    GyroCoeffsCoords const m_gyrocoeffs_coords;
    GyroCoeffsWeights const m_gyrocoeffs_weights;

    size_t const m_nb_gyro_points;

public:
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

    void operator()(DField<IdxRangeGyroSpace> const& func) const
    {
        IdxRangeParams idx_range_params(get_idx_range(func));

        ddc::host_for_each(idx_range_params, [&](IdxParams const idx_params) {
            (*this)(func[idx_params], idx_params);
        });
    }

    void operator()(DField<IdxRangeInterpSpace> const& func, IdxParams idx_param) const
    {
        // Ensure that there are coefficients for the points in the input data
        assert(IdxRangeInterpSpace(get_idx_range(m_gyrocoeffs_weights))
                       .contains(get_idx_range(func).front()));
        assert(IdxRangeInterpSpace(get_idx_range(m_gyrocoeffs_weights))
                       .contains(get_idx_range(func).back()));
        // Ensure that there are coefficients for the requested parameter set
        assert(IdxRangeParams(get_idx_range(m_gyrocoeffs_weights)).contains(idx_param));

        // Initialise the interpolator
        DFieldMemInterpCoeffs scoef_alloc(
                "scoef (GyroAvrageOperator::operator())",
                batched_basis_idx_range(m_builder, get_idx_range(func)));
        DFieldInterpCoeffs scoef = get_field(scoef_alloc);
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
                KOKKOS_LAMBDA(IdxInterpSpace const idx_interpspace) {
                    // Average over gyro points
                    double gavg = 0.0;
                    for (std::size_t igyro = 0; igyro < nb_gyro_points_proxy; igyro++) {
                        IdxList const idx_gyro = idx_gyro0 + igyro;
                        double const value = evaluator(
                                gyrocoeffs_coords_proxy(idx_param, idx_interpspace, idx_gyro),
                                get_const_field(scoef));
                        gavg += gyrocoeffs_weights_proxy(idx_param, idx_interpspace, idx_gyro)
                                * value;
                    }
                    func(idx_interpspace) = gavg;
                });
    }
};
```


