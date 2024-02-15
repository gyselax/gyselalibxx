// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/kernels/splines.hpp>

#include "i_interpolator_batched.hpp"

/**
 * @brief A class for interpolating a function using splines.
 *
 * The class is parametrised by multiple template parameters. Please note that CTAD will deduce all these
 * template parameters from the Builder and Evaluator passed as constructor arguments.
 *
 * @tparam DDimI The dimension of interest.
 * @tparam BSplines The BSplines along the dimension of interest.
 * @tparam BcMin The boundary condition at the lower boundary.
 * @tparam BcMax The boundary condition at the upper boundary.
 * @tparam DDim... All the dimensions of the interpolation problem (batched + interpolated).
 */
template <
        class DDimI,
        class BSplines,
        ddc::BoundCond BcMin,
        ddc::BoundCond BcMax,
        class LeftExtrapolationRule,
        class RightExtrapolationRule,
        ddc::SplineSolver Solver,
        class... DDim>
class SplineInterpolatorBatched : public IInterpolatorBatched<DDimI, DDim...>
{
    using BuilderType = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            DDimI,
            BcMin,
            BcMax,
            Solver,
            DDim...>;
    using EvaluatorType = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            DDimI,
            LeftExtrapolationRule,
            RightExtrapolationRule,
            DDim...>;

private:
    BuilderType const& m_builder;

    EvaluatorType const& m_evaluator;

    mutable ddc::
            Chunk<double, typename BuilderType::spline_domain_type, ddc::DeviceAllocator<double>>
                    m_coefs;

    ddc::Chunk<double, typename BuilderType::derivs_domain_type, ddc::DeviceAllocator<double>>
            m_derivs_min_alloc;
    ddc::Chunk<double, typename BuilderType::derivs_domain_type, ddc::DeviceAllocator<double>>
            m_derivs_max_alloc;

public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    SplineInterpolatorBatched(BuilderType const& builder, EvaluatorType const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(builder.spline_domain())
        , m_derivs_min_alloc(builder.derivs_xmin_domain())
        , m_derivs_max_alloc(builder.derivs_xmax_domain())
    {
        ddc::fill(m_derivs_min_alloc, 0.);
        ddc::fill(m_derivs_max_alloc, 0.);
    }

    ~SplineInterpolatorBatched() override = default;

    device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> operator()(
            device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> const inout_data,
            device_t<ddc::ChunkSpan<
                    const ddc::Coordinate<typename DDimI::continuous_dimension_type>,
                    ddc::DiscreteDomain<DDim...>>> const coordinates) const override
    {
        m_builder(
                m_coefs.span_view(),
                inout_data.span_cview(),
                std::optional(m_derivs_min_alloc.span_cview()),
                std::optional(m_derivs_max_alloc.span_cview()));
        m_evaluator(inout_data, coordinates, m_coefs.span_cview());
        return inout_data;
    }
};

/**
 * @brief A class which stores information necessary to create an instance of the SplineInterpolatorBatched class.
 *
 * This class allows an instance of the SplineInterpolatorBatched class where necessary. This allows the
 * memory allocated in the private members of the SplineInterpolatorBatched to be freed when the object is not in use.
 * These objects are: m_coefs, m_derivs_min_alloc, m_derivs_max_alloc.
 */
template <
        class DDimI,
        class BSplines,
        ddc::BoundCond BcMin,
        ddc::BoundCond BcMax,
        class LeftExtrapolationRule,
        class RightExtrapolationRule,
        ddc::SplineSolver Solver,
        class... DDim>
class PreallocatableSplineInterpolatorBatched
    : public IPreallocatableInterpolatorBatched<DDimI, DDim...>
{
    using BuilderType = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            DDimI,
            BcMin,
            BcMax,
            Solver,
            DDim...>;
    using EvaluatorType = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            DDimI,
            LeftExtrapolationRule,
            RightExtrapolationRule,
            DDim...>;

    BuilderType const& m_builder;

    EvaluatorType const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolatorBatched objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSplineInterpolatorBatched(
            BuilderType const& builder,
            EvaluatorType const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    ~PreallocatableSplineInterpolatorBatched() override = default;

    /**
     * Create an instance of the SplineInterpolatorBatched class.
     *
     * @return A unique pointer to an instance of the SplineInterpolatorBatched class.
     */
    std::unique_ptr<IInterpolatorBatched<DDimI, DDim...>> preallocate() const override
    {
        return std::make_unique<SplineInterpolatorBatched<
                DDimI,
                BSplines,
                BcMin,
                BcMax,
                LeftExtrapolationRule,
                RightExtrapolationRule,
                Solver,
                DDim...>>(m_builder, m_evaluator);
    }
};
