// SPDX-License-Identifier: MIT

#pragma once

#include "i_interpolator_batched.hpp"

// TODO: create splines.hpp
#include <ddc/kernels/splines/spline_builder.hpp>
#include <ddc/kernels/splines/spline_builder_batched.hpp>
#include <ddc/kernels/splines/spline_evaluator.hpp>
#include <ddc/kernels/splines/spline_evaluator_batched.hpp>

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
template <class DDimI, class BSplines, ddc::BoundCond BcMin, ddc::BoundCond BcMax, class... DDim>
class SplineInterpolatorBatched : public IInterpolatorBatched<DDimI, DDim...>
{
    using BuilderType = ddc::SplineBuilderBatched<
            ddc::SplineBuilder<
                    Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultExecutionSpace::memory_space,
                    BSplines,
                    DDimI,
                    BcMin,
                    BcMax>,
            DDim...>;
    using EvaluatorType = ddc::SplineEvaluatorBatched<
            ddc::SplineEvaluator<
                    Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultExecutionSpace::memory_space,
                    BSplines,
                    DDimI>,
            DDim...>;

private:
    BuilderType const& m_builder;

    EvaluatorType const& m_evaluator;

    mutable ddc::
            Chunk<double, typename BuilderType::spline_domain_type, ddc::DeviceAllocator<double>>
                    m_coefs;

    std::vector<double> m_derivs_min_alloc;

    std::vector<double> m_derivs_max_alloc;

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
        , m_derivs_min_alloc(BuilderType::builder_type::s_nbe_xmin, 0.)
        , m_derivs_max_alloc(BuilderType::builder_type::s_nbe_xmax, 0.)
    {
    }

    ~SplineInterpolatorBatched() override = default;

    device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> operator()(
            device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> const inout_data,
            device_t<ddc::ChunkSpan<
                    const ddc::Coordinate<typename DDimI::continuous_dimension_type>,
                    ddc::DiscreteDomain<DDim...>>> const coordinates) const override
    {
        std::optional<ddc::CDSpan1D> derivs_min;
        std::optional<ddc::CDSpan1D> derivs_max;
        if constexpr (BcMin == ddc::BoundCond::HERMITE) {
            derivs_min = ddc::CDSpan1D(m_derivs_min_alloc.data(), m_derivs_min_alloc.size());
        }
        if constexpr (BcMax == ddc::BoundCond::HERMITE) {
            derivs_max = ddc::CDSpan1D(m_derivs_max_alloc.data(), m_derivs_max_alloc.size());
        }
        // m_builder(m_coefs.span_view(), inout_data, derivs_min, derivs_max);
        m_builder(m_coefs.span_view(), inout_data);
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
template <class DDimI, class BSplines, ddc::BoundCond BcMin, ddc::BoundCond BcMax, class... DDim>
class PreallocatableSplineInterpolatorBatched
    : public IPreallocatableInterpolatorBatched<DDimI, DDim...>
{
    using BuilderType = ddc::SplineBuilderBatched<
            ddc::SplineBuilder<
                    Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultExecutionSpace::memory_space,
                    BSplines,
                    DDimI,
                    BcMin,
                    BcMax>,
            DDim...>;
    using EvaluatorType = ddc::SplineEvaluatorBatched<
            ddc::SplineEvaluator<
                    Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultExecutionSpace::memory_space,
                    BSplines,
                    DDimI>,
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
                DDim...>>(m_builder, m_evaluator);
    }
};
