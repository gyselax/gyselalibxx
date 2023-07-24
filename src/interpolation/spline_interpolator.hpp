// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "i_interpolator.hpp"

/**
 * @brief A class for interpolating a function using splines.
 *
 */
template <class DDim, class BSplines, BoundCond BcMin, BoundCond BcMax>
class SplineInterpolator : public IInterpolator<DDim>
{
    using CDim = typename DDim::continuous_dimension_type;

private:
    SplineBuilder<BSplines, DDim, BcMin, BcMax> const& m_builder;

    SplineEvaluator<BSplines> const& m_evaluator;

    mutable ddc::Chunk<double, ddc::DiscreteDomain<BSplines>> m_coefs;

    std::vector<double> m_derivs_min_alloc;

    std::vector<double> m_derivs_max_alloc;

public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    SplineInterpolator(
            SplineBuilder<BSplines, DDim, BcMin, BcMax> const& builder,
            SplineEvaluator<BSplines> const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(builder.spline_domain())
        , m_derivs_min_alloc(BSplines::degree() / 2, 0.)
        , m_derivs_max_alloc(BSplines::degree() / 2, 0.)
    {
    }

    ~SplineInterpolator() override = default;

    ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> const inout_data,
            ddc::ChunkSpan<const ddc::Coordinate<CDim>, ddc::DiscreteDomain<DDim>> const
                    coordinates) const override
    {
        std::optional<CDSpan1D> derivs_min;
        std::optional<CDSpan1D> derivs_max;
        if constexpr (BcMin == BoundCond::HERMITE) {
            derivs_min = CDSpan1D(m_derivs_min_alloc.data(), m_derivs_min_alloc.size());
        }
        if constexpr (BcMax == BoundCond::HERMITE) {
            derivs_max = CDSpan1D(m_derivs_max_alloc.data(), m_derivs_max_alloc.size());
        }
        m_builder(m_coefs, inout_data, derivs_min, derivs_max);
        m_evaluator(inout_data, coordinates, m_coefs);
        return inout_data;
    }
};

/**
 * @brief A class which stores information necessary to create an instance of the SplineInterpolator class.
 *
 * This class allows an instance of the SplineInterpolator class where necessary. This allows the
 * memory allocated in the private members of the SplineInterpolator to be freed when the object is not in use.
 * These objects are: m_coefs, m_derivs_min_alloc, m_derivs_max_alloc.
 */
template <class DDim, class BSplines, BoundCond BcMin, BoundCond BcMax>
class PreallocatableSplineInterpolator : public IPreallocatableInterpolator<DDim>
{
    SplineBuilder<BSplines, DDim, BcMin, BcMax> const& m_builder;

    SplineEvaluator<BSplines> const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolator objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSplineInterpolator(
            SplineBuilder<BSplines, DDim, BcMin, BcMax> const& builder,
            SplineEvaluator<BSplines> const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    ~PreallocatableSplineInterpolator() override = default;

    /**
     * Create an instance of the SplineInterpolator class.
     *
     * @return A unique pointer to an instance of the SplineInterpolator class.
     */
    std::unique_ptr<IInterpolator<DDim>> preallocate() const override
    {
        return std::make_unique<
                SplineInterpolator<DDim, BSplines, BcMin, BcMax>>(m_builder, m_evaluator);
    }
};
