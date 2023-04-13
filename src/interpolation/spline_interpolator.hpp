// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "i_interpolator.hpp"

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

template <class DDim, class BSplines, BoundCond BcMin, BoundCond BcMax>
class PreallocatableSplineInterpolator : public IPreallocatableInterpolator<DDim>
{
    using CDim = typename DDim::continuous_dimension_type;

    SplineBuilder<BSplines, DDim, BcMin, BcMax> const& m_builder;

    SplineEvaluator<BSplines> const& m_evaluator;

public:
    PreallocatableSplineInterpolator(
            SplineBuilder<BSplines, DDim, BcMin, BcMax> const& builder,
            SplineEvaluator<BSplines> const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    ~PreallocatableSplineInterpolator() override = default;

    InterpolatorProxy<DDim> preallocate() const override
    {
        return InterpolatorProxy<DDim>(
                std::make_unique<
                        SplineInterpolator<DDim, BSplines, BcMin, BcMax>>(m_builder, m_evaluator));
    }

    ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> const inout_data,
            ddc::ChunkSpan<const ddc::Coordinate<CDim>, ddc::DiscreteDomain<DDim>> const
                    coordinates) const override
    {
        return preallocate()(inout_data, coordinates);
    }
};
