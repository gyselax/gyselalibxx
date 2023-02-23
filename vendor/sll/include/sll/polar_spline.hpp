#pragma once

#include <ddc/ddc.hpp>

template <class PolarBSplinesType>
struct PolarSpline
{
public:
    using BSplineR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplineP = typename PolarBSplinesType::BSplinesP_tag;

public:
    ddc::Chunk<double, ddc::DiscreteDomain<BSplineR, BSplineP>> spline_coef;
    ddc::Chunk<double, ddc::DiscreteDomain<PolarBSplinesType>> singular_spline_coef;

public:
    PolarSpline<PolarBSplinesType>(
            ddc::DiscreteDomain<PolarBSplinesType> singular_domain,
            ddc::DiscreteDomain<BSplineR, BSplineP> domain)
        : spline_coef(domain)
        , singular_spline_coef(singular_domain.take_first(
                  ddc::DiscreteVector<PolarBSplinesType>(PolarBSplinesType::n_singular_basis())))
    {
    }
};

template <class PolarBSplinesType>
struct PolarSplineSpan
{
public:
    using BSplineR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplineP = typename PolarBSplinesType::BSplinesP_tag;

public:
    ddc::ChunkSpan<double, ddc::DiscreteDomain<BSplineR, BSplineP>> spline_coef;
    ddc::ChunkSpan<double, ddc::DiscreteDomain<PolarBSplinesType>> singular_spline_coef;

public:
    PolarSplineSpan<PolarBSplinesType>(PolarSpline<PolarBSplinesType>& spl)
        : spline_coef(spl.spline_coef.span_view())
        , singular_spline_coef(spl.singular_spline_coef.span_view())
    {
    }
};

template <class PolarBSplinesType>
struct PolarSplineView
{
public:
    using BSplineR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplineP = typename PolarBSplinesType::BSplinesP_tag;

public:
    ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplineR, BSplineP>> const spline_coef;
    ddc::ChunkSpan<double const, ddc::DiscreteDomain<PolarBSplinesType>> const singular_spline_coef;

public:
    PolarSplineView<PolarBSplinesType>(PolarSpline<PolarBSplinesType> const& spl)
        : spline_coef(spl.spline_coef.span_cview())
        , singular_spline_coef(spl.singular_spline_coef.span_cview())
    {
    }
};
