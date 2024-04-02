#pragma once

#include <ddc/ddc.hpp>

template <class PolarBSplinesType>
struct PolarSplineSpan;

template <class PolarBSplinesType>
struct PolarSplineView;

/**
 * @brief A structure containing the two Chunks necessary to define a spline on a set of
 * polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType>
struct PolarSpline
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplineR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplineP = typename PolarBSplinesType::BSplinesP_tag;

public:
    /**
     * A Chunk containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    ddc::Chunk<double, ddc::DiscreteDomain<BSplineR, BSplineP>> spline_coef;

    /**
     * A Chunk containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    ddc::Chunk<double, ddc::DiscreteDomain<PolarBSplinesType>> singular_spline_coef;

public:
    /**
     * @brief A constructor for the PolarSpline.
     *
     * A constructor for the PolarSpline which takes a 2D domain of bsplines. This domain is then
     * split into the domains which are relevant for a polar bspline constructed from these 2D bsplines.
     * These new domains are then used to create the chunks which define the spline.
     *
     * @param domain A 2D domain of bsplines.
     */
    PolarSpline<PolarBSplinesType>(ddc::DiscreteDomain<BSplineR, BSplineP> domain)
        : spline_coef(ddc::DiscreteDomain<BSplineR, BSplineP>(
                ddc::select<BSplineR>(domain).remove_first(
                        ddc::DiscreteVector<BSplineR>(PolarBSplinesType::continuity + 1)),
                ddc::select<BSplineP>(domain)))
        , singular_spline_coef(ddc::DiscreteDomain<PolarBSplinesType>(
                  ddc::DiscreteElement<PolarBSplinesType>(0),
                  ddc::DiscreteVector<PolarBSplinesType>(PolarBSplinesType::n_singular_basis())))
    {
    }

    /**
     * @brief Initialise the Chunks containing the coefficients of the spline representation.
     *
     * @param singular_domain The domain for the coefficients in front of the b-spline elements
     *                        near the singular point. These are the elements which cannot be
     *                        expressed as a tensor product of 1D bsplines.
     * @param domain          The domain for the coefficients in front of the b-spline elements
     *                        which can be expressed as a tensor product of 1D bsplines.
     */
    PolarSpline<PolarBSplinesType>(
            ddc::DiscreteDomain<PolarBSplinesType> singular_domain,
            ddc::DiscreteDomain<BSplineR, BSplineP> domain)
        : spline_coef(domain)
        , singular_spline_coef(singular_domain.take_first(
                  ddc::DiscreteVector<PolarBSplinesType>(PolarBSplinesType::n_singular_basis())))
    {
    }

    /**
     * Get a modifiable reference to this polar spline.
     *
     * @return A modifiable reference to this polar spline.
     */
    PolarSplineSpan<PolarBSplinesType> span_view()
    {
        return PolarSplineSpan<PolarBSplinesType>(*this);
    }

    /**
     * Get a constant reference to this polar spline view.
     *
     * @return A constant reference to this polar spline.
     */
    PolarSplineView<PolarBSplinesType> span_cview()
    {
        return PolarSplineView<PolarBSplinesType>(*this);
    }
};

/**
 * @brief A structure containing the two ChunkSpans necessary to define a reference to a spline
 * on a set of polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType>
struct PolarSplineSpan
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplineR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplineP = typename PolarBSplinesType::BSplinesP_tag;

public:
    /**
     * A ChunkSpan containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    ddc::ChunkSpan<double, ddc::DiscreteDomain<BSplineR, BSplineP>> spline_coef;

    /**
     * A ChunkSpan containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    ddc::ChunkSpan<double, ddc::DiscreteDomain<PolarBSplinesType>> singular_spline_coef;

public:
    /**
     * Construct a reference to a PolarSpline.
     *
     * @param spl The PolarSpline being referenced.
     */
    PolarSplineSpan<PolarBSplinesType>(PolarSpline<PolarBSplinesType>& spl)
        : spline_coef(spl.spline_coef.span_view())
        , singular_spline_coef(spl.singular_spline_coef.span_view())
    {
    }

    /**
     * Get a modifiable reference to the polar spline referenced by this polar spline view.
     *
     * @return A modifiable reference to a polar spline.
     */
    PolarSplineSpan<PolarBSplinesType> span_view()
    {
        return *this;
    }

    /**
     * Get a constant reference to the polar spline referenced by this polar spline view.
     *
     * @return A constant reference to a polar spline.
     */
    PolarSplineView<PolarBSplinesType> span_cview()
    {
        return PolarSplineView<PolarBSplinesType>(*this);
    }
};

/**
 * @brief A structure containing the two ChunkViews necessary to define a constant reference to a
 * spline on a set of polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType>
struct PolarSplineView
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplineR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplineP = typename PolarBSplinesType::BSplinesP_tag;

public:
    /**
     * A ChunkView containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplineR, BSplineP>> const spline_coef;

    /**
     * A ChunkView containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    ddc::ChunkSpan<double const, ddc::DiscreteDomain<PolarBSplinesType>> const singular_spline_coef;

public:
    /**
     * Construct a constant reference to a PolarSpline.
     *
     * @param spl The PolarSpline being referenced.
     */
    PolarSplineView<PolarBSplinesType>(PolarSpline<PolarBSplinesType> const& spl)
        : spline_coef(spl.spline_coef.span_cview())
        , singular_spline_coef(spl.singular_spline_coef.span_cview())
    {
    }

    /**
     * Construct a constant reference to a PolarSpline from a PolarSplineSpan
     *
     * @param spl The PolarSpline being referenced.
     */
    PolarSplineView<PolarBSplinesType>(PolarSplineSpan<PolarBSplinesType> const& spl)
        : spline_coef(spl.spline_coef.span_cview())
        , singular_spline_coef(spl.singular_spline_coef.span_cview())
    {
    }

    /**
     * Get a reference to the polar spline referenced by this polar spline view.
     *
     * @return A reference to a polar spline.
     */
    PolarSplineSpan<PolarBSplinesType> span_view()
    {
        return *this;
    }

    /**
     * Get a constant reference to the polar spline referenced by this polar spline view.
     *
     * @return A constant reference to a polar spline.
     */
    PolarSplineView<PolarBSplinesType> span_cview()
    {
        return *this;
    }
};
