#pragma once

#include <ddc/ddc.hpp>

template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct PolarSplineSpan;

template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct PolarSplineView;


/**
 * @brief A structure containing the two Chunks necessary to define a spline on a set of
 * polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct PolarSpline
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    /**
     * A Chunk containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    ddc::Chunk<
            double,
            ddc::DiscreteDomain<BSplinesR, BSplinesTheta>,
            ddc::KokkosAllocator<double, MemSpace>>
            spline_coef;

    /**
     * A Chunk containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    ddc::Chunk<
            double,
            ddc::DiscreteDomain<PolarBSplinesType>,
            ddc::KokkosAllocator<double, MemSpace>>
            singular_spline_coef;

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
    PolarSpline<PolarBSplinesType, MemSpace>(ddc::DiscreteDomain<BSplinesR, BSplinesTheta> domain)
        : spline_coef(ddc::DiscreteDomain<BSplinesR, BSplinesTheta>(
                ddc::select<BSplinesR>(domain).remove_first(
                        ddc::DiscreteVector<BSplinesR>(PolarBSplinesType::continuity + 1)),
                ddc::select<BSplinesTheta>(domain)))
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
    PolarSpline<PolarBSplinesType, MemSpace>(
            ddc::DiscreteDomain<PolarBSplinesType> singular_domain,
            ddc::DiscreteDomain<BSplinesR, BSplinesTheta> domain)
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
    PolarSplineSpan<PolarBSplinesType, MemSpace> span_view()
    {
        return PolarSplineSpan<PolarBSplinesType, MemSpace>(*this);
    }

    /**
     * Get a constant reference to this polar spline view.
     *
     * @return A constant reference to this polar spline.
     */
    PolarSplineView<PolarBSplinesType, MemSpace> span_cview() const
    {
        return PolarSplineView<PolarBSplinesType, MemSpace>(*this);
    }
};

/**
 * @brief A structure containing the two ChunkSpans necessary to define a reference to a spline
 * on a set of polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType, class MemSpace>
struct PolarSplineSpan
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    /**
     * A ChunkSpan containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    ddc::ChunkSpan<
            double,
            ddc::DiscreteDomain<BSplinesR, BSplinesTheta>,
            Kokkos::layout_right,
            MemSpace>
            spline_coef;

    /**
     * A ChunkSpan containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    ddc::ChunkSpan<double, ddc::DiscreteDomain<PolarBSplinesType>, Kokkos::layout_right, MemSpace>
            singular_spline_coef;

public:
    /**
     * Construct a reference to a PolarSpline.
     *
     * @param spl The PolarSpline being referenced.
     */
    PolarSplineSpan<PolarBSplinesType, MemSpace>(PolarSpline<PolarBSplinesType, MemSpace>& spl)
        : spline_coef(spl.spline_coef.span_view())
        , singular_spline_coef(spl.singular_spline_coef.span_view())
    {
    }

    /**
     * Get a modifiable reference to the polar spline referenced by this polar spline view.
     *
     * @return A modifiable reference to a polar spline.
     */
    PolarSplineSpan<PolarBSplinesType, MemSpace> span_view()
    {
        return *this;
    }

    /**
     * Get a constant reference to the polar spline referenced by this polar spline view.
     *
     * @return A constant reference to a polar spline.
     */
    PolarSplineView<PolarBSplinesType, MemSpace> span_cview() const
    {
        return PolarSplineView<PolarBSplinesType, MemSpace>(*this);
    }
};

/**
 * @brief A structure containing the two ChunkViews necessary to define a constant reference to a
 * spline on a set of polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType, class MemSpace>
struct PolarSplineView
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    /**
     * A ChunkView containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    ddc::ChunkSpan<
            double const,
            ddc::DiscreteDomain<BSplinesR, BSplinesTheta>,
            Kokkos::layout_right,
            MemSpace> const spline_coef;

    /**
     * A ChunkView containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    ddc::ChunkSpan<
            double const,
            ddc::DiscreteDomain<PolarBSplinesType>,
            Kokkos::layout_right,
            MemSpace> const singular_spline_coef;

public:
    /**
     * Construct a constant reference to a PolarSpline.
     *
     * @param spl The PolarSpline being referenced.
     */
    PolarSplineView<PolarBSplinesType, MemSpace>(
            PolarSpline<PolarBSplinesType, MemSpace> const& spl)
        : spline_coef(spl.spline_coef.span_cview())
        , singular_spline_coef(spl.singular_spline_coef.span_cview())
    {
    }

    /**
     * Construct a constant reference to a PolarSpline from a PolarSplineSpan
     *
     * @param spl The PolarSpline being referenced.
     */
    PolarSplineView<PolarBSplinesType, MemSpace>(
            PolarSplineSpan<PolarBSplinesType, MemSpace> const& spl)
        : spline_coef(spl.spline_coef.span_cview())
        , singular_spline_coef(spl.singular_spline_coef.span_cview())
    {
    }

    /**
     * Get a reference to the polar spline referenced by this polar spline view.
     *
     * @return A reference to a polar spline.
     */
    PolarSplineSpan<PolarBSplinesType, MemSpace> span_view() const
    {
        return *this;
    }

    /**
     * Get a constant reference to the polar spline referenced by this polar spline view.
     *
     * @return A constant reference to a polar spline.
     */
    PolarSplineView<PolarBSplinesType, MemSpace> span_cview() const
    {
        return *this;
    }
};

template <class T>
inline constexpr bool is_polar_spline_v = false;

template <class PolarBSplinesType, class MemSpace>
inline constexpr bool is_polar_spline_v<PolarSpline<PolarBSplinesType, MemSpace>> = true;

template <class PolarBSplinesType, class MemSpace>
inline constexpr bool is_polar_spline_v<PolarSplineSpan<PolarBSplinesType, MemSpace>> = true;

template <class PolarBSplinesType, class MemSpace>
inline constexpr bool is_polar_spline_v<PolarSplineView<PolarBSplinesType, MemSpace>> = true;

/**
* @brief A function to create a PolarSpline instance, which has the same attributes as src,
* except the memory space which is related to ExecSpace.
*
* @param[in] exec_space Execution space from which the result must be accessible.
* @param[in] src A reference to PolarSpline.
*/
template <class ExecSpace, class PolarBSplinesType, class MemSpace>
PolarSpline<PolarBSplinesType, typename ExecSpace::memory_space> create_mirror(
        ExecSpace const& exec_space,
        PolarSplineSpan<PolarBSplinesType, MemSpace> const& src)
{
    PolarSpline<PolarBSplinesType, typename ExecSpace::memory_space>
            dst(src.singular_spline_coef.domain(), src.spline_coef.domain());
    return dst;
}

/**
* @brief A function to create a PolarSpline instance, which has the same attributes as src,
* This function allocates memory on the host.
* @param[in] src A reference to PolarSpline.
*/
template <class PolarBSplinesType, class MemSpace>
PolarSpline<PolarBSplinesType, Kokkos::HostSpace> create_mirror(
        PolarSplineSpan<PolarBSplinesType, MemSpace> const& src)
{
    return create_mirror(Kokkos::DefaultHostExecutionSpace(), src);
}

/**
* @brief A function to create copies of PolarSpline class instances. It first creates an instance
* on a memory space accessible from the specified execution space, then it copies the data to the new instance.
*
* @param[in] exec_space Execution space for allocation.
* @param[in] src A reference to PolarSpline.
*/
template <class ExecSpace, class PolarBSplinesType, class MemSpace>
PolarSpline<PolarBSplinesType, typename ExecSpace::memory_space> create_mirror_and_copy(
        ExecSpace const& exec_space,
        PolarSplineSpan<PolarBSplinesType, MemSpace> const& src)
{
    PolarSpline<PolarBSplinesType, typename ExecSpace::memory_space> dst
            = create_mirror(exec_space, src);
    ddc::parallel_deepcopy(dst.spline_coef, src.spline_coef);
    ddc::parallel_deepcopy(dst.singular_spline_coef, src.singular_spline_coef);
    return dst;
}

/**
* @brief A function to create copies for PolarSpline class.It first creates an instance of the class on host,
* then it operates a copy of the data.
*
* @param[in] src A reference to PolarSpline.
*/
template <class PolarBSplinesType, class MemSpace>
PolarSpline<PolarBSplinesType, Kokkos::HostSpace> create_mirror_and_copy(
        PolarSplineSpan<PolarBSplinesType, MemSpace> const& src)
{
    return create_mirror_and_copy(Kokkos::DefaultHostExecutionSpace(), src);
}

/**
* @brief A function to create a mirror view instance for PolarSpline class on the specified execution space.
*
* @param[in] exec_space Execution space from which the result must be accessible.
* @param[in] src A reference to PolarSpline.
*/
template <class ExecSpace, class PolarBSplinesType, class MemSpace>
auto create_mirror_view(
        ExecSpace const& exec_space,
        PolarSplineSpan<PolarBSplinesType, MemSpace> const& src)
{
    if constexpr (std::is_same_v<MemSpace, typename ExecSpace::memory_space>) {
        return src;
    } else {
        PolarSpline<PolarBSplinesType, typename ExecSpace::memory_space>
                dst(src.singular_spline_coef.domain(), src.spline_coef.domain());
        return dst;
    }
}

/**
* @brief A function to create a host allocated mirror view for PolarSpline class.
*
* @param[in] src A reference to PolarSpline.
*/
template <class PolarBSplinesType, class MemSpace>
auto create_mirror_view(PolarSplineSpan<PolarBSplinesType, MemSpace> const& src)
{
    return create_mirror_view(Kokkos::DefaultHostExecutionSpace(), src);
}

/**
* @brief A function to create copies for PolarSpline class.
* If src is accessible from exec_space, src is returned, else, it first creates a mirror view of the class,
* then it copies the data to the new instance. 
*
* @param[in] exec_space Execution space for allocation.
* @param[in] src A reference to PolarSpline.
*/
template <class ExecSpace, class PolarBSplinesType, class MemSpace>
auto create_mirror_view_and_copy(
        ExecSpace const& exec_space,
        PolarSplineSpan<PolarBSplinesType, MemSpace> const& src)
{
    auto dst = create_mirror_view(exec_space, src);
    ddc::parallel_deepcopy(dst.spline_coef, src.spline_coef);
    ddc::parallel_deepcopy(dst.singular_spline_coef, src.singular_spline_coef);
    return dst;
}

/**
* @brief A function to create host allocated view for PolarSpline class.
* @see create_mirror_view_and_copy.
*
* @param[in] src A reference to PolarSpline.
*/
template <class PolarBSplinesType, class MemSpace>
auto create_mirror_view_and_copy(PolarSplineSpan<PolarBSplinesType, MemSpace> const& src)
{
    return create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), src);
}