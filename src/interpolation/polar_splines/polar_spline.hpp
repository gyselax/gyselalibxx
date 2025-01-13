// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"

template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct PolarSpline;

template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct ConstPolarSpline;


/**
 * @brief A structure containing the two FieldMems necessary to define a spline on a set of
 * polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct PolarSplineMem
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    /**
     * A FieldMem containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    DFieldMem<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> spline_coef;

    /**
     * A FieldMem containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    DFieldMem<IdxRange<PolarBSplinesType>, MemSpace> singular_spline_coef;

public:
    /**
     * @brief A constructor for the PolarSplineMem.
     *
     * A constructor for the PolarSplineMem which takes a 2D domain of bsplines. This domain is then
     * split into the domains which are relevant for a polar bspline constructed from these 2D bsplines.
     * These new domains are then used to create the chunks which define the spline.
     *
     * @param domain A 2D domain of bsplines.
     */
    explicit PolarSplineMem<PolarBSplinesType, MemSpace>(IdxRange<BSplinesR, BSplinesTheta> domain)
        : spline_coef(IdxRange<BSplinesR, BSplinesTheta>(
                ddc::select<BSplinesR>(domain).remove_first(
                        IdxStep<BSplinesR>(PolarBSplinesType::continuity + 1)),
                ddc::select<BSplinesTheta>(domain)))
        , singular_spline_coef(IdxRange<PolarBSplinesType>(
                  Idx<PolarBSplinesType>(0),
                  IdxStep<PolarBSplinesType>(PolarBSplinesType::n_singular_basis())))
    {
    }

    /**
     * @brief Initialise the FieldMems containing the coefficients of the spline representation.
     *
     * @param singular_domain The domain for the coefficients in front of the b-spline elements
     *                        near the singular point. These are the elements which cannot be
     *                        expressed as a tensor product of 1D bsplines.
     * @param domain          The domain for the coefficients in front of the b-spline elements
     *                        which can be expressed as a tensor product of 1D bsplines.
     */
    PolarSplineMem<PolarBSplinesType, MemSpace>(
            IdxRange<PolarBSplinesType> singular_domain,
            IdxRange<BSplinesR, BSplinesTheta> domain)
        : spline_coef(domain)
        , singular_spline_coef(singular_domain.take_first(
                  IdxStep<PolarBSplinesType>(PolarBSplinesType::n_singular_basis())))
    {
    }

    /**
     * Get a modifiable reference to this polar spline.
     *
     * @return A modifiable reference to this polar spline.
     */
    PolarSpline<PolarBSplinesType, MemSpace> span_view()
    {
        return PolarSpline<PolarBSplinesType, MemSpace>(*this);
    }

    /**
     * Get a constant reference to this polar spline view.
     *
     * @return A constant reference to this polar spline.
     */
    ConstPolarSpline<PolarBSplinesType, MemSpace> span_cview() const
    {
        return ConstPolarSpline<PolarBSplinesType, MemSpace>(*this);
    }
};

/**
 * @brief A structure containing the two Fields necessary to define a reference to a spline
 * on a set of polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType, class MemSpace>
struct PolarSpline
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    /**
     * A Field containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    DField<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> spline_coef;

    /**
     * A Field containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    DField<IdxRange<PolarBSplinesType>, MemSpace> singular_spline_coef;

public:
    /**
     * Construct a reference to a PolarSplineMem.
     *
     * @param spl The PolarSplineMem being referenced.
     */
    explicit PolarSpline<PolarBSplinesType, MemSpace>(
            PolarSplineMem<PolarBSplinesType, MemSpace>& spl)
        : spline_coef(get_field(spl.spline_coef))
        , singular_spline_coef(get_field(spl.singular_spline_coef))
    {
    }

    /**
     * Get a modifiable reference to the polar spline referenced by this polar spline view.
     *
     * @return A modifiable reference to a polar spline.
     */
    PolarSpline<PolarBSplinesType, MemSpace> span_view()
    {
        return *this;
    }

    /**
     * Get a constant reference to the polar spline referenced by this polar spline view.
     *
     * @return A constant reference to a polar spline.
     */
    ConstPolarSpline<PolarBSplinesType, MemSpace> span_cview() const
    {
        return ConstPolarSpline<PolarBSplinesType, MemSpace>(*this);
    }
};

/**
 * @brief A structure containing the two ConstFields necessary to define a constant reference to a
 * spline on a set of polar basis splines.
 *
 * @tparam PolarBSplinesType The type of the polar bsplines on which this spline is defined.
 */
template <class PolarBSplinesType, class MemSpace>
struct ConstPolarSpline
{
public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    /**
     * A ConstField containing the coefficients in front of the b-spline elements which can be
     * expressed as a tensor product of 1D bsplines.
     */
    DConstField<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> const spline_coef;

    /**
     * A ConstField containing the coefficients in front of the b-spline elements near the
     * singular point which cannot be expressed as a tensor product of 1D bsplines.
     */
    DConstField<IdxRange<PolarBSplinesType>, MemSpace> const singular_spline_coef;

public:
    /**
     * Construct a constant reference to a PolarSplineMem.
     *
     * @param spl The PolarSplineMem being referenced.
     */
    explicit ConstPolarSpline<PolarBSplinesType, MemSpace>(
            PolarSplineMem<PolarBSplinesType, MemSpace> const& spl)
        : spline_coef(get_const_field(spl.spline_coef))
        , singular_spline_coef(get_const_field(spl.singular_spline_coef))
    {
    }

    /**
     * Construct a constant reference to a PolarSplineMem from a PolarSpline
     *
     * @param spl The PolarSplineMem being referenced.
     */
    explicit ConstPolarSpline<PolarBSplinesType, MemSpace>(
            PolarSpline<PolarBSplinesType, MemSpace> const& spl)
        : spline_coef(get_const_field(spl.spline_coef))
        , singular_spline_coef(get_const_field(spl.singular_spline_coef))
    {
    }

    /**
     * Get a reference to the polar spline referenced by this polar spline view.
     *
     * @return A reference to a polar spline.
     */
    PolarSpline<PolarBSplinesType, MemSpace> span_view() const
    {
        return *this;
    }

    /**
     * Get a constant reference to the polar spline referenced by this polar spline view.
     *
     * @return A constant reference to a polar spline.
     */
    ConstPolarSpline<PolarBSplinesType, MemSpace> span_cview() const
    {
        return *this;
    }
};

template <class T>
inline constexpr bool is_polar_spline_v = false;

template <class PolarBSplinesType, class MemSpace>
inline constexpr bool is_polar_spline_v<PolarSplineMem<PolarBSplinesType, MemSpace>> = true;

template <class PolarBSplinesType, class MemSpace>
inline constexpr bool is_polar_spline_v<PolarSpline<PolarBSplinesType, MemSpace>> = true;

template <class PolarBSplinesType, class MemSpace>
inline constexpr bool is_polar_spline_v<ConstPolarSpline<PolarBSplinesType, MemSpace>> = true;

/**
* @brief A function to create a PolarSplineMem instance, which has the same attributes as src,
* except the memory space which is related to ExecSpace.
*
* @param[in] exec_space Execution space from which the result must be accessible.
* @param[in] src A reference to PolarSplineMem.
*/
template <class ExecSpace, class PolarBSplinesType, class MemSpace>
PolarSplineMem<PolarBSplinesType, typename ExecSpace::memory_space> create_mirror(
        ExecSpace const& exec_space,
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    PolarSplineMem<PolarBSplinesType, typename ExecSpace::memory_space>
            dst(get_idx_range(src.singular_spline_coef), get_idx_range(src.spline_coef));
    return dst;
}

/**
* @brief A function to create a PolarSplineMem instance, which has the same attributes as src,
* This function allocates memory on the host.
* @param[in] src A reference to PolarSplineMem.
*/
template <class PolarBSplinesType, class MemSpace>
PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace> create_mirror(
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace>
            dst(get_idx_range(src.singular_spline_coef), get_idx_range(src.spline_coef));
    return dst;
}

/**
* @brief A function to create copies of PolarSplineMem class instances. It first creates an instance
* on a memory space accessible from the specified execution space, then it copies the data to the new instance.
*
* @param[in] exec_space Execution space for allocation.
* @param[in] src A reference to PolarSplineMem.
*/
template <class ExecSpace, class PolarBSplinesType, class MemSpace>
PolarSplineMem<PolarBSplinesType, typename ExecSpace::memory_space> create_mirror_and_copy(
        ExecSpace const& exec_space,
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    PolarSplineMem<PolarBSplinesType, typename ExecSpace::memory_space> dst
            = create_mirror(exec_space, src);
    ddc::parallel_deepcopy(dst.spline_coef, src.spline_coef);
    ddc::parallel_deepcopy(dst.singular_spline_coef, src.singular_spline_coef);
    return dst;
}

/**
* @brief A function to create copies for PolarSplineMem class.It first creates an instance of the class on host,
* then it operates a copy of the data.
*
* @param[in] src A reference to PolarSplineMem.
*/
template <class PolarBSplinesType, class MemSpace>
PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace> create_mirror_and_copy(
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace> dst = create_mirror(src);
    ddc::parallel_deepcopy(dst.spline_coef, src.spline_coef);
    ddc::parallel_deepcopy(dst.singular_spline_coef, src.singular_spline_coef);
    return dst;
}

/**
* @brief A function to create a mirror view instance for PolarSplineMem class on the specified execution space.
*
* @param[in] exec_space Execution space from which the result must be accessible.
* @param[in] src A reference to PolarSplineMem.
*/
template <class ExecSpace, class PolarBSplinesType, class MemSpace>
auto create_mirror_view(
        ExecSpace const& exec_space,
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    if constexpr (std::is_same_v<MemSpace, typename ExecSpace::memory_space>) {
        return src;
    } else {
        PolarSplineMem<PolarBSplinesType, typename ExecSpace::memory_space>
                dst(get_idx_range(src.singular_spline_coef), get_idx_range(src.spline_coef));
        return dst;
    }
}

/**
* @brief A function to create a host allocated mirror view for PolarSplineMem class.
*
* @param[in] src A reference to PolarSplineMem.
*/
template <class PolarBSplinesType, class MemSpace>
auto create_mirror_view(PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    if constexpr (std::is_same_v<MemSpace, Kokkos::HostSpace>) {
        return src;
    } else {
        PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace>
                dst(get_idx_range(src.singular_spline_coef), get_idx_range(src.spline_coef));
        return dst;
    }
}

/**
* @brief A function to create copies for PolarSplineMem class.
* If src is accessible from exec_space, src is returned, else, it first creates a mirror view of the class,
* then it copies the data to the new instance. 
*
* @param[in] exec_space Execution space for allocation.
* @param[in] src A reference to PolarSplineMem.
*/
template <class ExecSpace, class PolarBSplinesType, class MemSpace>
auto create_mirror_view_and_copy(
        ExecSpace const& exec_space,
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    auto dst = create_mirror_view(exec_space, src);
    ddc::parallel_deepcopy(dst.spline_coef, src.spline_coef);
    ddc::parallel_deepcopy(dst.singular_spline_coef, src.singular_spline_coef);
    return dst;
}

/**
* @brief A function to create host allocated view for PolarSplineMem class.
* @see create_mirror_view_and_copy.
*
* @param[in] src A reference to PolarSplineMem.
*/
template <class PolarBSplinesType, class MemSpace>
auto create_mirror_view_and_copy(PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    auto dst = create_mirror_view(src);
    ddc::parallel_deepcopy(dst.spline_coef, src.spline_coef);
    ddc::parallel_deepcopy(dst.singular_spline_coef, src.singular_spline_coef);
    return dst;
}

template <class PolarBSplines, class MemorySpace>
inline constexpr bool enable_data_access_methods<PolarSplineMem<PolarBSplines, MemorySpace>> = true;

template <class PolarBSplines, class MemorySpace>
inline constexpr bool
        enable_data_access_methods<ConstPolarSpline<PolarBSplines, MemorySpace>> = true;

template <class PolarBSplines, class MemorySpace>
inline constexpr bool enable_data_access_methods<PolarSpline<PolarBSplines, MemorySpace>> = true;

template <class PolarBSplines, class MemorySpace>
inline constexpr bool enable_mem_type<PolarSplineMem<PolarBSplines, MemorySpace>> = true;

namespace detail {

/**
 * @brief Get a new `PolarSplineMem` type with the same parametrisation
 * except in the memory space which is set to NewMemorySpace.
 * @tparam NewMemorySpace The new memory space. 
 * @tparam PolarSplineType Type of the B-splines.
 * @tparam MemorySpace The original memory space of the chunk of the VectorFieldMem.
 * @see VectorField
 */
template <class NewMemorySpace, class PolarSplineType, class MemorySpace>
struct OnMemorySpace<NewMemorySpace, PolarSplineMem<PolarSplineType, MemorySpace>>
{
    using type = PolarSplineMem<PolarSplineType, NewMemorySpace>;
};

template <class NewMemorySpace, class PolarSplineType, class MemorySpace>
struct OnMemorySpace<NewMemorySpace, PolarSpline<PolarSplineType, MemorySpace>>
{
    using type = PolarSpline<PolarSplineType, NewMemorySpace>;
};

template <class NewMemorySpace, class PolarSplineType, class MemorySpace>
struct OnMemorySpace<NewMemorySpace, ConstPolarSpline<PolarSplineType, MemorySpace>>
{
    using type = ConstPolarSpline<PolarSplineType, NewMemorySpace>;
};
} // namespace detail
