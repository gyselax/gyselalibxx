

# File polar\_spline.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**polar\_splines**](dir_a6779ae02b71d57f488d261458bab1ce.md) **>** [**polar\_spline.hpp**](polar__spline_8hpp.md)

[Go to the documentation of this file](polar__spline_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"

template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct PolarSpline;

template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct ConstPolarSpline;


template <class PolarBSplinesType, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct PolarSplineMem
{
public:
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    DFieldMem<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> spline_coef;

    DFieldMem<IdxRange<PolarBSplinesType>, MemSpace> singular_spline_coef;

public:
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

    PolarSplineMem<PolarBSplinesType, MemSpace>(
            IdxRange<PolarBSplinesType> singular_domain,
            IdxRange<BSplinesR, BSplinesTheta> domain)
        : spline_coef(domain)
        , singular_spline_coef(singular_domain.take_first(
                  IdxStep<PolarBSplinesType>(PolarBSplinesType::n_singular_basis())))
    {
    }

    PolarSpline<PolarBSplinesType, MemSpace> span_view()
    {
        return PolarSpline<PolarBSplinesType, MemSpace>(*this);
    }

    ConstPolarSpline<PolarBSplinesType, MemSpace> span_cview() const
    {
        return ConstPolarSpline<PolarBSplinesType, MemSpace>(*this);
    }
};

template <class PolarBSplinesType, class MemSpace>
struct PolarSpline
{
public:
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    DField<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> spline_coef;

    DField<IdxRange<PolarBSplinesType>, MemSpace> singular_spline_coef;

public:
    explicit PolarSpline<PolarBSplinesType, MemSpace>(
            PolarSplineMem<PolarBSplinesType, MemSpace>& spl)
        : spline_coef(get_field(spl.spline_coef))
        , singular_spline_coef(get_field(spl.singular_spline_coef))
    {
    }

    PolarSpline<PolarBSplinesType, MemSpace> span_view()
    {
        return *this;
    }

    ConstPolarSpline<PolarBSplinesType, MemSpace> span_cview() const
    {
        return ConstPolarSpline<PolarBSplinesType, MemSpace>(*this);
    }
};

template <class PolarBSplinesType, class MemSpace>
struct ConstPolarSpline
{
public:
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;

public:
    DConstField<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> const spline_coef;

    DConstField<IdxRange<PolarBSplinesType>, MemSpace> const singular_spline_coef;

public:
    explicit ConstPolarSpline<PolarBSplinesType, MemSpace>(
            PolarSplineMem<PolarBSplinesType, MemSpace> const& spl)
        : spline_coef(get_const_field(spl.spline_coef))
        , singular_spline_coef(get_const_field(spl.singular_spline_coef))
    {
    }

    explicit ConstPolarSpline<PolarBSplinesType, MemSpace>(
            PolarSpline<PolarBSplinesType, MemSpace> const& spl)
        : spline_coef(get_const_field(spl.spline_coef))
        , singular_spline_coef(get_const_field(spl.singular_spline_coef))
    {
    }

    PolarSpline<PolarBSplinesType, MemSpace> span_view() const
    {
        return *this;
    }

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

template <class ExecSpace, class PolarBSplinesType, class MemSpace>
PolarSplineMem<PolarBSplinesType, typename ExecSpace::memory_space> create_mirror(
        ExecSpace const& exec_space,
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    PolarSplineMem<PolarBSplinesType, typename ExecSpace::memory_space>
            dst(get_idx_range(src.singular_spline_coef), get_idx_range(src.spline_coef));
    return dst;
}

template <class PolarBSplinesType, class MemSpace>
PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace> create_mirror(
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace>
            dst(get_idx_range(src.singular_spline_coef), get_idx_range(src.spline_coef));
    return dst;
}

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

template <class PolarBSplinesType, class MemSpace>
PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace> create_mirror_and_copy(
        PolarSpline<PolarBSplinesType, MemSpace> const& src)
{
    PolarSplineMem<PolarBSplinesType, Kokkos::HostSpace> dst = create_mirror(src);
    ddc::parallel_deepcopy(dst.spline_coef, src.spline_coef);
    ddc::parallel_deepcopy(dst.singular_spline_coef, src.singular_spline_coef);
    return dst;
}

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
```


