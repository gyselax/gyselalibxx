

# File i\_interpolation\_builder.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**i\_interpolation\_builder.hpp**](i__interpolation__builder_8hpp.md)

[Go to the documentation of this file](i__interpolation__builder_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <optional>
#include <type_traits>

#include "ddc_aliases.hpp"

template <class Builder>
struct InterpolationBuilderTraits
{
    using data_type = typename Builder::data_type;

    using interpolation_idx_range_type = typename Builder::interpolation_idx_range_type;

    using coeff_idx_range_type = typename Builder::coeff_idx_range_type;

    static constexpr std::size_t rank()
    {
        return interpolation_idx_range_type::rank();
    }

    template <class IdxRangeBatchedInterpolation>
    using batched_basis_idx_range_type =
            typename Builder::template batched_basis_idx_range_type<IdxRangeBatchedInterpolation>;

    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type =
            typename Builder::template batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>;
};

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationDDim,
        ddc::BoundCond BcLower,
        ddc::BoundCond BcUpper,
        ddc::SplineSolver Solver>
struct InterpolationBuilderTraits<ddc::SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        InterpolationDDim,
        BcLower,
        BcUpper,
        Solver>>
{
private:
    using Builder = ddc::SplineBuilder<
            ExecSpace,
            MemorySpace,
            BSplines,
            InterpolationDDim,
            BcLower,
            BcUpper,
            Solver>;

public:
    using data_type = double;

    using interpolation_grid_type = typename Builder::interpolation_discrete_dimension_type;

    using interpolation_idx_range_type = typename Builder::interpolation_domain_type;

    using coeff_idx_range_type = IdxRange<typename Builder::bsplines_type>;

    static constexpr std::size_t rank()
    {
        return 1;
    }

    using basis_domain_type = typename Builder::bsplines_type;

    template <class IdxRangeBatchedInterpolation>
    using batched_basis_idx_range_type =
            typename Builder::template batched_spline_domain_type<IdxRangeBatchedInterpolation>;

    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type =
            typename Builder::template batched_derivs_domain_type<IdxRangeBatchedInterpolation>;
};

template <class Builder, class IdxRangeBatchedInterpolation>
auto batched_basis_idx_range(
        Builder const& builder,
        IdxRangeBatchedInterpolation const& batched_interpolation_domain)
{
    if constexpr (requires { builder.batched_spline_domain(batched_interpolation_domain); }) {
        return builder.batched_spline_domain(batched_interpolation_domain);
    } else {
        return builder.batched_basis_idx_range(batched_interpolation_domain);
    }
}

namespace concepts {

template <class Builder>
concept InterpolationBuilder = requires
{
    typename Builder::exec_space;
    typename Builder::memory_space;
    typename InterpolationBuilderTraits<Builder>::data_type;
    typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type;
}
&&requires(
        Builder const& b,
        Field<typename InterpolationBuilderTraits<Builder>::data_type,
              typename InterpolationBuilderTraits<Builder>::template batched_basis_idx_range_type<
                      typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type>,
              typename Builder::memory_space> coeffs,
        ConstField<
                typename InterpolationBuilderTraits<Builder>::data_type,
                typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type,
                typename Builder::memory_space> vals)
{
    {b(coeffs, vals)};
};

template <class Builder>
concept InterpolationBuilder1D
        = InterpolationBuilder<Builder> &&(InterpolationBuilderTraits<Builder>::rank() == 1)
          && requires(
                  Builder const& b,
                  typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type domain,
                  Field<typename InterpolationBuilderTraits<Builder>::data_type,
                        typename InterpolationBuilderTraits<Builder>::
                                template batched_basis_idx_range_type<
                                        typename InterpolationBuilderTraits<
                                                Builder>::interpolation_idx_range_type>,
                        typename Builder::memory_space> coeffs,
                  ConstField<
                          typename InterpolationBuilderTraits<Builder>::data_type,
                          typename InterpolationBuilderTraits<
                                  Builder>::interpolation_idx_range_type,
                          typename Builder::memory_space> vals,
                  std::optional<ConstField<
                          typename InterpolationBuilderTraits<Builder>::data_type,
                          typename InterpolationBuilderTraits<Builder>::
                                  template batched_derivs_idx_range_type<
                                          typename InterpolationBuilderTraits<
                                                  Builder>::interpolation_idx_range_type>,
                          typename Builder::memory_space>> derivs)
{
    {b(coeffs, vals, derivs, derivs)};
    {
        b.batched_derivs_xmin_domain(domain)
        } -> std::same_as<
                typename InterpolationBuilderTraits<Builder>::
                        template batched_derivs_idx_range_type<typename InterpolationBuilderTraits<
                                Builder>::interpolation_idx_range_type>>;
    {
        b.batched_derivs_xmax_domain(domain)
        } -> std::same_as<
                typename InterpolationBuilderTraits<Builder>::
                        template batched_derivs_idx_range_type<typename InterpolationBuilderTraits<
                                Builder>::interpolation_idx_range_type>>;
};

} // namespace concepts

template <concepts::InterpolationBuilder1D BuilderType>
using interpolation_grid_t = ddc::type_seq_element_t<
        0,
        ddc::to_type_seq_t<
                typename InterpolationBuilderTraits<BuilderType>::interpolation_idx_range_type>>;
```


