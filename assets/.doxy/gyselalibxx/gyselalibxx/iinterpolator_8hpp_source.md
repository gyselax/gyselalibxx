

# File iinterpolator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**iinterpolator.hpp**](iinterpolator_8hpp.md)

[Go to the documentation of this file](iinterpolator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

// TODO: Generalise (IDimI -> Tags...) and make it usable for all Gysela operators ?
template <template <class...> class Interp, class GridInterp, class IdxRange>
struct interpolator_on_idx_range
{
};

template <template <class...> class Interp, class GridInterp, class... Grid1D>
struct interpolator_on_idx_range<Interp, GridInterp, IdxRange<Grid1D...>>
{
    using type = Interp<GridInterp, Grid1D...>;
};

template <template <class...> class Interp, class GridInterp, class IdxRange>
using interpolator_on_idx_range_t =
        typename interpolator_on_idx_range<Interp, GridInterp, IdxRange>::type;

template <class GridInterp, class... Grid1D>
class IInterpolator
{
public:
    virtual ~IInterpolator() = default;

    using deriv_type = ddc::Deriv<typename GridInterp::continuous_dimension_type>;
    using batched_derivs_idx_range_type
            = ddc::replace_dim_of_t<IdxRange<Grid1D...>, GridInterp, deriv_type>;

    virtual batched_derivs_idx_range_type batched_derivs_idx_range_xmin(
            IdxRange<Grid1D...> idx_range) const = 0;

    virtual batched_derivs_idx_range_type batched_derivs_idx_range_xmax(
            IdxRange<Grid1D...> idx_range) const = 0;

    virtual Field<double, IdxRange<Grid1D...>> operator()(
            Field<double, IdxRange<Grid1D...>> inout_data,
            Field<const Coord<typename GridInterp::continuous_dimension_type>, IdxRange<Grid1D...>>
                    coordinates,
            std::optional<Field<double const, batched_derivs_idx_range_type>> derivs_xmin
            = std::nullopt,
            std::optional<Field<double const, batched_derivs_idx_range_type>> derivs_xmax
            = std::nullopt) const = 0;
};

template <class GridInterp, class... Grid1D>
class IPreallocatableInterpolator : public IInterpolator<GridInterp, Grid1D...>
{
public:
    ~IPreallocatableInterpolator() override = default;

    using deriv_type = typename IInterpolator<GridInterp, Grid1D...>::deriv_type;
    using batched_derivs_idx_range_type =
            typename IInterpolator<GridInterp, Grid1D...>::batched_derivs_idx_range_type;

    virtual std::unique_ptr<IInterpolator<GridInterp, Grid1D...>> preallocate() const = 0;

    batched_derivs_idx_range_type batched_derivs_idx_range_xmin(
            IdxRange<Grid1D...> idx_range) const override
    {
        return (*preallocate()).batched_derivs_idx_range_xmin(idx_range);
    }

    batched_derivs_idx_range_type batched_derivs_idx_range_xmax(
            IdxRange<Grid1D...> idx_range) const override
    {
        return (*preallocate()).batched_derivs_idx_range_xmax(idx_range);
    }

    Field<double, IdxRange<Grid1D...>> operator()(
            Field<double, IdxRange<Grid1D...>> const inout_data,
            Field<const Coord<typename GridInterp::continuous_dimension_type>,
                  IdxRange<Grid1D...>> const coordinates,
            std::optional<Field<double const, batched_derivs_idx_range_type>> derivs_xmin
            = std::nullopt,
            std::optional<Field<double const, batched_derivs_idx_range_type>> derivs_xmax
            = std::nullopt) const override
    {
        return (*preallocate())(inout_data, coordinates);
    }
};
```


