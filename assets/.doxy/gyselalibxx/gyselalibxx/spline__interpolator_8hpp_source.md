

# File spline\_interpolator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**spline\_interpolator.hpp**](spline__interpolator_8hpp.md)

[Go to the documentation of this file](spline__interpolator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "iinterpolator.hpp"

template <
        class GridInterp,
        class BSplines,
        ddc::BoundCond BcMin,
        ddc::BoundCond BcMax,
        class LeftExtrapolationRule,
        class RightExtrapolationRule,
        ddc::SplineSolver Solver,
        class... Grid1D>
class SplineInterpolator : public IInterpolator<GridInterp, Grid1D...>
{
    using BuilderType = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            GridInterp,
            BcMin,
            BcMax,
            Solver>;
    using EvaluatorType = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            GridInterp,
            LeftExtrapolationRule,
            RightExtrapolationRule>;
    using deriv_type = typename IInterpolator<GridInterp, Grid1D...>::deriv_type;
    using batched_spline_domain_type =
            typename BuilderType::batched_spline_domain_type<IdxRange<Grid1D...>>;
    using batched_derivs_idx_range_type =
            typename IInterpolator<GridInterp, Grid1D...>::batched_derivs_idx_range_type;
    using batched_deriv_field_type = ConstField<double, batched_derivs_idx_range_type>;

private:
    BuilderType const& m_builder;

    EvaluatorType const& m_evaluator;

    mutable DFieldMem<batched_spline_domain_type> m_coefs;


public:
    SplineInterpolator(
            BuilderType const& builder,
            EvaluatorType const& evaluator,
            IdxRange<Grid1D...> idx_range_batched)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(builder.batched_spline_domain(idx_range_batched))
    {
    }

    ~SplineInterpolator() override = default;

    batched_derivs_idx_range_type batched_derivs_idx_range_xmin(
            IdxRange<Grid1D...> idx_range) const override
    {
        return ddc::replace_dim_of<GridInterp, deriv_type>(
                idx_range,
                IdxRange<deriv_type>(
                        Idx<deriv_type>(1),
                        IdxStep<deriv_type>(BuilderType::s_nbc_xmin)));
    }

    batched_derivs_idx_range_type batched_derivs_idx_range_xmax(
            IdxRange<Grid1D...> idx_range) const override
    {
        return ddc::replace_dim_of<GridInterp, deriv_type>(
                idx_range,
                IdxRange<deriv_type>(
                        Idx<deriv_type>(1),
                        IdxStep<deriv_type>(BuilderType::s_nbc_xmax)));
    }

    Field<double, IdxRange<Grid1D...>> operator()(
            Field<double, IdxRange<Grid1D...>> const inout_data,
            ConstField<
                    Coord<typename GridInterp::continuous_dimension_type>,
                    IdxRange<Grid1D...>> const coordinates,
            std::optional<batched_deriv_field_type> derivs_xmin = std::nullopt,
            std::optional<batched_deriv_field_type> derivs_xmax = std::nullopt) const override
    {
        m_builder(get_field(m_coefs), get_const_field(inout_data), derivs_xmin, derivs_xmax);
        m_evaluator(inout_data, coordinates, get_const_field(m_coefs));
        return inout_data;
    }
};

template <
        class GridInterp,
        class BSplines,
        ddc::BoundCond BcMin,
        ddc::BoundCond BcMax,
        class LeftExtrapolationRule,
        class RightExtrapolationRule,
        ddc::SplineSolver Solver,
        class... Grid1D>
class PreallocatableSplineInterpolator : public IPreallocatableInterpolator<GridInterp, Grid1D...>
{
    using BuilderType = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            GridInterp,
            BcMin,
            BcMax,
            Solver>;
    using EvaluatorType = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            GridInterp,
            LeftExtrapolationRule,
            RightExtrapolationRule>;

    BuilderType const& m_builder;

    EvaluatorType const& m_evaluator;

    IdxRange<Grid1D...> m_idx_range_batched;

public:
    PreallocatableSplineInterpolator(
            BuilderType const& builder,
            EvaluatorType const& evaluator,
            IdxRange<Grid1D...> idx_range_batched)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_idx_range_batched(idx_range_batched)
    {
    }

    ~PreallocatableSplineInterpolator() override = default;

    std::unique_ptr<IInterpolator<GridInterp, Grid1D...>> preallocate() const override
    {
        return std::make_unique<SplineInterpolator<
                GridInterp,
                BSplines,
                BcMin,
                BcMax,
                LeftExtrapolationRule,
                RightExtrapolationRule,
                Solver,
                Grid1D...>>(m_builder, m_evaluator, m_idx_range_batched);
    }
};
```


