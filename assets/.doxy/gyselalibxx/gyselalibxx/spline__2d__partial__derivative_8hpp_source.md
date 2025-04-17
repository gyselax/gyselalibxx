

# File spline\_2d\_partial\_derivative.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**spline\_2d\_partial\_derivative.hpp**](spline__2d__partial__derivative_8hpp.md)

[Go to the documentation of this file](spline__2d__partial__derivative_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"
#include "spline_builder_2d_cache.hpp"


template <
        class SplineBuilder2D,
        class SplineEvaluator2D,
        class DerivativeDimension,
        class IdxRangeBatched>
class Spline2DPartialDerivative : public IPartialDerivative<IdxRangeBatched, DerivativeDimension>
{
    static_assert(
            (std::is_same_v<
                    typename SplineEvaluator2D::continuous_dimension_type1,
                    DerivativeDimension>)
            || (std::is_same_v<
                    typename SplineEvaluator2D::continuous_dimension_type2,
                    DerivativeDimension>));

private:
    using base_type = IPartialDerivative<IdxRangeBatched, DerivativeDimension>;

    using typename base_type::DConstFieldType;
    using typename base_type::DFieldType;

    using IdxRangeBS = IdxRangeBatched;
    using DFieldBSMem = DFieldMem<IdxRangeBS>;
    using DFieldBS = DField<IdxRangeBS>;

    SplineBuilder2DCache<SplineBuilder2D, IdxRangeBatched>& m_builder_cache;
    SplineEvaluator2D const& m_evaluator;

public:
    explicit Spline2DPartialDerivative(
            SplineBuilder2DCache<SplineBuilder2D, IdxRangeBatched>& builder_cache,
            SplineEvaluator2D const& evaluator,
            DConstFieldType const field)
        : m_builder_cache(builder_cache)
        , m_evaluator(evaluator)
    {
        m_builder_cache.template compute_coeffs<DerivativeDimension>(field);
    }

    void operator()(DFieldType differentiated_field) const final
    {
        if constexpr (std::is_same_v<
                              DerivativeDimension,
                              typename SplineEvaluator2D::continuous_dimension_type1>) {
            m_evaluator.deriv_dim_1(differentiated_field, m_builder_cache());
        } else {
            m_evaluator.deriv_dim_2(differentiated_field, m_builder_cache());
        }
    }
};

template <
        class SplineBuilder2D,
        class SplineEvaluator2D,
        class DerivativeDimension,
        class IdxRangeBatched>
class Spline2DPartialDerivativeCreator
    : public IPartialDerivativeCreator<IdxRangeBatched, DerivativeDimension>
{
private:
    using DConstFieldType = DConstField<IdxRangeBatched>;

    SplineBuilder2DCache<SplineBuilder2D, IdxRangeBatched>& m_builder_cache;
    SplineEvaluator2D const& m_evaluator;

public:
    Spline2DPartialDerivativeCreator(
            SplineBuilder2DCache<SplineBuilder2D, IdxRangeBatched>& builder_cache,
            SplineEvaluator2D const& evaluator)
        : m_builder_cache(builder_cache)
        , m_evaluator(evaluator)
    {
    }

    std::unique_ptr<IPartialDerivative<IdxRangeBatched, DerivativeDimension>> create_instance(
            DConstFieldType field) const
    {
        return std::make_unique<Spline2DPartialDerivative<
                SplineBuilder2D,
                SplineEvaluator2D,
                DerivativeDimension,
                IdxRangeBatched>>(m_builder_cache, m_evaluator, field);
    }
};
```


