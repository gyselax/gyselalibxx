

# File spline\_interpolator\_2d.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**spline\_interpolator\_2d.hpp**](spline__interpolator__2d_8hpp.md)

[Go to the documentation of this file](spline__interpolator__2d_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolator_2d.hpp"

template <class Spline2DBuilder, class Spline2DEvaluator>
class SplineInterpolator2D
    : public IInterpolator2D<
              typename Spline2DBuilder::interpolation_domain_type,
              typename Spline2DBuilder::batched_interpolation_domain_type>
{
    using base_type = IInterpolator2D<
            typename Spline2DBuilder::interpolation_domain_type,
            typename Spline2DBuilder::batched_interpolation_domain_type>;

public:
    using typename base_type::CConstFieldType;
    using typename base_type::CoordType;
    using typename base_type::DFieldType;

private:
    Spline2DBuilder const& m_builder;

    Spline2DEvaluator const& m_evaluator;

    using IdxRangeBSRTheta = typename Spline2DBuilder::batched_spline_domain_type;

    mutable DFieldMem<IdxRangeBSRTheta> m_coefs;

    using r_deriv_type = DConstField<typename Spline2DBuilder::batched_derivs_domain_type1>;
    using theta_deriv_type = DConstField<typename Spline2DBuilder::batched_derivs_domain_type2>;
    using mixed_deriv_type = DConstField<typename Spline2DBuilder::batched_derivs_domain_type>;

public:
    SplineInterpolator2D(Spline2DBuilder const& builder, Spline2DEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(get_spline_idx_range(builder))
    {
    }

    DFieldType operator()(DFieldType const inout_data, CConstFieldType const coordinates)
            const override
    {
        m_builder(get_field(m_coefs), get_const_field(inout_data));
        m_evaluator(get_field(inout_data), coordinates, get_const_field(m_coefs));

        return inout_data;
    }
};



template <class Spline2DBuilder, class Spline2DEvaluator>
class PreallocatableSplineInterpolator2D
    : public IPreallocatableInterpolator2D<
              typename Spline2DBuilder::interpolation_domain_type,
              typename Spline2DBuilder::batched_interpolation_domain_type>
{
private:
    Spline2DBuilder const& m_builder;

    Spline2DEvaluator const& m_evaluator;

public:
    PreallocatableSplineInterpolator2D(
            Spline2DBuilder const& builder,
            Spline2DEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    std::unique_ptr<IInterpolator2D<
            typename Spline2DBuilder::interpolation_domain_type,
            typename Spline2DBuilder::batched_interpolation_domain_type>>
    preallocate() const override
    {
        return std::make_unique<
                SplineInterpolator2D<Spline2DBuilder, Spline2DEvaluator>>(m_builder, m_evaluator);
    }
};
```


