

# File spline\_builder\_2d\_cache.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**spline\_builder\_2d\_cache.hpp**](spline__builder__2d__cache_8hpp.md)

[Go to the documentation of this file](spline__builder__2d__cache_8hpp.md)


```C++

// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"


template <class SplineBuilder2D, class IdxRangeBatched>
class SplineBuilder2DCache
{
private:
    using Dim1 = typename SplineBuilder2D::continuous_dimension_type1;
    using Dim2 = typename SplineBuilder2D::continuous_dimension_type2;

    using IdxRangeBSField = typename SplineBuilder2D::batched_spline_domain_type<IdxRangeBatched>;
    using IdxRangeField = IdxRangeBatched;

    using DFieldSplineCoeffMem = DFieldMem<IdxRangeBSField>;
    using DConstFieldSplineCoeffs = DConstField<IdxRangeBSField>;

    SplineBuilder2D const& m_spline_builder;
    DFieldSplineCoeffMem m_spline_coeffs;
    bool m_compute_coeffs_dim1;
    bool m_compute_coeffs_dim2;

public:
    explicit SplineBuilder2DCache(
            SplineBuilder2D const& spline_builder,
            IdxRangeBatched idx_range_batched)
        : m_spline_builder(spline_builder)
        , m_spline_coeffs(spline_builder.batched_spline_domain(idx_range_batched))
        , m_compute_coeffs_dim1(true)
        , m_compute_coeffs_dim2(true)
    {
    }

    template <class DimOfInterest>
    DConstFieldSplineCoeffs compute_coeffs(DConstField<IdxRangeField> field_values)
    {
        static_assert(std::is_same_v<DimOfInterest, Dim1> || std::is_same_v<DimOfInterest, Dim2>);
        if (std::is_same_v<DimOfInterest, Dim1> && m_compute_coeffs_dim1) {
            m_compute_coeffs_dim2 = false;
            m_compute_coeffs_dim1 = true;
            m_spline_builder(get_field(m_spline_coeffs), field_values);

        } else if (std::is_same_v<DimOfInterest, Dim2> && m_compute_coeffs_dim2) {
            m_compute_coeffs_dim1 = false;
            m_compute_coeffs_dim2 = true;
            m_spline_builder(get_field(m_spline_coeffs), field_values);
        }

        return get_const_field(m_spline_coeffs);
    }

    DConstFieldSplineCoeffs operator()() const
    {
        return get_const_field(m_spline_coeffs);
    }
};
```


