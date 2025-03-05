
// SPDX-License-Identifier: MIT
/**
 * @file spline_builder_2d_cache.hpp
 * File containing a class to store the spline builder coefficients
 */

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"


/**
 * @brief A class which represents a cache for spline builder coefficients
 */
template <class SplineBuilder2D, class Dim1, class Dim2>
class SplineBuilder2DCache
{
public:
private:
    // Type for spline representation of the field
    using IdxRangeBSField = typename SplineBuilder2D::batched_spline_domain_type;
    using IdxRangeField = typename SplineBuilder2D::batched_interpolation_domain_type;
    using DFieldSplineCoeffMem = DFieldMem<IdxRangeBSField>;
    using DConstFieldSplineCoeffs = DConstField<IdxRangeBSField>;

    SplineBuilder2D const& m_spline_builder;
    DFieldSplineCoeffMem m_spline_coeffs;
    bool m_compute_coeffs_dim1;
    bool m_compute_coeffs_dim2;

public:
    /**
    * @brief Construct an instance of the class SplineBuilderCache.
    *
    * @param spline_builder A 2D spline builder.
    */
    SplineBuilderCache(
            SplineBuilder2D const& spline_builder)
        : m_spline_builder(spline_builder)
        , m_spline_coeffs(spline_builder.batched_interpolation_domain())
        , m_compute_coeffs_dim1(true)
        , m_compute_coeffs_dim2(true)
    {
    }

    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[out] dfieldval_dxi Partial derivatives in Xi direction.
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    template <class DimOfInterest>
    DConstFieldSplineCoeffs compute_coeffs(DConstFieldMem<IdxRangeField> field_values)
    {
        static_assert( std::is_same_v<DimOfInterest, Dim1> || std::is_same_v<DimOfInterest, Dim2> )

        if ( std::is_same_v<DimOfInterest, Dim1> && m_compute_coeffs_dim1 ) {
            m_compute_coeffs_dim2 = false;
            m_compute_coeffs_dim1 = true;
            m_spline_builder(get_field(m_spline_coeffs), field_values));

      } else if ( std::is_same_v<DimOfInterest, Dim2> && m_compute_coeffs_dim2 ) { 
            m_compute_coeffs_dim1 = false;
            m_compute_coeffs_dim2 = true;
            m_spline_builder(get_field(m_spline_coeffs), field_values));

    }
        return get_const_field(m_spline_coeffs)
    }
};
