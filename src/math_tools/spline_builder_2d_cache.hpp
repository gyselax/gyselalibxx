
// SPDX-License-Identifier: MIT
/**
 * @file spline_builder_2d_cache.hpp
 * File containing a class to store the spline builder coefficients and recompute them when required..
 */

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"


/**
 * @brief A class that stores spline builder coefficients and recomputes them when required.
 */
template <class SplineBuilder2D>
class SplineBuilder2DCache
{
private:
    using Dim1 = typename SplineBuilder2D::continuous_dimension_type1;
    using Dim2 = typename SplineBuilder2D::continuous_dimension_type2;

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
    * @brief Construct an instance of the class SplineBuilder2DCache.
    *
    * @param spline_builder A 2D spline builder.
    */
    SplineBuilder2DCache(SplineBuilder2D const& spline_builder)
        : m_spline_builder(spline_builder)
        , m_spline_coeffs(spline_builder.batched_spline_domain())
        , m_compute_coeffs_dim1(true)
        , m_compute_coeffs_dim2(true)
    {
    }

    /**
    * @brief Compute the spline coefficients of the spline representation of field_values
    * if required, i.e. if not already computed by the other dimension than DimOfInterest.
    *
    * @param[in] field_values The field to be used to compute spline coefficients.
    *
    * @return The spline coefficients updated if required.
    */
    template <class DimOfInterest>
    DConstFieldSplineCoeffs compute_coeffs(DConstField<IdxRangeField> field_values)
    {
        static_assert(std::is_same_v<DimOfInterest, Dim1> || std::is_same_v<DimOfInterest, Dim2>);
        /**
        * The compute_coeff function can be called with either Dim1 or Dim2 as a template parameter. 
        * When a call is made with one of these two dimensions, the corresponding boolean is set to 
        * false so that when the function is called with the other dimension, the builder is not 
        * called again. 
        *
        * To summarise with X and Y as the two dimensions: 
        * ITERATION N 
        *   compute_coeff<X>(field) --> calls spline_builder
        *   compute_coeff<Y>(field) --> does not call spline_builder
        *
        * ITERATION N+1 
        *   compute_coeff<X>(field_updated) --> calls spline_builder
        *   compute_coeff<Y>(field_updated) --> does not call spline_builder
        */
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
};
