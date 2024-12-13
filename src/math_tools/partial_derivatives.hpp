
// SPDX-License-Identifier: MIT
/**
 * @file partial_derivatives.hpp
 * File containing functions to compute the partial derivatives
 */

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"


/**
 * @brief Compute the partial derivative of a 
 *
 * For a given field @f$F(X1,X2)@f$ , compute
 * @f$\partial_{Xi} F(X1,X2)@f$ with i=1 or 2.
 *
 * @param[in] Xi
 *      The given coordinate for the partial derivative.
 *
 * @return @f$\partial_{Xi} F(X1,X2)@f$ 
 */

template <
        class GridXi,
        class GridX1,
        class GridX2,
        class FieldX1Builder,
        class FieldX1Evaluator>
class PartialDerivative
{
private:
    using IdxRangeX1X2 = IdxRange<GridX1, GridX2>;
    using IdxX1X2 = typename IdxRangeX1X2::discrete_element_type;
    using DFieldX1X2 = DField<IdxRangeX1X2>;
    using DConstFieldX1X2 = DConstField<IdxRangeX1X2>;

    // Type for spline representation of the field
    using IdxRangeBSFieldX1 = typename FieldX1Builder::batched_spline_domain_type;
    using FieldX1SplineMem = DFieldMem<IdxRangeBSFieldX1>;
    using FieldX1SplineCoeffs = DField<IdxRangeBSFieldX1>;

    FieldX1Builder const& m_fieldx1_builder;
    FieldX1Evaluator const& m_fieldx1_evaluator;

public:
    explicit PartialDerivative(
        FieldX1Builder const& fieldx1_builder,
        FieldX1Evaluator const& fieldx1_evaluator)
        : m_fieldx1_builder(fieldx1_builder)
        , m_fieldx1_evaluator(fieldx1_evaluator)             
    {    
    }

    ~PartialDerivative() = default;

    DFieldX1X2 operator()(DFieldX1X2 dfield_dx1_x1x2, DConstFieldX1X2 field_x1x2)
    {
        // Build spline representation of the field ....................................
        FieldX1SplineMem fieldx1_coefs_alloc(
            m_fieldx1_builder.batched_spline_domain());
        FieldX1SplineCoeffs fieldx1_coefs = get_field(fieldx1_coefs_alloc);

        m_fieldx1_builder(fieldx1_coefs, get_const_field(field_x1x2));
        m_fieldx1_evaluator.deriv(dfield_dx1_x1x2, get_const_field(fieldx1_coefs));

        return dfield_dx1_x1x2;
    }
};
