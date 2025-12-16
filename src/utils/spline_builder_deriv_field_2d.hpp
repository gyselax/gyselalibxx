// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"

/**
 * @brief [Temporary] Apply a SplineBuilder2D to a DerivField.
 * 
 * DerivField stores the values and the derivatives of a function. 
 * The inputs of the SplineBuilder2D are on different layouts than the fields
 * we can get from the DerivField. SplineBuliderDerivField2D allows to 
 * directly apply a stored SplineBuilder2D to a DerivField by copying 
 * data in fields with correct layout. 
 * 
 * Implemented only for 2D case. 
 */
template <
        class ExecSpace,
        class BSplines1,
        class BSplines2,
        class Grid1,
        class Grid2,
        ddc::BoundCond BoundCond1min,
        ddc::BoundCond BoundCond1max,
        ddc::BoundCond BoundCond2min,
        ddc::BoundCond BoundCond2max>
class SplineBuliderDerivField2D
{
    using MemorySpace = typename ExecSpace::memory_space;

    using Dim1 = typename BSplines1::continuous_dimension_type;
    using Dim2 = typename BSplines2::continuous_dimension_type;
    using Deriv1 = ddc::Deriv<Dim1>;
    using Deriv2 = ddc::Deriv<Dim2>;

    using Builder2D = ddc::SplineBuilder2D<
            ExecSpace,
            MemorySpace,
            BSplines1,
            BSplines2,
            Grid1,
            Grid2,
            BoundCond1min,
            BoundCond1max,
            BoundCond2min,
            BoundCond2max,
            ddc::SplineSolver::LAPACK>;

    using SplineType = DField<IdxRange<BSplines1, BSplines2>, MemorySpace>;
    using DerivFieldType = DerivField<double, IdxRange<Deriv1, Grid1, Deriv2, Grid2>, MemorySpace>;

    using FunctFieldMem = DFieldMem<IdxRange<Grid1, Grid2>, MemorySpace>;
    using FunctField = DField<IdxRange<Grid1, Grid2>, MemorySpace>;

    using CrossDerivFieldMem = DFieldMem<IdxRange<Deriv1, Deriv2>, MemorySpace>;
    using CrossDerivField = DField<IdxRange<Deriv1, Deriv2>, MemorySpace>;
    using CrossDerivConstField = DConstField<IdxRange<Deriv1, Deriv2>, MemorySpace>;

    using Deriv1FieldMem = DFieldMem<IdxRange<Deriv1, Grid2>, MemorySpace>;
    using Deriv1Field = DField<IdxRange<Deriv1, Grid2>, MemorySpace>;
    using Deriv1ConstField = DConstField<IdxRange<Deriv1, Grid2>, MemorySpace>;

    using Deriv2FieldMem = DFieldMem<IdxRange<Grid1, Deriv2>, MemorySpace>;
    using Deriv2Field = DField<IdxRange<Grid1, Deriv2>, MemorySpace>;
    using Deriv2ConstField = DConstField<IdxRange<Grid1, Deriv2>, MemorySpace>;

private:
    Builder2D const& m_builder;

public:
    /**
     * @brief Instantiate the class by storing a reference to the spline builder 
     * we can to use. 
     * 
     * @param[in] builder A reference to SplineBuilder2D from DDC that we store
     * in the class. 
     */
    explicit SplineBuliderDerivField2D(Builder2D const& builder) : m_builder(builder) {}

    /**
     * @brief Build the spline representation of the function given in a 
     * DerivField applying the referenced spline builder stored in the class. 
     * 
     * @param[out] spline Spline coefficients on a 2D grid. 
     * @param[in] function_and_derivs Data defining the function on a 2D grid. 
     */
    void operator()(SplineType spline, DerivFieldType function_and_derivs) const
    {
        // Check that the DerivField contains the necessary derivatives for the builder.
        IdxRange<Deriv1> idx_range_d1_min(
                Idx<Deriv1>(1),
                IdxStep<Deriv1>(Builder2D::builder_type1::s_nbc_xmin));
        IdxRange<Deriv1> idx_range_d1_max(
                Idx<Deriv1>(1),
                IdxStep<Deriv1>(Builder2D::builder_type1::s_nbc_xmax));
        IdxRange<Deriv2> idx_range_d2_min(
                Idx<Deriv2>(1),
                IdxStep<Deriv2>(Builder2D::builder_type2::s_nbc_xmin));
        IdxRange<Deriv2> idx_range_d2_max(
                Idx<Deriv2>(1),
                IdxStep<Deriv2>(Builder2D::builder_type2::s_nbc_xmax));

        IdxRange<Deriv1, Deriv2> idx_range_d1d2_min_min(idx_range_d1_min, idx_range_d2_min);
        IdxRange<Deriv1, Deriv2> idx_range_d1d2_max_min(idx_range_d1_max, idx_range_d2_min);
        IdxRange<Deriv1, Deriv2> idx_range_d1d2_min_max(idx_range_d1_min, idx_range_d2_max);
        IdxRange<Deriv1, Deriv2> idx_range_d1d2_max_max(idx_range_d1_max, idx_range_d2_max);

        // Get the fields on the right layout.
        // --- get the index ranges.
        IdxRange<Grid1, Grid2> idx_range = get_idx_range(function_and_derivs);
        IdxRange<Grid1> idx_range_1(idx_range);
        IdxRange<Grid2> idx_range_2(idx_range);

        // --- allocate memory for fields on the correct layout.
        FunctFieldMem function_alloc(idx_range);

        IdxRange<Deriv1, Grid2> idx_range_d1_2_min(idx_range_d1_min, idx_range_2);
        IdxRange<Deriv1, Grid2> idx_range_d1_2_max(idx_range_d1_max, idx_range_2);
        Deriv1FieldMem deriv1_min_alloc(idx_range_d1_2_min);
        Deriv1FieldMem deriv1_max_alloc(idx_range_d1_2_max);

        IdxRange<Grid1, Deriv2> idx_range_1_d2_min(idx_range_1, idx_range_d2_min);
        IdxRange<Grid1, Deriv2> idx_range_1_d2_max(idx_range_1, idx_range_d2_max);
        Deriv2FieldMem deriv2_min_alloc(idx_range_1_d2_min);
        Deriv2FieldMem deriv2_max_alloc(idx_range_1_d2_max);

        CrossDerivFieldMem cross_min_min_alloc(idx_range_d1d2_min_min);
        CrossDerivFieldMem cross_max_min_alloc(idx_range_d1d2_max_min);
        CrossDerivFieldMem cross_min_max_alloc(idx_range_d1d2_min_max);
        CrossDerivFieldMem cross_max_max_alloc(idx_range_d1d2_max_max);

        // --- fill in the new fields with the data from the DerivField.
        fill_in_function(get_field(function_alloc), function_and_derivs);

        // --- fill in the first derivatives with the data from the DerivField.
        fill_in_deriv1(get_field(deriv1_min_alloc), function_and_derivs, true);
        fill_in_deriv1(get_field(deriv1_max_alloc), function_and_derivs, false);
        fill_in_deriv2(get_field(deriv2_min_alloc), function_and_derivs, true);
        fill_in_deriv2(get_field(deriv2_max_alloc), function_and_derivs, false);

        // --- fill in the cross-derivatives with the data from the DerivField.
        fill_in_cross_deriv(get_field(cross_min_min_alloc), function_and_derivs, true, true);
        fill_in_cross_deriv(get_field(cross_max_min_alloc), function_and_derivs, false, true);
        fill_in_cross_deriv(get_field(cross_min_max_alloc), function_and_derivs, true, false);
        fill_in_cross_deriv(get_field(cross_max_max_alloc), function_and_derivs, false, false);

        // Build the spline with the fields on the correct layout.
        m_builder(
                spline,
                get_const_field(function_alloc),
                std::optional(get_const_field(deriv1_min_alloc)),
                std::optional(get_const_field(deriv1_max_alloc)),
                std::optional(get_const_field(deriv2_min_alloc)),
                std::optional(get_const_field(deriv2_max_alloc)),
                std::optional(get_const_field(cross_min_min_alloc)),
                std::optional(get_const_field(cross_max_min_alloc)),
                std::optional(get_const_field(cross_min_max_alloc)),
                std::optional(get_const_field(cross_max_max_alloc)));
    };

public:
    /**
     * @brief Fill in the function field with the values stored in the function_and_derivs.
     * @param[out] function Field with layout_right where we copy the function values.
     * @param[in] function_and_derivs DerivField from where the function values are copied.
     */
    void fill_in_function(FunctField function, DerivFieldType function_and_derivs) const
    {
        // Fill the field with correct layout.
        ddc::parallel_for_each(
                ExecSpace(),
                get_idx_range(function_and_derivs),
                KOKKOS_LAMBDA(Idx<Grid1, Grid2> const idx) {
                    function(idx) = function_and_derivs(idx);
                });
    }

    /**
     *  @brief Fill in the deriv1 field with the derivatives along the first dimension stored 
     * in the function_and_derivs. 
     * @param[out] deriv1 Field with layout_right where we copy the derivatives.
     * @param[in] function_and_derivs DerivField from where the derivatives are copied.
     * @param[in] is_min Boolean which helps to determine which bound (mon/max) we select for the derivative field. 
     */
    void fill_in_deriv1(Deriv1Field deriv1, DerivFieldType function_and_derivs, bool const is_min)
            const
    {
        IdxRangeSlice<Grid1> idx_range_slice_1
                = function_and_derivs.template idx_range_for_deriv<Grid1>();
        Idx<Grid1> idx_slice = is_min ? idx_range_slice_1.front() : idx_range_slice_1.back();

        // Fill the field with correct layout.
        ddc::parallel_for_each(
                ExecSpace(),
                get_idx_range(deriv1),
                KOKKOS_LAMBDA(Idx<Deriv1, Grid2> const idx) {
                    deriv1(idx) = function_and_derivs(idx, idx_slice);
                });
    }

    /**
     *  @brief Fill in the deriv2 field with the derivatives along the second dimension stored 
     * in the function_and_derivs. 
     * @param[out] deriv2 Field with layout_right where we copy the derivatives.
     * @param[in] function_and_derivs DerivField from where the derivatives are copied.
     * @param[in] is_min Boolean which helps to determine which bound (mon/max) we select for the derivative field. 
     */
    void fill_in_deriv2(Deriv2Field deriv2, DerivFieldType function_and_derivs, bool const is_min)
            const
    {
        IdxRangeSlice<Grid2> idx_range_slice_2
                = function_and_derivs.template idx_range_for_deriv<Grid2>();
        Idx<Grid2> idx_slice = is_min ? idx_range_slice_2.front() : idx_range_slice_2.back();

        // Fill the field with correct layout.
        ddc::parallel_for_each(
                ExecSpace(),
                get_idx_range(deriv2),
                KOKKOS_LAMBDA(Idx<Grid1, Deriv2> const idx) {
                    deriv2(idx) = function_and_derivs(idx, idx_slice);
                });
    }

    /**
     *  @brief Fill in the cross_deriv field with the cross derivatives stored 
     * in the function_and_derivs. 
     * @param[out] cross_deriv Field with layout_right where we copy the cross-derivatives.
     * @param[in] function_and_derivs DerivField from where the cross-derivatives are copied.
     * @param[in] is_1min Boolean which helps to determine which bound (mon/max) we select for the derivative field
     * on the first dimension. 
     * @param[in] is_2min Boolean which helps to determine which bound (mon/max) we select for the derivative field
     * on the second dimension. 
     */
    void fill_in_cross_deriv(
            CrossDerivField cross_deriv,
            DerivFieldType function_and_derivs,
            bool const is_1min,
            bool const is_2min) const
    {
        IdxRangeSlice<Grid1> idx_range_slice_1
                = function_and_derivs.template idx_range_for_deriv<Grid1>();
        IdxRangeSlice<Grid2> idx_range_slice_2
                = function_and_derivs.template idx_range_for_deriv<Grid2>();

        Idx<Grid1> idx_slice_1 = is_1min ? idx_range_slice_1.front() : idx_range_slice_1.back();
        Idx<Grid2> idx_slice_2 = is_2min ? idx_range_slice_2.front() : idx_range_slice_2.back();

        // Fill the field with correct layout.
        ddc::parallel_for_each(
                get_idx_range(cross_deriv),
                KOKKOS_LAMBDA(Idx<Deriv1, Deriv2> idx_derivs) {
                    cross_deriv(idx_derivs)
                            = function_and_derivs(idx_derivs, idx_slice_1, idx_slice_2);
                });
    }
};
