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
        ddc::BoundCond BoundCond2max,
        class LayoutSpline = Kokkos::layout_right,
        class LayoutFunction = Kokkos::layout_right>
class SplineBuliderDerivField2D
{
    static_assert(
            std::is_same_v<ExecSpace, Kokkos::DefaultHostExecutionSpace>,
            "Implemented for DefaultHostExecutionSpace only yet.");

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

    using SplineType = DField<IdxRange<BSplines1, BSplines2>, MemorySpace, LayoutSpline>;
    using DerivFieldType = DerivField<
            double,
            IdxRange<Deriv1, Grid1, Deriv2, Grid2>,
            MemorySpace,
            LayoutFunction>;

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
        // Get the fields on the right layout.
        // --- get the index ranges.
        IdxRange<Grid1, Grid2> idx_range = get_idx_range(function_and_derivs.get_values_field());
        IdxRange<Grid1> idx_range_1(idx_range);
        IdxRange<Grid2> idx_range_2(idx_range);

        IdxRange<Deriv1> idx_range_d1(Idx<Deriv1>(1), IdxStep<Deriv1>(1));
        IdxRange<Deriv2> idx_range_d2(Idx<Deriv2>(1), IdxStep<Deriv2>(1));
        IdxRange<Deriv1, Deriv2> idx_range_d1d2(idx_range_d1, idx_range_d2);

        // --- allocate memory for fields on the correct layout.
        FunctFieldMem function_alloc(idx_range);

        IdxRange<Deriv1, Grid2> idx_range_d1_2(idx_range_d1, idx_range_2);
        Deriv1FieldMem deriv1_min_alloc(idx_range_d1_2);
        Deriv1FieldMem deriv1_max_alloc(idx_range_d1_2);

        IdxRange<Grid1, Deriv2> idx_range_1_d2(idx_range_1, idx_range_d2);
        Deriv2FieldMem deriv2_min_alloc(idx_range_1_d2);
        Deriv2FieldMem deriv2_max_alloc(idx_range_1_d2);

        CrossDerivFieldMem cross_min_min_alloc(idx_range_d1d2);
        CrossDerivFieldMem cross_max_min_alloc(idx_range_d1d2);
        CrossDerivFieldMem cross_min_max_alloc(idx_range_d1d2);
        CrossDerivFieldMem cross_max_max_alloc(idx_range_d1d2);

        // --- fill in the new fields with the data from the DerivField.
        fill_in_function(get_field(function_alloc), function_and_derivs);

        // If the boundary is not a ddc::BoundCond::HERMITE, we don't use derivatives.
        std::optional<Deriv1ConstField> deriv1_min_optional;
        std::optional<Deriv1ConstField> deriv1_max_optional;
        std::optional<Deriv2ConstField> deriv2_min_optional;
        std::optional<Deriv2ConstField> deriv2_max_optional;
        std::optional<CrossDerivConstField> cross_min_min_optional;
        std::optional<CrossDerivConstField> cross_max_min_optional;
        std::optional<CrossDerivConstField> cross_min_max_optional;
        std::optional<CrossDerivConstField> cross_max_max_optional;

        // --- fill in the first derivatives with the data from the DerivField.
        if constexpr (BoundCond1min == ddc::BoundCond::HERMITE) {
            fill_in_deriv1(get_field(deriv1_min_alloc), function_and_derivs, true);
            deriv1_min_optional = std::optional(get_const_field(deriv1_min_alloc));
        } else {
            deriv1_min_optional = std::optional<Deriv1ConstField> {std::nullopt};
        }

        if constexpr (BoundCond1max == ddc::BoundCond::HERMITE) {
            fill_in_deriv1(get_field(deriv1_max_alloc), function_and_derivs, false);
            deriv1_max_optional = std::optional(get_const_field(deriv1_max_alloc));
        } else {
            deriv1_max_optional = std::optional<Deriv1ConstField> {std::nullopt};
        }

        if constexpr (BoundCond2min == ddc::BoundCond::HERMITE) {
            fill_in_deriv2(get_field(deriv2_min_alloc), function_and_derivs, true);
            deriv2_min_optional = std::optional(get_const_field(deriv2_min_alloc));
        } else {
            deriv2_min_optional = std::optional<Deriv2ConstField> {std::nullopt};
        }

        if constexpr (BoundCond2max == ddc::BoundCond::HERMITE) {
            fill_in_deriv2(get_field(deriv2_max_alloc), function_and_derivs, false);
            deriv2_max_optional = std::optional(get_const_field(deriv2_max_alloc));
        } else {
            deriv2_max_optional = std::optional<Deriv2ConstField> {std::nullopt};
        }

        // --- fill in the cross-derivatives with the data from the DerivField.
        if constexpr (
                (BoundCond1min == ddc::BoundCond::HERMITE)
                && (BoundCond2min == ddc::BoundCond::HERMITE)) {
            fill_in_cross_deriv(get_field(cross_min_min_alloc), function_and_derivs, true, true);
            cross_min_min_optional = std::optional(get_const_field(cross_min_min_alloc));
        } else {
            cross_min_min_optional = std::optional<CrossDerivConstField> {std::nullopt};
        }

        if constexpr (
                (BoundCond1max == ddc::BoundCond::HERMITE)
                && (BoundCond2min == ddc::BoundCond::HERMITE)) {
            fill_in_cross_deriv(get_field(cross_max_min_alloc), function_and_derivs, false, true);
            cross_max_min_optional = std::optional(get_const_field(cross_max_min_alloc));
        } else {
            cross_max_min_optional = std::optional<CrossDerivConstField> {std::nullopt};
        }

        if constexpr (
                (BoundCond1min == ddc::BoundCond::HERMITE)
                && (BoundCond2max == ddc::BoundCond::HERMITE)) {
            fill_in_cross_deriv(get_field(cross_min_max_alloc), function_and_derivs, true, false);
            cross_min_max_optional = std::optional(get_const_field(cross_min_max_alloc));
        } else {
            cross_min_max_optional = std::optional<CrossDerivConstField> {std::nullopt};
        }

        if constexpr (
                (BoundCond1max == ddc::BoundCond::HERMITE)
                && (BoundCond2max == ddc::BoundCond::HERMITE)) {
            fill_in_cross_deriv(get_field(cross_max_max_alloc), function_and_derivs, false, false);
            cross_max_max_optional = std::optional(get_const_field(cross_max_max_alloc));
        } else {
            cross_max_max_optional = std::optional<CrossDerivConstField> {std::nullopt};
        }

        // Build the spline with the fields on the correct layout.
        m_builder(
                spline,
                get_const_field(function_alloc),
                deriv1_min_optional,
                deriv1_max_optional,
                deriv2_min_optional,
                deriv2_max_optional,
                cross_min_min_optional,
                cross_max_min_optional,
                cross_min_max_optional,
                cross_max_max_optional);
    };

private:
    /// @brief Fill in the function field with the values stored in the function_and_derivs.
    void fill_in_function(FunctField function, DerivFieldType function_and_derivs) const
    {
        // Extract data.
        DField<IdxRange<Grid1, Grid2>, MemorySpace, Kokkos::layout_stride> function_extracted
                = function_and_derivs.get_values_field();
        IdxRange<Grid1, Grid2> idx_range = get_idx_range(function_extracted);
        // Fill the field with correct layout.
        ddc::for_each(idx_range, [&](Idx<Grid1, Grid2> const idx) {
            function(idx) = function_extracted(idx);
        });
    }

    /**
     *  @brief Fill in the deriv1 field with the derivatives along the first dimension stored 
     * in the function_and_derivs. 
     * The boolean is_min helps to determine which bound we select for the derivative field. 
     */
    void fill_in_deriv1(Deriv1Field deriv1, DerivFieldType function_and_derivs, bool const is_min)
            const
    {
        // Extract data.
        IdxRange<Grid1, Grid2> idx_range = get_idx_range(function_and_derivs.get_values_field());
        IdxRange<Grid2> idx_range_2(idx_range);
        IdxRangeSlice<Grid1> idx_range_slice_1
                = function_and_derivs.template idx_range_for_deriv<Grid1>();

        Idx<Deriv1> idx_d1(Idx<Deriv1>(1));
        Idx<Grid1> idx_slice = is_min ? idx_range_slice_1.front() : idx_range_slice_1.back();
        Idx<Deriv1, Grid1> idx_deriv1(idx_d1, idx_slice);

        DField<IdxRange<Grid2>, MemorySpace, Kokkos::layout_stride> deriv1_extracted
                = function_and_derivs[idx_deriv1];

        // Fill the field with correct layout.
        ddc::for_each(idx_range_2, [&](Idx<Grid2> const idx) {
            deriv1(idx_d1, idx) = deriv1_extracted(idx);
        });
    }

    /**
     *  @brief Fill in the deriv2 field with the derivatives along the second dimension stored 
     * in the function_and_derivs. 
     * The boolean is_min helps to determine which bound we select for the derivative field. 
     */
    void fill_in_deriv2(Deriv2Field deriv2, DerivFieldType function_and_derivs, bool const is_min)
            const
    {
        // Extract data.
        IdxRange<Grid1, Grid2> idx_range = get_idx_range(function_and_derivs.get_values_field());
        IdxRange<Grid1> idx_range_1(idx_range);
        IdxRangeSlice<Grid2> idx_range_slice_2
                = function_and_derivs.template idx_range_for_deriv<Grid2>();

        Idx<Deriv2> idx_d2(Idx<Deriv2>(1));
        Idx<Grid2> idx_slice = is_min ? idx_range_slice_2.front() : idx_range_slice_2.back();
        Idx<Deriv2, Grid2> idx_deriv2(idx_d2, idx_slice);

        DField<IdxRange<Grid1>, MemorySpace, Kokkos::layout_stride> deriv2_extracted
                = function_and_derivs[idx_deriv2];

        // Fill the field with correct layout.
        ddc::for_each(idx_range_1, [&](Idx<Grid1> const idx) {
            deriv2(idx, idx_d2) = deriv2_extracted(idx);
        });
    }

    /**
     *  @brief Fill in the cross_deriv field with the cross derivatives stored 
     * in the function_and_derivs. 
     * The booleans is_1min and is_2min help to determine which corner we select 
     * for the derivative field. 
     */
    void fill_in_cross_deriv(
            CrossDerivField cross_deriv,
            DerivFieldType function_and_derivs,
            bool const is_1min,
            bool const is_2min) const
    {
        // Extract data.
        IdxRange<Grid1, Grid2> idx_range = get_idx_range(function_and_derivs.get_values_field());
        IdxRange<Grid1> idx_range_1(idx_range);
        IdxRange<Grid2> idx_range_2(idx_range);

        IdxRangeSlice<Grid1> idx_range_slice_1
                = function_and_derivs.template idx_range_for_deriv<Grid1>();
        IdxRangeSlice<Grid2> idx_range_slice_2
                = function_and_derivs.template idx_range_for_deriv<Grid2>();

        Idx<Deriv1> idx_d1(Idx<Deriv1>(1));
        Idx<Deriv2> idx_d2(Idx<Deriv2>(1));
        Idx<Grid1> idx_slice_1 = is_1min ? idx_range_slice_1.front() : idx_range_slice_1.back();
        Idx<Grid2> idx_slice_2 = is_2min ? idx_range_slice_2.front() : idx_range_slice_2.back();

        Idx<Deriv1, Grid1, Deriv2, Grid2> idx_cross_deriv(idx_d1, idx_slice_1, idx_d2, idx_slice_2);

        // Fill the field with correct layout.
        cross_deriv(idx_d1, idx_d2) = function_and_derivs(idx_cross_deriv);
    }
};
