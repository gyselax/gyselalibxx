

# File bsl\_advection\_1d.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**bsl\_advection\_1d.hpp**](bsl__advection__1d_8hpp.md)

[Go to the documentation of this file](bsl__advection__1d_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "euler.hpp"
#include "iinterpolator.hpp"


template <
        class GridInterest,
        class IdxRangeAdvection,
        class IdxRangeFunction,
        class AdvectionFieldBuilder,
        class AdvectionFieldEvaluator,
        class TimeStepper
        = Euler<FieldMem<
                        Coord<typename GridInterest::continuous_dimension_type>,
                        IdxRangeAdvection>,
                DFieldMem<IdxRangeAdvection>>>
class BslAdvection1D
{
private:
    // Advection index range element:
    using IdxAdvection = typename IdxRangeAdvection::discrete_element_type;

    // Full index range element:
    using IdxFunction = typename IdxRangeFunction::discrete_element_type;

    // Advection dimension (or Interest dimension):
    using DimInterest = typename GridInterest::continuous_dimension_type;
    using CoordInterest = Coord<DimInterest>;
    using IdxRangeInterest = IdxRange<GridInterest>;
    using IdxInterest = typename IdxRangeInterest::discrete_element_type;

    // Type for the feet and advection field:
    using FeetFieldMem = FieldMem<CoordInterest, IdxRangeAdvection>;
    using FeetField = typename FeetFieldMem::span_type;
    using FeetConstField = typename FeetFieldMem::view_type;

    using AdvecFieldMem = DFieldMem<IdxRangeAdvection>;
    using AdvecField = typename AdvecFieldMem::span_type;

    using FunctionField = DField<IdxRangeFunction>;

    // Type for spline representation of the advection field
    using IdxRangeBSAdvection = typename AdvectionFieldBuilder::batched_spline_domain_type;
    using AdvecFieldSplineMem = DFieldMem<IdxRangeBSAdvection>;
    using AdvecFieldSplineCoeffs = DField<IdxRangeBSAdvection>;

    // Type for the derivatives of the advection field
    using DerivDim = ddc::Deriv<DimInterest>;
    using IdxRangeAdvecFieldDeriv
            = ddc::replace_dim_of_t<IdxRangeAdvection, GridInterest, DerivDim>;
    using AdvecFieldDerivConstField = Field<const double, IdxRangeAdvecFieldDeriv>;


    // Interpolators:
    using FunctionPreallocatableInterpolatorType = interpolator_on_idx_range_t<
            IPreallocatableInterpolator,
            GridInterest,
            IdxRangeFunction>;
    using FunctionInterpolatorType
            = interpolator_on_idx_range_t<IInterpolator, GridInterest, IdxRangeFunction>;

    // Type for the derivatives of the function
    using IdxRangeFunctionDeriv = typename FunctionInterpolatorType::batched_derivs_idx_range_type;
    using FunctionDerivFieldMem = DFieldMem<IdxRangeFunctionDeriv>;

    FunctionPreallocatableInterpolatorType const& m_function_interpolator;

    AdvectionFieldBuilder const& m_adv_field_builder;
    AdvectionFieldEvaluator const& m_adv_field_evaluator;

    TimeStepper const& m_time_stepper;

public:
    explicit BslAdvection1D(
            FunctionPreallocatableInterpolatorType const& function_interpolator,
            AdvectionFieldBuilder const& adv_field_builder,
            AdvectionFieldEvaluator const& adv_field_evaluator,
            TimeStepper const& time_stepper)
        : m_function_interpolator(function_interpolator)
        , m_adv_field_builder(adv_field_builder)
        , m_adv_field_evaluator(adv_field_evaluator)
        , m_time_stepper(time_stepper)
    {
    }

    ~BslAdvection1D() = default;

    FunctionField operator()(
            FunctionField const allfdistribu,
            AdvecField const advection_field,
            double const dt,
            std::optional<AdvecFieldDerivConstField> const advection_field_derivatives_min
            = std::nullopt,
            std::optional<AdvecFieldDerivConstField> const advection_field_derivatives_max
            = std::nullopt) const
    {
        Kokkos::Profiling::pushRegion("BslAdvection1D");

        // Get index ranges and operators .............................................................
        IdxRangeFunction const idx_range_function = get_idx_range(allfdistribu);
        IdxRangeAdvection const idx_range_advection = get_idx_range(advection_field);

        std::unique_ptr<FunctionInterpolatorType> const function_interpolator_ptr
                = m_function_interpolator.preallocate();
        FunctionInterpolatorType const& function_interpolator = *function_interpolator_ptr;


        // Build spline representation of the advection field ....................................
        AdvecFieldSplineMem advection_field_coefs_alloc(
                m_adv_field_builder.batched_spline_domain());
        AdvecFieldSplineCoeffs advection_field_coefs = get_field(advection_field_coefs_alloc);

        m_adv_field_builder(
                advection_field_coefs,
                get_const_field(advection_field),
                advection_field_derivatives_min,
                advection_field_derivatives_max);

        // Build derivatives on boundaries and fill with zeros....................................
        FunctionDerivFieldMem function_derivatives_min(
                m_function_interpolator.batched_derivs_idx_range_xmin(idx_range_function));
        FunctionDerivFieldMem function_derivatives_max(
                m_function_interpolator.batched_derivs_idx_range_xmax(idx_range_function));
        ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), function_derivatives_min, 0.);
        ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), function_derivatives_max, 0.);

        // Initialise the characteristics on the mesh points .....................................
        /*
            For the time integration solver, the function we advect (here the characteristics)
            need to be defined on the same index range as the advection field. We then work on space
            slices of the characteristic feet.  
        */
        FeetFieldMem slice_feet_alloc(idx_range_advection);
        FeetField slice_feet = get_field(slice_feet_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_advection,
                KOKKOS_LAMBDA(IdxAdvection const idx) {
                    slice_feet(idx) = ddc::coordinate(IdxInterest(idx));
                });


        // Compute the characteristic feet .......................................................
        /*
            We use a time stepper method to solve the charateristic equation. 
            A TimeStepper needs a function to compute the updated advection field and a function to 
            compute the updated feet. 
                * update_adv_field: evaluate the advection field spline at the updated feet. 
        */
        // The function describing how the derivative of the evolve function is calculated.
        std::function<void(AdvecField, FeetConstField)> update_adv_field
                = [&](AdvecField updated_advection_field, FeetConstField slice_feet) {
                      m_adv_field_evaluator(
                              updated_advection_field,
                              slice_feet,
                              get_const_field(advection_field_coefs));
                  };


        // Solve the characteristic equation with a time integration method
        m_time_stepper
                .update(Kokkos::DefaultExecutionSpace(),
                        get_field(slice_feet),
                        -dt,
                        update_adv_field);


        // Interpolate the function ..............................................................
        /*
            To interpolate the function we want to advect, we build for the feet a Field defined 
            on the index range where the function is defined. 
        */
        FieldMem<CoordInterest, IdxRangeFunction> feet_alloc(idx_range_function);
        Field<CoordInterest, IdxRangeFunction> feet = get_field(feet_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_function,
                KOKKOS_LAMBDA(IdxFunction const idx) {
                    IdxAdvection slice_foot_index(idx);
                    feet(idx) = slice_feet(slice_foot_index);
                });


        // Interpolate the function at the characteristic feet
        function_interpolator(
                allfdistribu,
                get_const_field(feet),
                get_const_field(function_derivatives_min),
                get_const_field(function_derivatives_max));


        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
```


