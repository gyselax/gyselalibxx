

# File bsl\_advection\_1d.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**bsl\_advection\_1d.hpp**](bsl__advection__1d_8hpp.md)

[Go to the documentation of this file](bsl__advection__1d_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "euler.hpp"
#include "i_interpolation.hpp"
#include "itimestepper.hpp"


template <
        class GridInterest,
        class IdxRangeAdvection,
        class IdxRangeFunction,
        concepts::Interpolation FunctionInterpolator,
        concepts::Interpolation AdvectionFieldInterpolator,
        class TimeStepperBuilder = EulerBuilder,
        class DataType = double>
class BslAdvection1D
{
    static_assert(is_timestepper_builder_v<TimeStepperBuilder>);
    static_assert(std::is_floating_point_v<DataType>);

public:
    using FunctionBuilder = typename FunctionInterpolator::BuilderType;
    using FunctionEvaluator = typename FunctionInterpolator::EvaluatorType;
    using AdvectionFieldBuilder = typename AdvectionFieldInterpolator::BuilderType;
    using AdvectionFieldEvaluator = typename AdvectionFieldInterpolator::EvaluatorType;

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

    using AdvecFieldMem = FieldMem<DataType, IdxRangeAdvection>;
    using AdvecField = typename AdvecFieldMem::span_type;

    using FunctionField = Field<DataType, IdxRangeFunction>;

    // Type for spline representation of the advection field
    using IdxRangeBSAdvection = typename InterpolationBuilderTraits<
            AdvectionFieldBuilder>::template batched_basis_idx_range_type<IdxRangeAdvection>;
    using AdvecFieldSplineMem = FieldMem<DataType, IdxRangeBSAdvection>;
    using AdvecFieldSplineCoeffs = Field<DataType, IdxRangeBSAdvection>;

    // Type for the derivatives of the advection field
    using DerivDim = ddc::Deriv<DimInterest>;
    using IdxRangeAdvecFieldDeriv
            = ddc::replace_dim_of_t<IdxRangeAdvection, GridInterest, DerivDim>;
    using AdvecFieldDerivConstField = Field<const DataType, IdxRangeAdvecFieldDeriv>;

    // Type for the spline representation of the function
    using IdxRangeFunctionBasis = typename InterpolationBuilderTraits<
            FunctionBuilder>::template batched_basis_idx_range_type<IdxRangeFunction>;
    using FunctionBasisFieldMem = FieldMem<DataType, IdxRangeFunctionBasis>;

    // Type for the derivatives of the function
    using IdxRangeFunctionDeriv = typename InterpolationBuilderTraits<
            FunctionBuilder>::template batched_derivs_idx_range_type<IdxRangeFunction>;
    using FunctionDerivFieldMem = FieldMem<DataType, IdxRangeFunctionDeriv>;

    using TimeStepper = typename TimeStepperBuilder::template time_stepper_t<
            FieldMem<Coord<typename GridInterest::continuous_dimension_type>, IdxRangeAdvection>,
            FieldMem<DataType, IdxRangeAdvection>>;

    FunctionBuilder const& m_function_builder;
    FunctionEvaluator const& m_function_evaluator;

    AdvectionFieldBuilder const& m_adv_field_builder;
    AdvectionFieldEvaluator const& m_adv_field_evaluator;

    TimeStepperBuilder const& m_time_stepper_builder;

public:
    [[deprecated]] explicit BslAdvection1D(
            FunctionBuilder const& function_builder,
            FunctionEvaluator const& function_evaluator,
            AdvectionFieldBuilder const& adv_field_builder,
            AdvectionFieldEvaluator const& adv_field_evaluator,
            TimeStepperBuilder const& time_stepper_builder)
        : m_function_builder(function_builder)
        , m_function_evaluator(function_evaluator)
        , m_adv_field_builder(adv_field_builder)
        , m_adv_field_evaluator(adv_field_evaluator)
        , m_time_stepper_builder(time_stepper_builder)
    {
    }

    explicit BslAdvection1D(
            FunctionInterpolator const& function_interpolator,
            AdvectionFieldInterpolator const& adv_field_interpolator,
            TimeStepperBuilder const& time_stepper_builder)
        : m_function_builder(function_interpolator.get_builder())
        , m_function_evaluator(function_interpolator.get_evaluator())
        , m_adv_field_builder(adv_field_interpolator.get_builder())
        , m_adv_field_evaluator(adv_field_interpolator.get_evaluator())
        , m_time_stepper_builder(time_stepper_builder)
    {
    }

    explicit BslAdvection1D(
            FunctionInterpolator const& interpolator,
            TimeStepperBuilder const& time_stepper_builder)
        : m_function_builder(interpolator.get_builder())
        , m_function_evaluator(interpolator.get_evaluator())
        , m_adv_field_builder(interpolator.get_builder())
        , m_adv_field_evaluator(interpolator.get_evaluator())
        , m_time_stepper_builder(time_stepper_builder)
    {
        static_assert(std::is_same_v<FunctionInterpolator, AdvectionFieldInterpolator>);
    }

    ~BslAdvection1D() = default;

    FunctionField operator()(
            FunctionField const allfdistribu,
            AdvecField const advection_field,
            DataType const dt,
            std::optional<AdvecFieldDerivConstField> const advection_field_derivatives_min
            = std::nullopt,
            std::optional<AdvecFieldDerivConstField> const advection_field_derivatives_max
            = std::nullopt) const
    {
        Kokkos::Profiling::pushRegion("BslAdvection1D");

        // Get index ranges and operators ........................................................
        IdxRangeFunction const idx_range_function = get_idx_range(allfdistribu);
        IdxRangeAdvection const idx_range_advection = get_idx_range(advection_field);


        // Build spline representation of the advection field ....................................
        AdvecFieldSplineMem advection_field_coefs_alloc(
                "advection_field_coefs (BslAdvection1D::operator())",
                batched_basis_idx_range(m_adv_field_builder, get_idx_range(advection_field)));
        AdvecFieldSplineCoeffs advection_field_coefs = get_field(advection_field_coefs_alloc);

        m_adv_field_builder(
                advection_field_coefs,
                get_const_field(advection_field),
                advection_field_derivatives_min,
                advection_field_derivatives_max);

        // Allocate buffer for the function interpolation coefficients ...........................
        FunctionBasisFieldMem function_coefs_alloc(
                "function_coefs (BslAdvection1D::operator())",
                batched_basis_idx_range(m_function_builder, idx_range_function));

        // Build derivatives on boundaries and fill with zeros....................................
        FunctionDerivFieldMem function_derivatives_min(
                m_function_builder.batched_derivs_xmin_domain(idx_range_function));
        FunctionDerivFieldMem function_derivatives_max(
                m_function_builder.batched_derivs_xmax_domain(idx_range_function));
        ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), function_derivatives_min, 0.);
        ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), function_derivatives_max, 0.);

        // Initialise the characteristics on the mesh points .....................................
        /*
            For the time integration solver, the function we advect (here the characteristics)
            need to be defined on the same index range as the advection field. We then work on space
            slices of the characteristic feet.
        */
        FeetFieldMem
                slice_feet_alloc("slice_feet (BslAdvection1D::operator())", idx_range_advection);
        FeetField slice_feet = get_field(slice_feet_alloc);
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                Kokkos::DefaultExecutionSpace(),
                idx_range_advection,
                KOKKOS_LAMBDA(IdxAdvection const idx) {
                    slice_feet(idx) = ddc::coordinate(IdxInterest(idx));
                });


        // Compute the characteristic feet .......................................................
        /*
            We use a time stepper method to solve the characteristic equation.
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

        TimeStepper time_stepper
                = m_time_stepper_builder.template preallocate<TimeStepper>(idx_range_advection);

        // Solve the characteristic equation with a time integration method
        time_stepper
                .update(Kokkos::DefaultExecutionSpace(),
                        get_field(slice_feet),
                        -dt,
                        update_adv_field);


        // Interpolate the function ..............................................................
        /*
            To interpolate the function we want to advect, we build for the feet a Field defined
            on the index range where the function is defined.
        */
        FieldMem<CoordInterest, IdxRangeFunction>
                feet_alloc("feet (BslAdvection1D::operator())", idx_range_function);
        Field<CoordInterest, IdxRangeFunction> feet = get_field(feet_alloc);
        ddc::parallel_for_each(
                location.function_name(),
                Kokkos::DefaultExecutionSpace(),
                idx_range_function,
                KOKKOS_LAMBDA(IdxFunction const idx) {
                    IdxAdvection slice_foot_index(idx);
                    feet(idx) = slice_feet(slice_foot_index);
                });


        // Build interpolation coefficients from the function values
        m_function_builder(
                get_field(function_coefs_alloc),
                get_const_field(allfdistribu),
                std::optional(get_const_field(function_derivatives_min)),
                std::optional(get_const_field(function_derivatives_max)));

        // Evaluate the function at the characteristic feet
        m_function_evaluator(
                allfdistribu,
                get_const_field(feet),
                get_const_field(function_coefs_alloc));


        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
```


