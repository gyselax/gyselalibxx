// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "euler.hpp"
#include "iinterpolator.hpp"


/**
 * @brief A class which computes the advection along the dimension of interest GridInterest. 
 * 
 * This operator solves the following equation type
 * 
 * @f$ \partial_t f_s(t,x) + A_{s, x_i} (x') \cdot \partial_{x_i} f_s (t, x) = 0, \qquad x\in \Omega, x'\in\Omega'@f$
 * 
 * with 
 *  * @f$ f @f$, a function defined on an index range @f$ \Omega @f$; 
 *  * @f$ A @f$, an advection field defined on subindex range @f$ \Omega'\subset \Omega @f$; 
 *  * @f$ x_i @f$, an advection dimension.
 * 
 * 
 * The characteristic equation is solved on the advection index range @f$ \Omega'@f$. 
 * Then the feet on @f$ \Omega@f$ are computed from the characteristic feet on @f$ \Omega'@f$ and the function 
 * @f$ f @f$ is interpolated at the feet in @f$ \Omega @f$. 
 * 
 * The characteristic equation is solved using a time integration method (ITimeStepper). 
 * 
 * 
 * @tparam GridInterest
 *          The dimension along which the advection is computed. 
 *          It refers to the dimension of @f$ x_i @f$ in the equation. 
 * @tparam AdvectionIdxRange
 *          The index range @f$ \Omega' @f$ where the characteristic equation is solved. 
 *          It also refers to the index range of the advection field. 
 *          It had to also be defined on the GridInterest for the time integration method. 
 * @tparam FunctionIdxRange
 *          The index range @f$ \Omega @f$ where allfdistribu is defined. 
 * @tparam AdvectionFieldBuilder
 *          The type of the spline builder for the advection field (see SplineBuilder). 
 * @tparam AdvectionFieldEvaluator
 *          The type of the spline evalutor for the advection field (see SplineEvaluator).  
 * @tparam TimeStepper 
 *          The time integration method applied to solve the characteristic equation. 
 *          The method is picked among the child classes of ITimeStepper. 
 */
template <
        class GridInterest,
        class AdvectionIdxRange,
        class FunctionIdxRange,
        class AdvectionFieldBuilder,
        class AdvectionFieldEvaluator,
        class TimeStepper
        = Euler<FieldMem<
                        Coord<typename GridInterest::continuous_dimension_type>,
                        AdvectionIdxRange>,
                FieldMem<double, AdvectionIdxRange>>>
class BslAdvection1D
{
private:
    // Advection index range element:
    using AdvectionIdx = typename AdvectionIdxRange::discrete_element_type;

    // Full index range element:
    using FunctionIdx = typename FunctionIdxRange::discrete_element_type;

    // Advection dimension (or Interest dimension):
    using DimInterest = typename GridInterest::continuous_dimension_type;
    using CoordInterest = Coord<DimInterest>;
    using IdxRangeInterest = IdxRange<GridInterest>;
    using IdxInterest = typename IdxRangeInterest::discrete_element_type;

    // Type for the feet and advection field:
    using FeetFieldMem = FieldMem<CoordInterest, AdvectionIdxRange>;
    using FeetField = typename FeetFieldMem::span_type;
    using FeetConstField = typename FeetFieldMem::view_type;

    using AdvecFieldMem = FieldMem<double, AdvectionIdxRange>;
    using AdvecField = typename AdvecFieldMem::span_type;

    using FunctionField = Field<double, FunctionIdxRange>;

    // Type for spline representation of the advection field
    using BSAdvectionIdxRange = typename AdvectionFieldBuilder::batched_spline_domain_type;
    using AdvecFieldSplineMem = FieldMem<double, BSAdvectionIdxRange>;
    using AdvecFieldSplineCoeffs = Field<double, BSAdvectionIdxRange>;

    // Type for the derivatives of the advection field
    using DerivDim = ddc::Deriv<DimInterest>;
    using AdvecFieldDerivIdxRange
            = ddc::replace_dim_of_t<AdvectionIdxRange, GridInterest, DerivDim>;
    using AdvecFieldDerivConstField = Field<const double, AdvecFieldDerivIdxRange>;


    // Interpolators:
    using FunctionPreallocatableInterpolatorType = interpolator_on_idx_range_t<
            IPreallocatableInterpolator,
            GridInterest,
            FunctionIdxRange>;
    using FunctionInterpolatorType
            = interpolator_on_idx_range_t<IInterpolator, GridInterest, FunctionIdxRange>;

    // Type for the derivatives of the function
    using FunctionDerivIdxRange = typename FunctionInterpolatorType::batched_derivs_idx_range_type;
    using FunctionDerivFieldMem = FieldMem<double, FunctionDerivIdxRange>;

    FunctionPreallocatableInterpolatorType const& m_function_interpolator;

    AdvectionFieldBuilder const& m_adv_field_builder;
    AdvectionFieldEvaluator const& m_adv_field_evaluator;

    TimeStepper const& m_time_stepper;

public:
    /**
     * @brief Constructor when the advection index range and the function index range are different. 
     * 
     * When AdvectionIdxRange and FunctionIdxRange are different, we need one interpolator for 
     * each index range. 
     * 
     * We can also use it when we want two differents interpolators but defined on the same 
     * index range (e.g. different boundary conditions for the evaluators).
     * 
     * @param[in] function_interpolator interpolator along the GridInterest direction to interpolate 
     *          the advected function (allfdistribu) on the index range of the function.
     * @param[in] adv_field_builder builder along the GridInterest direction to build a spline representation
     *          of the advection field on the function index range. 
     * @param[in] adv_field_evaluator evaluator along the GridInterest direction to evaluate 
     *          the advection field spline representation on the function index range.  
     * @param[in] time_stepper time integration method for the characteristic equation. 
     */
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

    /**
     * @brief Advects allfdistribu along the advection dimension GridInterest for a duration dt.
     * 
     * @param[in, out] allfdistribu Reference to the advected function, allocated on the device 
     * @param[in] advection_field Reference to the advection field, allocated on the device.
     * @param[in] advection_field_derivatives_min Reference to the advection field 
     *              derivatives at the left side of the interest dimension, allocated on the device.
     * @param[in] advection_field_derivatives_max Reference to the advection field 
     *              derivatives at the right side of the interest dimension, allocated on the device.
     * @param[in] dt Time step.
     * 
     * @return A reference to the allfdistribu array after advection on dt. 
     */
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
        FunctionIdxRange const function_dom = get_idx_range(allfdistribu);
        AdvectionIdxRange const advection_dom = get_idx_range(advection_field);
        auto const batch_dom = ddc::remove_dims_of(function_dom, advection_dom);
        IdxRangeInterest const interest_dom(advection_dom);

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
                m_function_interpolator.batched_derivs_idx_range_xmin(function_dom));
        FunctionDerivFieldMem function_derivatives_max(
                m_function_interpolator.batched_derivs_idx_range_xmax(function_dom));
        ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), function_derivatives_min, 0.);
        ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), function_derivatives_max, 0.);

        // Initialise the characteristics on the mesh points .....................................
        /*
            For the time integration solver, the function we advect (here the characteristics)
            need to be defined on the same index range as the advection field. We then work on space
            slices of the characteristic feet.  
        */
        FeetFieldMem slice_feet_alloc(advection_dom);
        FeetField slice_feet = get_field(slice_feet_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                advection_dom,
                KOKKOS_LAMBDA(AdvectionIdx const idx) {
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
            To interpolate the function we want to advect, we build for the feet a Chunk defined 
            on the index range where the function is defined. 
        */
        FieldMem<CoordInterest, FunctionIdxRange> feet_alloc(function_dom);
        ddc::ChunkSpan feet = get_field(feet_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                function_dom,
                KOKKOS_LAMBDA(FunctionIdx const idx) {
                    AdvectionIdx slice_foot_index(idx);
                    feet(idx) = slice_feet(slice_foot_index);
                });


        // Interpolate the function at the characteristic feet
        function_interpolator(
                allfdistribu,
                feet,
                function_derivatives_min,
                function_derivatives_max);


        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
