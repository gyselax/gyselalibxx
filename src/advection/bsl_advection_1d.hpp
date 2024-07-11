// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include "ddc_helper.hpp"
#include "euler.hpp"
#include "iinterpolator.hpp"


/**
 * @brief A class which computes the advection along the dimension of interest IDimInterest. 
 * 
 * This operator solves the following equation type
 * 
 * @f$ \partial_t f_s(t,x) + A_{s, x_i} (x') \cdot \partial_{x_i} f_s (t, x) = 0, \qquad x\in \Omega, x'\in\Omega'@f$
 * 
 * with 
 *  * @f$ f @f$, a function defined on a domain @f$ \Omega @f$; 
 *  * @f$ A @f$, an advection field defined on subdomain @f$ \Omega'\subset \Omega @f$; 
 *  * @f$ x_i @f$, an advection dimension.
 * 
 * 
 * The characteristic equation is solved on the advection domain @f$ \Omega'@f$. 
 * Then the feet on @f$ \Omega@f$ are computed from the characteristic feet on @f$ \Omega'@f$ and the function 
 * @f$ f @f$ is interpolated at the feet in @f$ \Omega @f$. 
 * 
 * The characteristic equation is solved using a time integration method (ITimeStepper). 
 * 
 * 
 * @tparam IDimInterest
 *          The dimension along which the advection is computed. 
 *          It refers to the dimension of @f$ x_i @f$ in the equation. 
 * @tparam AdvectionDomain
 *          The domain @f$ \Omega' @f$ where the characteristic equation is solved. 
 *          It also refers to the domain of the advection field. 
 *          It had to also be defined on the IDimInterest for the time integration method. 
 * @tparam FunctionDomain
 *          The domain @f$ \Omega @f$ where allfdistribu is defined. 
 * @tparam AdvectionFieldBuilder
 *          The type of the spline builder for the advection field (see SplineBuilder). 
 * @tparam AdvectionFieldEvaluator
 *          The type of the spline evalutor for the advection field (see SplineEvaluator).  
 * @tparam TimeStepper 
 *          The time integration method applied to solve the characteristic equation. 
 *          The method is picked among the child classes of ITimeStepper. 
 */
template <
        class IDimInterest,
        class AdvectionDomain,
        class FunctionDomain,
        class AdvectionFieldBuilder,
        class AdvectionFieldEvaluator,
        class TimeStepper
        = Euler<device_t<ddc::Chunk<
                        ddc::Coordinate<typename IDimInterest::continuous_dimension_type>,
                        AdvectionDomain>>,
                device_t<ddc::Chunk<double, AdvectionDomain>>>>
class BslAdvection1D
{
private:
    // Advection domain element:
    using AdvectionIndex = typename AdvectionDomain::discrete_element_type;

    // Full domain element:
    using FunctionIndex = typename FunctionDomain::discrete_element_type;

    // Advection dimension (or Interest dimension):
    using RDimInterest = typename IDimInterest::continuous_dimension_type;
    using CoordInterest = ddc::Coordinate<RDimInterest>;
    using DomainInterest = ddc::DiscreteDomain<IDimInterest>;
    using IndexInterest = typename DomainInterest::discrete_element_type;

    // Type for the feet and advection field:
    using FeetChunk = device_t<ddc::Chunk<CoordInterest, AdvectionDomain>>;
    using FeetSpan = typename FeetChunk::span_type;
    using FeetView = typename FeetChunk::view_type;

    using AdvecFieldChunk = device_t<ddc::Chunk<double, AdvectionDomain>>;
    using AdvecFieldSpan = typename AdvecFieldChunk::span_type;

    using FunctionSpan = device_t<ddc::ChunkSpan<double, FunctionDomain>>;

    // Type for spline representation of the advection field
    using BSAdvectionDomain = typename AdvectionFieldBuilder::batched_spline_domain_type;
    using AdvecFieldSplineChunk = device_t<ddc::Chunk<double, BSAdvectionDomain>>;
    using AdvecFieldSplineSpan = device_t<ddc::ChunkSpan<double, BSAdvectionDomain>>;

    // Type for the derivatives of the advection field
    using DerivDim = ddc::Deriv<RDimInterest>;
    using DomainInterestDeriv = ddc::DiscreteDomain<DerivDim>;
    using AdvecFieldDerivDomain = decltype(ddc::replace_dim_of<IDimInterest, DerivDim>(
            std::declval<AdvectionDomain>(),
            std::declval<DomainInterestDeriv>()));
    using AdvecFieldDerivView = device_t<ddc::ChunkSpan<const double, AdvecFieldDerivDomain>>;


    // Interpolators:
    using FunctionPreallocatableInterpolatorType
            = interpolator_on_domain_t<IPreallocatableInterpolator, IDimInterest, FunctionDomain>;
    using FunctionInterpolatorType
            = interpolator_on_domain_t<IInterpolator, IDimInterest, FunctionDomain>;


    FunctionPreallocatableInterpolatorType const& m_function_interpolator;

    AdvectionFieldBuilder const& m_adv_field_builder;
    AdvectionFieldEvaluator const& m_adv_field_evaluator;

    TimeStepper const& m_time_stepper;

public:
    /**
     * @brief Constructor when the advection domain and the function domain are different. 
     * 
     * When AdvectionDomain and FunctionDomain are different, we need one interpolator for 
     * each domain. 
     * 
     * We can also use it when we want two differents interpolators but defined on the same 
     * domain (e.g. different boundary conditions for the evaluators).
     * 
     * @param[in] function_interpolator interpolator along the IDimInterest direction to interpolate 
     *          the advected function (allfdistribu) on the domain of the function.
     * @param[in] adv_field_builder builder along the IDimInterest direction to build a spline representation
     *          of the advection field on the function domain. 
     * @param[in] adv_field_evaluator evaluator along the IDimInterest direction to evaluate 
     *          the advection field spline representation on the function domain.  
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
     * @brief Advects allfdistribu along the advection dimension IDimInterest for a duration dt.
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
    FunctionSpan operator()(
            FunctionSpan const allfdistribu,
            AdvecFieldSpan const advection_field,
            double const dt,
            std::optional<AdvecFieldDerivView> const advection_field_derivatives_min = std::nullopt,
            std::optional<AdvecFieldDerivView> const advection_field_derivatives_max
            = std::nullopt) const
    {
        Kokkos::Profiling::pushRegion("BslAdvection1D");

        // Get domains and operators .............................................................
        FunctionDomain const function_dom = allfdistribu.domain();
        AdvectionDomain const advection_dom = advection_field.domain();
        auto const batch_dom = ddc::remove_dims_of(function_dom, advection_dom);
        DomainInterest const interest_dom(advection_dom);

        std::unique_ptr<FunctionInterpolatorType> const function_interpolator_ptr
                = m_function_interpolator.preallocate();
        FunctionInterpolatorType const& function_interpolator = *function_interpolator_ptr;


        // Build spline representation of the advection field ....................................
        AdvecFieldSplineChunk advection_field_coefs_alloc(
                m_adv_field_builder.batched_spline_domain());
        AdvecFieldSplineSpan advection_field_coefs = advection_field_coefs_alloc.span_view();

        m_adv_field_builder(
                advection_field_coefs,
                advection_field.span_cview(),
                advection_field_derivatives_min,
                advection_field_derivatives_max);


        // Initialise the characteristics on the mesh points .....................................
        /*
            For the time integration solver, the function we advect (here the characteristics)
            need to be defined on the same domain as the advection field. We then work on space
            slices of the characteristic feet.  
        */
        FeetChunk slice_feet_alloc(advection_dom);
        FeetSpan slice_feet = slice_feet_alloc.span_view();
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                advection_dom,
                KOKKOS_LAMBDA(AdvectionIndex const idx) {
                    slice_feet(idx) = ddc::coordinate(IndexInterest(idx));
                });


        // Compute the characteristic feet .......................................................
        /*
            We use a time stepper method to solve the charateristic equation. 
            A TimeStepper needs a function to compute the updated advection field and a function to 
            compute the updated feet. 
                * update_adv_field: evaluate the advection field spline at the updated feet. 
        */
        // The function describing how the derivative of the evolve function is calculated.
        std::function<void(AdvecFieldSpan, FeetView)> update_adv_field
                = [&](AdvecFieldSpan updated_advection_field, FeetView slice_feet) {
                      m_adv_field_evaluator(
                              updated_advection_field,
                              slice_feet,
                              advection_field_coefs.span_cview());
                  };


        // Solve the characteristic equation with a time integration method
        m_time_stepper
                .update(Kokkos::DefaultExecutionSpace(),
                        slice_feet.span_view(),
                        -dt,
                        update_adv_field);


        // Interpolate the function ..............................................................
        /*
            To interpolate the function we want to advect, we build for the feet a Chunk defined 
            on the domain where the function is defined. 
        */
        device_t<ddc::Chunk<CoordInterest, FunctionDomain>> feet_alloc(function_dom);
        ddc::ChunkSpan feet = feet_alloc.span_view();
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                function_dom,
                KOKKOS_LAMBDA(FunctionIndex const idx) {
                    AdvectionIndex slice_foot_index(idx);
                    feet(idx) = slice_feet(slice_foot_index);
                });


        // Interpolate the function at the characteristic feet
        function_interpolator(allfdistribu, feet);


        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
