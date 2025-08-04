

# File bsl\_predcorr\_second\_order\_explicit.hpp

[**File List**](files.md) **>** [**geometryRTheta**](dir_e9f169004bcfe9f3cb1f8a27ce024e59.md) **>** [**time\_solver**](dir_4c2664fc2adc717d896afdb0f76e6fe5.md) **>** [**bsl\_predcorr\_second\_order\_explicit.hpp**](bsl__predcorr__second__order__explicit_8hpp.md)

[Go to the documentation of this file](bsl__predcorr__second__order__explicit_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "advection_field_rtheta.hpp"
#include "bsl_advection_polar.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "itimesolver.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "spline_polar_foot_finder.hpp"



template <class LogicalToPhysicalMapping, class LogicalToPseudoPhysicalMapping>
class BslExplicitPredCorrRTheta : public ITimeSolverRTheta
{
private:
    using SplinePolarFootFinderType = SplinePolarFootFinder<
            IdxRangeRTheta,
            EulerBuilder,
            LogicalToPhysicalMapping,
            LogicalToPseudoPhysicalMapping,
            SplineRThetaBuilder,
            SplineRThetaEvaluatorConstBound>;

    using SplinePolarFootFinderType_host = SplinePolarFootFinder<
            IdxRangeRTheta,
            EulerBuilder,
            LogicalToPhysicalMapping,
            LogicalToPseudoPhysicalMapping,
            SplineRThetaBuilder_host,
            SplineRThetaEvaluatorConstBound_host>;

    using BslAdvectionRTheta = BslAdvectionPolar<
            SplinePolarFootFinderType,
            LogicalToPhysicalMapping,
            PreallocatableSplineInterpolator2D<
                    SplineRThetaBuilder,
                    SplineRThetaEvaluatorNullBound,
                    IdxRangeRTheta>>;


    LogicalToPhysicalMapping const& m_logical_to_physical;

    BslAdvectionRTheta const& m_advection_solver;

    EulerBuilder const m_euler;
    SplinePolarFootFinderType_host const m_find_feet;

    PolarSplineFEMPoissonLikeSolver<
            GridR,
            GridTheta,
            PolarBSplinesRTheta,
            SplineRThetaEvaluatorNullBound> const& m_poisson_solver;

    SplineRThetaBuilder_host const& m_builder;
    SplineRThetaEvaluatorConstBound_host const& m_evaluator;



public:
    BslExplicitPredCorrRTheta(
            LogicalToPhysicalMapping const& logical_to_physical,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical,
            BslAdvectionRTheta const& advection_solver,
            IdxRangeRTheta const& grid,
            SplineRThetaBuilder_host const& builder,
            PolarSplineFEMPoissonLikeSolver<
                    GridR,
                    GridTheta,
                    PolarBSplinesRTheta,
                    SplineRThetaEvaluatorNullBound> const& poisson_solver,
            SplineRThetaEvaluatorConstBound_host const& advection_evaluator)
        : m_logical_to_physical(logical_to_physical)
        , m_advection_solver(advection_solver)
        , m_find_feet(
                  grid,
                  m_euler,
                  logical_to_physical,
                  logical_to_pseudo_physical,
                  builder,
                  advection_evaluator)
        , m_poisson_solver(poisson_solver)
        , m_builder(builder)
        , m_evaluator(advection_evaluator)
    {
    }



    host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> density_host,
            double const dt,
            int const steps) const final
    {
        std::chrono::time_point<std::chrono::system_clock> start_time;
        std::chrono::time_point<std::chrono::system_clock> end_time;

        // Grid. ------------------------------------------------------------------------------------------
        IdxRangeRTheta const grid(get_idx_range<GridR, GridTheta>(density_host));

        // --- Electrostatic potential (phi). -------------------------------------------------------------
        DFieldMemRTheta electrical_potential(grid);
        host_t<DFieldMemRTheta> electrical_potential_host(grid);

        host_t<PolarSplineMemRTheta> electrostatic_potential_coef_alloc_host(
                ddc::discrete_space<PolarBSplinesRTheta>().full_domain());

        host_t<Spline2DMem> density_coef_alloc_host(get_spline_idx_range(m_builder));
        DFieldMemRTheta density_predicted_alloc(grid);
        auto density_predicted_host_alloc
                = ddc::create_mirror_view(get_field(density_predicted_alloc));
        auto density_alloc = ddc::create_mirror_view(Kokkos::DefaultExecutionSpace(), density_host);
        host_t<FieldMemRTheta<CoordRTheta>> feet_coords_alloc_host(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_evaluated_alloc_host(grid);
        host_t<VectorSplineCoeffsMem2D<X, Y>> advection_field_coefs_alloc_host(
                get_spline_idx_range(m_builder));


        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_alloc_host(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_predicted_alloc_host(grid);
        auto advection_field_alloc = ddcHelper::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(advection_field_alloc_host));

        // Create fields
        host_t<DVectorFieldRTheta<X, Y>> advection_field_host(advection_field_alloc_host);
        host_t<DVectorFieldRTheta<X, Y>> advection_field_predicted_host(
                advection_field_predicted_alloc_host);
        DVectorFieldRTheta<X, Y> advection_field(advection_field_alloc);
        host_t<DFieldRTheta> density_predicted_host(density_predicted_host_alloc);
        DFieldRTheta density = get_field(density_alloc);

        // --- Operators ----------------------------------------------------------------------------------
        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<
                Kokkos::DefaultHostExecutionSpace,
                Kokkos::HostSpace,
                PolarBSplinesRTheta,
                ddc::NullExtrapolationRule>
                polar_spline_evaluator(extrapolation_rule);

        AdvectionFieldFinder advection_field_computer(m_logical_to_physical);

        PoissonLikeRHSFunction const
                charge_density(get_const_field(density_coef_alloc_host), m_evaluator);


        // --- Parameter for linearisation of advection field: --------------------------------------------
        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            double const time = iter * dt;
            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            m_builder(get_field(density_coef_alloc_host), get_const_field(density_host));
            m_poisson_solver(charge_density, get_field(electrostatic_potential_coef_alloc_host));

            polar_spline_evaluator(
                    get_field(electrical_potential_host),
                    get_const_field(electrostatic_potential_coef_alloc_host));

            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .with("time", time)
                    .with("density", density_host)
                    .with("electrical_potential", electrical_potential_host);

            // STEP 2: From phi^n, we compute A^n:
            advection_field_computer(
                    get_field(electrostatic_potential_coef_alloc_host),
                    advection_field_host);

            // STEP 3: From rho^n and A^n, we compute rho^P: Vlasov equation
            // --- Copy rho^n because it will be modified:
            ddc::parallel_deepcopy(get_field(density_predicted_alloc), density_host);
            ddcHelper::deepcopy(advection_field, advection_field_host);
            m_advection_solver(get_field(density_predicted_alloc), advection_field, dt);

            // --- advect also the feet because it is needed for the next step
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                feet_coords_alloc_host(irtheta) = CoordRTheta(ddc::coordinate(irtheta));
            });
            m_find_feet(get_field(feet_coords_alloc_host), get_field(advection_field_host), dt);

            ddc::parallel_deepcopy(density_predicted_host, get_field(density_predicted_alloc));
            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            m_builder(get_field(density_coef_alloc_host), get_const_field(density_predicted_host));
            m_poisson_solver(charge_density, get_field(electrostatic_potential_coef_alloc_host));

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(
                    get_field(electrostatic_potential_coef_alloc_host),
                    advection_field_predicted_host);


            // ---  we evaluate the advection field A^n at the characteristic feet X^P
            m_builder(
                    ddcHelper::get<X>(advection_field_coefs_alloc_host),
                    ddcHelper::get<X>(get_const_field(advection_field_host)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs_alloc_host),
                    ddcHelper::get<Y>(get_const_field(advection_field_host)));

            m_evaluator(
                    get_field(ddcHelper::get<X>(advection_field_evaluated_alloc_host)),
                    get_const_field(feet_coords_alloc_host),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs_alloc_host)));
            m_evaluator(
                    get_field(ddcHelper::get<Y>(advection_field_evaluated_alloc_host)),
                    get_const_field(feet_coords_alloc_host),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs_alloc_host)));


            // STEP 6: From rho^n and (A^n(X^P) + A^P(X^n))/2, we compute rho^{n+1}: Vlasov equation
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                ddcHelper::assign_vector_field_element(
                        advection_field_host,
                        irtheta,
                        (advection_field_evaluated_alloc_host(irtheta)
                         + advection_field_predicted_host(irtheta))
                                / 2.);
            });

            ddcHelper::deepcopy(advection_field, advection_field_host);
            ddc::parallel_deepcopy(density, density_host);
            m_advection_solver(get_field(density), get_const_field(advection_field), dt);
            ddc::parallel_deepcopy(density_host, density);
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        m_builder(get_field(density_coef_alloc_host), get_const_field(density_host));
        m_poisson_solver(charge_density, get_field(electrical_potential));
        ddc::parallel_deepcopy(electrical_potential_host, electrical_potential);
        ddc::PdiEvent("last_iteration")
                .with("iter", steps)
                .with("time", steps * dt)
                .with("density", density_host)
                .with("electrical_potential", electrical_potential_host);


        end_time = std::chrono::system_clock::now();
        display_time_difference("Iterations time: ", start_time, end_time);


        return density_host;
    }
};
```


