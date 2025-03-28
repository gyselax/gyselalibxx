

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
#include "bsl_advection_rtheta.hpp"
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
    using EulerMethod
            = Euler<FieldMemRTheta<CoordRTheta>,
                    DVectorFieldMemRTheta<X, Y>,
                    Kokkos::DefaultExecutionSpace>;

    using EulerMethod_host
            = Euler<host_t<FieldMemRTheta<CoordRTheta>>,
                    host_t<DVectorFieldMemRTheta<X, Y>>,
                    Kokkos::DefaultHostExecutionSpace>;

    using SplinePolarFootFinderType = SplinePolarFootFinder<
            EulerMethod,
            LogicalToPhysicalMapping,
            LogicalToPseudoPhysicalMapping,
            SplineRThetaBuilder,
            SplineRThetaEvaluatorConstBound>;

    using SplinePolarFootFinderType_host = SplinePolarFootFinder<
            EulerMethod_host,
            LogicalToPhysicalMapping,
            LogicalToPseudoPhysicalMapping,
            SplineRThetaBuilder_host,
            SplineRThetaEvaluatorConstBound_host>;


    LogicalToPhysicalMapping const& m_logical_to_physical;

    BslAdvectionRTheta<SplinePolarFootFinderType, LogicalToPhysicalMapping> const&
            m_advection_solver;

    EulerMethod_host const m_euler;
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
            BslAdvectionRTheta<SplinePolarFootFinderType, LogicalToPhysicalMapping>&
                    advection_solver,
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
        , m_euler(grid)
        , m_find_feet(
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
            host_t<DFieldRTheta> allfdistribu_host,
            double const dt,
            int const steps) const final
    {
        std::chrono::time_point<std::chrono::system_clock> start_time;
        std::chrono::time_point<std::chrono::system_clock> end_time;

        // Grid. ------------------------------------------------------------------------------------------
        IdxRangeRTheta const grid(get_idx_range<GridR, GridTheta>(allfdistribu_host));

        host_t<FieldMemRTheta<CoordRTheta>> coords(grid);
        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            coords(irtheta) = ddc::coordinate(irtheta);
        });

        IdxRangeBSR radial_bsplines(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                IdxStep<BSplinesR> {PolarBSplinesRTheta::continuity + 1}));
        IdxRangeBSTheta polar_idx_range(ddc::discrete_space<BSplinesTheta>().full_domain());

        // --- Electrostatic potential (phi). -------------------------------------------------------------
        DFieldMemRTheta electrical_potential(grid);
        host_t<DFieldMemRTheta> electrical_potential_host(grid);

        host_t<PolarSplineMemRTheta> electrostatic_potential_coef(
                PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>(),
                IdxRangeBSRTheta(radial_bsplines, polar_idx_range));

        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule>
                polar_spline_evaluator(extrapolation_rule);

        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        host_t<DVectorFieldMemRTheta<X, Y>> electric_field_alloc(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> electric_field_predicted_alloc(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_alloc(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_predicted_alloc(grid);

        host_t<DVectorFieldRTheta<X, Y>> electric_field(electric_field_alloc);
        host_t<DVectorFieldRTheta<X, Y>> electric_field_predicted(electric_field_predicted_alloc);
        host_t<DVectorFieldRTheta<X, Y>> advection_field_host(advection_field_alloc);
        host_t<DVectorFieldRTheta<X, Y>> advection_field_predicted(advection_field_predicted_alloc);

        AdvectionFieldFinder advection_field_computer(m_logical_to_physical);



        // --- Parameter for linearisation of advection field: --------------------------------------------
        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            double const time = iter * dt;
            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            host_t<Spline2DMem> allfdistribu_coef(get_spline_idx_range(m_builder));
            m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu_host));
            PoissonLikeRHSFunction const
                    charge_density_coord_1(get_const_field(allfdistribu_coef), m_evaluator);
            m_poisson_solver(charge_density_coord_1, electrostatic_potential_coef);

            polar_spline_evaluator(
                    get_field(electrical_potential_host),
                    get_const_field(coords),
                    get_const_field(electrostatic_potential_coef));

            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .with("time", time)
                    .with("density", allfdistribu_host)
                    .with("electrical_potential", electrical_potential_host);

            // STEP 2: From phi^n, we compute A^n:
            advection_field_computer(electrostatic_potential_coef, advection_field_host);


            // STEP 3: From rho^n and A^n, we compute rho^P: Vlasov equation
            // --- Copy rho^n because it will be modified:
            DFieldMemRTheta allfdistribu_predicted(grid);
            ddc::parallel_deepcopy(get_field(allfdistribu_predicted), allfdistribu_host);
            auto advection_field = ddcHelper::create_mirror_view_and_copy(
                    Kokkos::DefaultExecutionSpace(),
                    advection_field_host);
            m_advection_solver(get_field(allfdistribu_predicted), get_field(advection_field), dt);

            // --- advect also the feet because it is needed for the next step
            host_t<FieldMemRTheta<CoordRTheta>> feet_coords(grid);
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                feet_coords(irtheta) = CoordRTheta(ddc::coordinate(irtheta));
            });
            m_find_feet(get_field(feet_coords), get_field(advection_field_host), dt);

            auto allfdistribu_predicted_host
                    = ddc::create_mirror_view_and_copy(get_field(allfdistribu_predicted));
            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu_predicted_host));
            PoissonLikeRHSFunction const
                    charge_density_coord_4(get_const_field(allfdistribu_coef), m_evaluator);
            m_poisson_solver(charge_density_coord_4, electrostatic_potential_coef);

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(electrostatic_potential_coef, advection_field_predicted);


            // ---  we evaluate the advection field A^n at the characteristic feet X^P
            host_t<DVectorFieldMemRTheta<X, Y>> advection_field_evaluated(grid);
            host_t<VectorSplineCoeffsMem2D<X, Y>> advection_field_coefs(
                    get_spline_idx_range(m_builder));

            m_builder(
                    ddcHelper::get<X>(advection_field_coefs),
                    ddcHelper::get<X>(get_const_field(advection_field_host)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs),
                    ddcHelper::get<Y>(get_const_field(advection_field_host)));

            m_evaluator(
                    get_field(ddcHelper::get<X>(advection_field_evaluated)),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs)));
            m_evaluator(
                    get_field(ddcHelper::get<Y>(advection_field_evaluated)),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs)));


            // STEP 6: From rho^n and (A^n(X^P) + A^P(X^n))/2, we compute rho^{n+1}: Vlasov equation
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                ddcHelper::get<X>(advection_field_host)(irtheta)
                        = (ddcHelper::get<X>(advection_field_evaluated)(irtheta)
                           + ddcHelper::get<X>(advection_field_predicted)(irtheta))
                          / 2.;
                ddcHelper::get<Y>(advection_field_host)(irtheta)
                        = (ddcHelper::get<Y>(advection_field_evaluated)(irtheta)
                           + ddcHelper::get<Y>(advection_field_predicted)(irtheta))
                          / 2.;
            });

            ddc::parallel_deepcopy(
                    ddcHelper::get<X>(advection_field),
                    ddcHelper::get<X>(advection_field_host));
            ddc::parallel_deepcopy(
                    ddcHelper::get<Y>(advection_field),
                    ddcHelper::get<Y>(advection_field_host));
            auto allfdistribu = ddc::
                    create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), allfdistribu_host);
            m_advection_solver(get_field(allfdistribu), get_const_field(advection_field), dt);
            ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        host_t<Spline2DMem> allfdistribu_coef(get_spline_idx_range(m_builder));
        m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu_host));
        PoissonLikeRHSFunction const
                charge_density_coord(get_const_field(allfdistribu_coef), m_evaluator);
        m_poisson_solver(charge_density_coord, get_field(electrical_potential));
        ddc::parallel_deepcopy(electrical_potential_host, electrical_potential);
        ddc::PdiEvent("last_iteration")
                .with("iter", steps)
                .with("time", steps * dt)
                .with("density", allfdistribu_host)
                .with("electrical_potential", electrical_potential_host);


        end_time = std::chrono::system_clock::now();
        display_time_difference("Iterations time: ", start_time, end_time);


        return allfdistribu_host;
    }
};
```


