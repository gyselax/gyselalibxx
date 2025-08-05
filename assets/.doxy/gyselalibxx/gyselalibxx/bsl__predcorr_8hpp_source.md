

# File bsl\_predcorr.hpp

[**File List**](files.md) **>** [**geometryRTheta**](dir_e9f169004bcfe9f3cb1f8a27ce024e59.md) **>** [**time\_solver**](dir_4c2664fc2adc717d896afdb0f76e6fe5.md) **>** [**bsl\_predcorr.hpp**](bsl__predcorr_8hpp.md)

[Go to the documentation of this file](bsl__predcorr_8hpp.md)


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
#include "geometry.hpp"
#include "itimesolver.hpp"
#include "l_norm_tools.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "rk2.hpp"

template <class Mapping, class FootFinder>
class BslPredCorrRTheta : public ITimeSolverRTheta
{
    using BslAdvectionRTheta = BslAdvectionPolar<
            FootFinder,
            Mapping,
            PreallocatableSplineInterpolator2D<
                    SplineRThetaBuilder,
                    SplineRThetaEvaluatorNullBound,
                    IdxRangeRTheta>>;

private:
    Mapping const& m_mapping;

    BslAdvectionRTheta const& m_advection_solver;

    PolarSplineFEMPoissonLikeSolver<
            GridR,
            GridTheta,
            PolarBSplinesRTheta,
            SplineRThetaEvaluatorNullBound> const& m_poisson_solver;

    SplineRThetaBuilder_host const& m_builder;
    SplineRThetaEvaluatorNullBound_host const& m_spline_evaluator;


public:
    BslPredCorrRTheta(
            Mapping const& mapping,
            BslAdvectionRTheta const& advection_solver,
            SplineRThetaBuilder_host const& builder,
            SplineRThetaEvaluatorNullBound_host const& rhs_evaluator,
            PolarSplineFEMPoissonLikeSolver<
                    GridR,
                    GridTheta,
                    PolarBSplinesRTheta,
                    SplineRThetaEvaluatorNullBound> const& poisson_solver)
        : m_mapping(mapping)
        , m_advection_solver(advection_solver)
        , m_poisson_solver(poisson_solver)
        , m_builder(builder)
        , m_spline_evaluator(rhs_evaluator)
    {
    }

    host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> density_host,
            double const dt,
            int const steps) const override
    {
        std::chrono::time_point<std::chrono::system_clock> start_time;
        std::chrono::time_point<std::chrono::system_clock> end_time;


        // Grid. ------------------------------------------------------------------------------------------
        IdxRangeRTheta grid(get_idx_range(density_host));

        // Data
        DFieldMemRTheta electrical_potential_alloc(grid);
        host_t<Spline2DMem> density_coef_alloc_host(get_spline_idx_range(m_builder));
        host_t<PolarSplineMemRTheta> electrostatic_potential_coef_alloc_host(
                ddc::discrete_space<PolarBSplinesRTheta>().full_domain());

        auto electrical_potential_alloc_host
                = ddc::create_mirror_view(get_field(electrical_potential_alloc));
        auto density_alloc = ddc::create_mirror_view(Kokkos::DefaultExecutionSpace(), density_host);

        // Fields
        DFieldRTheta electrical_potential = get_field(electrical_potential_alloc);
        host_t<DFieldRTheta> electrical_potential_host = get_field(electrical_potential_alloc_host);
        host_t<PolarSplineRTheta> electrostatic_potential_coef_host
                = get_field(electrostatic_potential_coef_alloc_host);
        DFieldRTheta density = get_field(density_alloc);
        host_t<Spline2D> density_coef_host(density_coef_alloc_host);

        // Operators
        AdvectionFieldFinder advection_field_computer(m_mapping);
        PoissonLikeRHSFunction const
                charge_density(get_const_field(density_coef_host), m_spline_evaluator);

        // Setup
        m_builder(density_coef_host, get_const_field(density_host));
        m_poisson_solver(charge_density, electrical_potential);
        ddc::parallel_deepcopy(electrical_potential_host, get_const_field(electrical_potential));
        ddc::PdiEvent("iteration")
                .with("iter", 0)
                .with("time", 0)
                .with("density", density_host)
                .with("electrical_potential", electrical_potential_host);

        // RK2 methods
        std::function<void(DVectorFieldRTheta<X, Y>, DConstFieldRTheta)> define_advection_field =
                [&](DVectorFieldRTheta<X, Y> advection_field, DConstFieldRTheta density) {
                    ddc::parallel_deepcopy(density_host, density);
                    // --- compute electrostatic potential:
                    m_builder(density_coef_host, get_const_field(density_host));
                    m_poisson_solver(charge_density, electrostatic_potential_coef_host);

                    auto advection_field_alloc_host = ddcHelper::create_mirror_view_and_copy(
                            Kokkos::DefaultHostExecutionSpace(),
                            advection_field);

                    // --- compute advection field:
                    advection_field_computer(
                            electrostatic_potential_coef_host,
                            get_field(advection_field_alloc_host));

                    ddcHelper::
                            deepcopy(advection_field, get_const_field(advection_field_alloc_host));
                };

        RK2<DFieldMemRTheta, DVectorFieldMemRTheta<X, Y>, Kokkos::DefaultExecutionSpace>
                time_stepper(grid);

        // Iteration loop
        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            ddc::parallel_deepcopy(density, get_const_field(density_host));
            time_stepper
                    .update(Kokkos::DefaultExecutionSpace(),
                            density,
                            dt,
                            define_advection_field,
                            m_advection_solver);

            ddc::parallel_deepcopy(density_host, get_const_field(density));
            m_builder(density_coef_host, get_const_field(density_host));
            m_poisson_solver(charge_density, electrical_potential);
            ddc::parallel_deepcopy(
                    electrical_potential_host,
                    get_const_field(electrical_potential));
            ddc::PdiEvent("iteration")
                    .with("iter", iter + 1)
                    .with("time", iter * dt)
                    .with("density", density_host)
                    .with("electrical_potential", electrical_potential_host);
        }
        end_time = std::chrono::system_clock::now();


        display_time_difference("Iterations time: ", start_time, end_time);


        return density_host;
    }
};
```


