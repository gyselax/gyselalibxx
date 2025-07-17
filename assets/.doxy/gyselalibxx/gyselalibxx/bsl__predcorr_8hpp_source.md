

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
        IdxRangeRTheta grid(get_idx_range<GridR, GridTheta>(density_host));
        AdvectionFieldFinder advection_field_computer(m_mapping);

        host_t<PolarSplineMemRTheta> electrostatic_potential_coef(
                ddc::discrete_space<PolarBSplinesRTheta>().full_domain());
        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule>
                polar_spline_evaluator(extrapolation_rule);


        host_t<DFieldMemRTheta> electrical_potential0_host(grid);
        DFieldMemRTheta electrical_potential0(grid);
        host_t<Spline2DMem> density_coef(get_spline_idx_range(m_builder));
        m_builder(get_field(density_coef), get_const_field(density_host));
        PoissonLikeRHSFunction const
                charge_density_coord(get_const_field(density_coef), m_spline_evaluator);
        m_poisson_solver(charge_density_coord, get_field(electrical_potential0));
        ddc::parallel_deepcopy(electrical_potential0, electrical_potential0_host);
        ddc::PdiEvent("iteration")
                .with("iter", 0)
                .with("time", 0)
                .with("density", density_host)
                .with("electrical_potential", electrical_potential0_host);


        std::function<void(host_t<DVectorFieldRTheta<X, Y>>, host_t<DConstFieldRTheta>)>
                define_advection_field = [&](host_t<DVectorFieldRTheta<X, Y>> advection_field_host,
                                             host_t<DConstFieldRTheta> density_host) {
                    // --- compute electrostatic potential:
                    host_t<Spline2DMem> density_coef(get_spline_idx_range(m_builder));
                    m_builder(get_field(density_coef), get_const_field(density_host));
                    m_poisson_solver(charge_density_coord, get_field(electrostatic_potential_coef));

                    // --- compute advection field:
                    advection_field_computer(
                            get_field(electrostatic_potential_coef),
                            advection_field_host);
                };

        std::function<void(host_t<DFieldRTheta>, host_t<DConstVectorFieldRTheta<X, Y>>, double)>
                advect_density = [&](host_t<DFieldRTheta> density_host,
                                     host_t<DConstVectorFieldRTheta<X, Y>> advection_field_host,
                                     double dt) {
                    auto density = ddc::create_mirror_view_and_copy(
                            Kokkos::DefaultExecutionSpace(),
                            density_host);
                    auto advection_field = ddcHelper::create_mirror_view_and_copy(
                            Kokkos::DefaultExecutionSpace(),
                            advection_field_host);
                    m_advection_solver(get_field(density), get_const_field(advection_field), dt);
                    ddc::parallel_deepcopy(density_host, density);
                };

        RK2<host_t<DFieldMemRTheta>,
            host_t<DVectorFieldMemRTheta<X, Y>>,
            Kokkos::DefaultHostExecutionSpace>
                time_stepper(grid);
        DFieldMemRTheta electrical_potential(grid);
        host_t<DFieldMemRTheta> electrical_potential_host(grid);
        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            time_stepper
                    .update(Kokkos::DefaultHostExecutionSpace(),
                            density_host,
                            dt,
                            define_advection_field,
                            advect_density);


            m_builder(get_field(density_coef), get_const_field(density_host));
            m_poisson_solver(charge_density_coord, get_field(electrical_potential));
            ddc::parallel_deepcopy(electrical_potential_host, electrical_potential);
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


