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

/**
 * @brief Predictor-corrector for the Vlasov-Poisson equations.
 *
 * It solves in time the following Vlasov-Poisson equations system:
 *
 * - @f$  - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho @f$,
 * - @f$ E = - \nabla \phi @f$,
 * - @f$ \partial_t \rho - E_y \partial_x \rho + E_x \partial_y \rho = 0@f$,
 *
 * we write @f$ (A_x, A_y) =  (-E_y, E_x)@f$.
 *
 * This method is mainly a Runge-Kutta 2 method:
 *
 * for @f$ n \geq 0 @f$,
 *
 * First, it advects on a half time step:
 * - 1. From @f$\rho^n@f$, it computes @f$\phi^n@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 2. From @f$\phi^n@f$, it computes @f$A^n@f$ with a AdvectionFieldFinder;
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes @f$\rho^{n+1/2}@f$ with a BslAdvectionPolar on @f$\frac{dt}{2}@f$;
 *
 * Secondly, it advects on a full time step:
 * - 4. From @f$\rho^{n+1/2}@f$, it computes @f$\phi^{n+1/2}@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 5. From @f$\phi^{n+1/2}@f$, it computes @f$A^{n+1/2}@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$A^{n+1/2}@f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionPolar on @f$dt@f$.
 *
 * @tparam Mapping
 *      A class describing a mapping from curvilinear coordinates to Cartesian coordinates.
 * @tparam FootFinder
 *      A IFootFinder class.
 *
 */
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

    SplineRThetaBuilder const& m_builder;
    SplineRThetaEvaluatorNullBound const& m_spline_evaluator;


public:
    /**
     * @brief Instantiate a BslPredCorrRTheta.
     *
     * @param[in] mapping
     *      The mapping function from the logical domain to the
     *      physical domain.
     * @param[in] advection_solver
     *      The advection operator.
     * @param[in] builder
     *      The spline builder for the computation of the RHS
     *      and the advection field.
     * @param[in] rhs_evaluator
     *      The evaluator of B-splines for the RHS.
     * @param[in] poisson_solver
     *      The PDE solver which computes the electrical
     *      potential.
     */
    BslPredCorrRTheta(
            Mapping const& mapping,
            BslAdvectionRTheta const& advection_solver,
            SplineRThetaBuilder const& builder,
            SplineRThetaEvaluatorNullBound const& rhs_evaluator,
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
        Spline2DMem density_coef_alloc(get_spline_idx_range(m_builder));
        PolarSplineMemRTheta electrostatic_potential_coef_alloc(
                ddc::discrete_space<PolarBSplinesRTheta>().full_domain());

        auto electrical_potential_alloc_host
                = ddc::create_mirror_view(get_field(electrical_potential_alloc));
        auto density_alloc = ddc::create_mirror_view(Kokkos::DefaultExecutionSpace(), density_host);

        // Fields
        DFieldRTheta electrical_potential = get_field(electrical_potential_alloc);
        host_t<DFieldRTheta> electrical_potential_host = get_field(electrical_potential_alloc_host);
        PolarSplineRTheta electrostatic_potential_coef
                = get_field(electrostatic_potential_coef_alloc);
        DFieldRTheta density = get_field(density_alloc);
        Spline2D density_coef(density_coef_alloc);

        // Operators
        AdvectionFieldFinder advection_field_computer(m_mapping);
        PoissonLikeRHSFunction const
                charge_density(get_const_field(density_coef), m_spline_evaluator);

        ddc::parallel_deepcopy(density, get_const_field(density_host));

        // Setup
        m_builder(density_coef, get_const_field(density));
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
                    // --- compute electrostatic potential:
                    m_builder(density_coef, get_const_field(density));
                    m_poisson_solver(charge_density, electrostatic_potential_coef);

                    auto advection_field_alloc_host = ddcHelper::create_mirror_view_and_copy(
                            Kokkos::DefaultHostExecutionSpace(),
                            advection_field);

                    auto electrostatic_potential_coef_alloc_host
                            = ddc::create_mirror_view_and_copy(electrostatic_potential_coef);

                    // --- compute advection field:
                    advection_field_computer(
                            get_field(electrostatic_potential_coef_alloc_host),
                            get_field(advection_field_alloc_host));

                    ddcHelper::
                            deepcopy(advection_field, get_const_field(advection_field_alloc_host));
                };

        RK2<DFieldMemRTheta, DVectorFieldMemRTheta<X, Y>, Kokkos::DefaultExecutionSpace>
                time_stepper(grid);

        // Iteration loop
        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            time_stepper
                    .update(Kokkos::DefaultExecutionSpace(),
                            density,
                            dt,
                            define_advection_field,
                            m_advection_solver);

            m_builder(density_coef, get_const_field(density));
            m_poisson_solver(charge_density, electrical_potential);
            ddc::parallel_deepcopy(density_host, get_const_field(density));
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
