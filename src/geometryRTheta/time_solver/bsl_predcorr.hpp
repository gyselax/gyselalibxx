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

    SplineRThetaBuilder_host const& m_builder;
    SplineRThetaEvaluatorNullBound_host const& m_spline_evaluator;


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
