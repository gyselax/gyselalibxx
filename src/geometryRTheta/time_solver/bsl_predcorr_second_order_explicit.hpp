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



/**
 * @brief A second order explicit predictor-corrector for the Vlasov-Poisson equations.
 *
 * It solves in time the following Vlasov-Poisson equations system:
 *
 * - @f$  - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho @f$,
 * - @f$ E = - \nabla \phi @f$,
 * - @f$ \partial_t \rho - E_y \partial_x \rho + E_x \partial_y \rho = 0@f$,
 *
 * we write @f$ (A_x, A_y) =  (-E_y, E_x)@f$.
 *
 * The second order explicit predictor-corrector is also detailed in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * for @f$ n \geq 0 @f$,
 *
 * First, it predicts:
 * - 1. From @f$\rho^n@f$, it computes @f$\phi^n@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 2. From @f$\phi^n@f$, it computes @f$A^n@f$ with a AdvectionFieldFinder;
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes @f$\rho^P@f$ with a BslAdvectionPolar on @f$ dt @f$;
 *
 * We write @f$X^P@f$ the characteristic feet such that @f$\partial_t X^P = A^n(X^n)@f$.
 *
 * Secondly, it corrects:
 * - 4. From @f$\rho^P@f$, it computes @f$\phi^P@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 5. From @f$\phi^P@f$, it computes @f$A^P@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$\frac{A^{P}(X^n) + A^n(X^P)}{2} @f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionPolar on @f$ dt @f$.
 *
 * (With @f$X^C@f$ the characteristic feet such that @f$\partial_t X^C = \frac{A^{P}(X^n) + A^n(X^P)}{2} @f$.)
 *
 * @tparam LogicalToPhysicalMapping
 *      A class describing a mapping from curvilinear coordinates to Cartesian coordinates.
 * @tparam LogicalToPseudoPhysicalMapping
 *      A class describing a mapping from curvilinear coordinates to pseudo-Cartesian coordinates.
 */
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
    SplinePolarFootFinderType const m_find_feet;

    PolarSplineFEMPoissonLikeSolver<
            GridR,
            GridTheta,
            PolarBSplinesRTheta,
            SplineRThetaEvaluatorNullBound> const& m_poisson_solver;

    SplineRThetaBuilder const& m_builder;
    SplineRThetaEvaluatorConstBound const& m_evaluator;



public:
    /**
     * @brief Instantiate a BslExplicitPredCorrRTheta.
     *
     * @param[in] logical_to_physical
     *      The mapping from the logical domain to the physical domain.
     * @param[in] logical_to_pseudo_physical
     *      The mapping from the logical domain to the pseudo-physical domain.
     * @param[in] advection_solver
     *      The advection operator with an Euler method.
     * @param[in] grid
     *      The index range on which the functions are defined.
     * @param[in] builder
     *      A spline builder to get the spline representation of the
     *      advection field and the RHS.
     * @param[in] poisson_solver
     *      The PDE solver which computes the electrical
     *      potential.
     * @param[in] advection_evaluator
     *      An evaluator of B-splines for the spline advection field.
     */
    BslExplicitPredCorrRTheta(
            LogicalToPhysicalMapping const& logical_to_physical,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical,
            BslAdvectionRTheta const& advection_solver,
            IdxRangeRTheta const& grid,
            SplineRThetaBuilder const& builder,
            PolarSplineFEMPoissonLikeSolver<
                    GridR,
                    GridTheta,
                    PolarBSplinesRTheta,
                    SplineRThetaEvaluatorNullBound> const& poisson_solver,
            SplineRThetaEvaluatorConstBound const& advection_evaluator)
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

        PolarSplineMemRTheta electrostatic_potential_coef_alloc(
                ddc::discrete_space<PolarBSplinesRTheta>().full_domain());

        auto electrostatic_potential_coef_alloc_host
                = ddc::create_mirror_view(get_field(electrostatic_potential_coef_alloc));

        Spline2DMem density_coef_alloc(get_spline_idx_range(m_builder));
        DFieldMemRTheta density_predicted_alloc(grid);
        auto density_predicted_host_alloc
                = ddc::create_mirror_view(get_field(density_predicted_alloc));
        auto density_alloc = ddc::create_mirror_view(Kokkos::DefaultExecutionSpace(), density_host);
        FieldMemRTheta<CoordRTheta> feet_coords_alloc(grid);
        DVectorFieldMemRTheta<X, Y> advection_field_evaluated_alloc(grid);
        VectorSplineCoeffsMem2D<X, Y> advection_field_coefs_alloc(get_spline_idx_range(m_builder));


        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_alloc_host(grid);
        DVectorFieldMemRTheta<X, Y> advection_field_predicted_alloc(grid);
        auto advection_field_alloc = ddcHelper::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(advection_field_alloc_host));

        // Create fields
        host_t<DVectorFieldRTheta<X, Y>> advection_field_host(advection_field_alloc_host);
        DVectorFieldRTheta<X, Y> advection_field_predicted(advection_field_predicted_alloc);
        DVectorFieldRTheta<X, Y> advection_field(advection_field_alloc);
        DVectorFieldRTheta<X, Y> advection_field_evaluated(advection_field_evaluated_alloc);

        FieldRTheta<CoordRTheta> feet_coords(feet_coords_alloc);

        Spline2D density_coef(density_coef_alloc);
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

        PoissonLikeRHSFunction const charge_density(get_const_field(density_coef), m_evaluator);

        ddc::parallel_deepcopy(density, get_const_field(density_host));


        // --- Parameter for linearisation of advection field: --------------------------------------------
        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            double const time = iter * dt;

            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            m_builder(density_coef, get_const_field(density));
            m_poisson_solver(charge_density, get_field(electrostatic_potential_coef_alloc));

            polar_spline_evaluator(
                    get_field(electrical_potential),
                    get_const_field(electrostatic_potential_coef_alloc));

            ddc::parallel_deepcopy(
                    get_field(electrical_potential_host),
                    get_const_field(electrical_potential));

            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .with("time", time)
                    .with("density", density_host)
                    .with("electrical_potential", electrical_potential_host);

            ddc::parallel_deepcopy(
                    get_field(electrostatic_potential_coef_alloc_host),
                    get_const_field(electrostatic_potential_coef_alloc));

            // STEP 2: From phi^n, we compute A^n:
            advection_field_computer(
                    get_field(electrostatic_potential_coef_alloc_host),
                    advection_field_host);

            ddcHelper::deepcopy(advection_field, advection_field_host);

            // STEP 3: From rho^n and A^n, we compute rho^P: Vlasov equation
            // --- Copy rho^n because it will be modified:
            ddc::parallel_deepcopy(get_field(density_predicted_alloc), density);
            m_advection_solver(get_field(density_predicted_alloc), advection_field, dt);

            // --- advect also the feet because it is needed for the next step
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        feet_coords(irtheta) = ddc::coordinate(irtheta);
                    });
            m_find_feet(feet_coords, advection_field, dt);

            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            m_builder(density_coef, get_const_field(density_predicted_alloc));
            m_poisson_solver(charge_density, get_field(electrostatic_potential_coef_alloc));

            ddc::parallel_deepcopy(
                    get_field(electrostatic_potential_coef_alloc_host),
                    get_const_field(electrostatic_potential_coef_alloc));

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(
                    get_field(electrostatic_potential_coef_alloc_host),
                    advection_field_host);

            ddcHelper::deepcopy(advection_field_predicted, advection_field_host);


            // ---  we evaluate the advection field A^n at the characteristic feet X^P
            m_builder(
                    ddcHelper::get<X>(advection_field_coefs_alloc),
                    ddcHelper::get<X>(get_const_field(advection_field_predicted)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs_alloc),
                    ddcHelper::get<Y>(get_const_field(advection_field_predicted)));

            m_evaluator(
                    ddcHelper::get<X>(advection_field_evaluated),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs_alloc)));
            m_evaluator(
                    ddcHelper::get<Y>(advection_field_evaluated),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs_alloc)));


            // STEP 6: From rho^n and (A^n(X^P) + A^P(X^n))/2, we compute rho^{n+1}: Vlasov equation
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        ddcHelper::assign_vector_field_element(
                                advection_field,
                                irtheta,
                                (advection_field_evaluated(irtheta)
                                 + advection_field_predicted(irtheta))
                                        / 2.);
                    });

            m_advection_solver(get_field(density), get_const_field(advection_field), dt);
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        m_builder(density_coef, get_const_field(density));
        m_poisson_solver(charge_density, get_field(electrical_potential));

        ddc::parallel_deepcopy(electrical_potential_host, electrical_potential);
        ddc::parallel_deepcopy(density_host, density);

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
