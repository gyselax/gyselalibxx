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
 * @brief A second order implicit predictor-corrector for the Vlasov-Poisson equations.
 *
 * It solves in time the following Vlasov-Poisson equations system:
 *
 * - @f$  - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho @f$,
 * - @f$ E = - \nabla \phi @f$,
 * - @f$ \partial_t \rho - E_y \partial_x \rho + E_x \partial_y \rho = 0@f$,
 *
 * we write @f$ (A_x, A_y) =  (-E_y, E_x)@f$.
 *
 * The second order implicit predictor-corrector is also detailed in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * for @f$ n \geq 0 @f$,
 *
 * First, it predicts:
 * - 1. From @f$\rho^n@f$, it computes @f$\phi^n@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 2. From @f$\phi^n@f$, it computes @f$A^n@f$ with a AdvectionFieldFinder;
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes implicitly @f$\rho^P@f$ with a BslAdvectionPolar on @f$ \frac{dt}{4} @f$:
 *      - the characteristic feet @f$X^P@f$ is such that @f$X^P = X^k@f$ with @f$X^k@f$ the result of the implicit method:
 *          - @f$ X^k = X^n - \frac{dt}{4} \partial_t X^k@f$.
 * 
 * Secondly, it corrects: 
 * - 4. From @f$\rho^P@f$, it computes @f$\phi^P@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 5. From @f$\phi^P@f$, it computes @f$A^P@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$ A^{P} @f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionPolar on @f$ \frac{dt}{2} @f$.
 *      - the characteristic feet @f$X^C@f$ is such that @f$X^C = X^k@f$ with @f$X^k@f$ the result of the implicit method:
 *          - @f$\partial_t X^k = A^P(X^n) + A^P(X^{k-1}) @f$,
 *
 *
 * @tparam LogicalToPhysicalMapping
 *      A class describing a mapping from curvilinear coordinates to Cartesian coordinates.
 * @tparam LogicalToPseudoPhysicalMapping
 *      A class describing a mapping from curvilinear coordinates to pseudo-Cartesian coordinates.
 */
template <class LogicalToPhysicalMapping, class LogicalToPseudoPhysicalMapping>
class BslImplicitPredCorrRTheta : public ITimeSolverRTheta
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
    SplinePolarFootFinderType_host const m_foot_finder;

    PolarSplineFEMPoissonLikeSolver<
            GridR,
            GridTheta,
            PolarBSplinesRTheta,
            SplineRThetaEvaluatorNullBound> const& m_poisson_solver;

    SplineRThetaBuilder_host const& m_builder;
    SplineRThetaEvaluatorConstBound_host const& m_evaluator;



public:
    /**
     * @brief Instantiate a BslImplicitPredCorrRTheta.
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
     *      advection field and the rhs.
     * @param[in] poisson_solver
     *      The PDE solver which computes the electrical
     *      potential.
     * @param[in] advection_evaluator
     *      An evaluator of B-splines for the spline advection field.
     */
    BslImplicitPredCorrRTheta(
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
        , m_foot_finder(
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
        IdxRangeRTheta const grid(get_idx_range(density_host));

        // --- Electrostatic potential (phi). -------------------------------------------------------------
        DFieldMemRTheta electrical_potential_alloc(grid);
        host_t<DFieldMemRTheta> electrical_potential_alloc_host(grid);

        host_t<PolarSplineMemRTheta> electrostatic_potential_coef_alloc_host(
                ddc::discrete_space<PolarBSplinesRTheta>().full_domain());

        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_alloc_host(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_k_alloc_host(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_k_tot_alloc_host(grid);
        host_t<VectorSplineCoeffsMem2D<X, Y>> advection_field_coefs_k_alloc_host(
                get_spline_idx_range(m_builder));
        host_t<FieldMemRTheta<CoordRTheta>> feet_coords_alloc_host(grid);
        DFieldMemRTheta density_predicted_alloc(grid);
        host_t<Spline2DMem> density_coef_alloc_host(get_spline_idx_range(m_builder));

        auto advection_field_k_tot_alloc = ddcHelper::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(advection_field_k_tot_alloc_host));
        auto density_predicted_alloc_host
                = ddc::create_mirror_view(get_field(density_predicted_alloc));
        auto density_alloc = ddc::create_mirror_view(Kokkos::DefaultExecutionSpace(), density_host);

        host_t<DVectorFieldRTheta<X, Y>> advection_field_host(advection_field_alloc_host);
        host_t<DVectorFieldRTheta<X, Y>> advection_field_k_host(advection_field_k_alloc_host);
        host_t<DVectorFieldRTheta<X, Y>> advection_field_k_tot_host(
                advection_field_k_tot_alloc_host);
        DFieldRTheta density_predicted = get_field(density_predicted_alloc);
        host_t<DFieldRTheta> density_predicted_host = get_field(density_predicted_alloc_host);
        DFieldRTheta density = get_field(density_alloc);
        host_t<Spline2D> density_coef_host(density_coef_alloc_host);

        // Operators
        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<
                Kokkos::DefaultHostExecutionSpace,
                Kokkos::HostSpace,
                PolarBSplinesRTheta,
                ddc::NullExtrapolationRule>
                polar_spline_evaluator(extrapolation_rule);

        AdvectionFieldFinder advection_field_computer(m_logical_to_physical);

        PoissonLikeRHSFunction const
                charge_density(get_const_field(density_coef_host), m_evaluator);

        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            m_builder(density_coef_host, get_const_field(density_host));
            m_poisson_solver(
                    charge_density,
                    get_field(electrostatic_potential_coef_alloc_host));

            polar_spline_evaluator(
                    get_field(electrical_potential_alloc_host),
                    get_const_field(get_field(electrostatic_potential_coef_alloc_host)));

            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .with("time", iter * dt)
                    .with("density", density_host)
                    .with("electrical_potential", electrical_potential_alloc_host);


            // STEP 2: From phi^n, we compute A^n:
            advection_field_computer(
                    get_field(electrostatic_potential_coef_alloc_host),
                    advection_field_host);


            // STEP 3: From rho^n and A^n, we compute rho^P: Vlasov equation
            m_builder(
                    ddcHelper::get<X>(advection_field_coefs_k_alloc_host),
                    ddcHelper::get<X>(get_const_field(advection_field_host)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs_k_alloc_host),
                    ddcHelper::get<Y>(get_const_field(advection_field_host)));

            // initialisation:
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                feet_coords_alloc_host(irtheta) = CoordRTheta(ddc::coordinate(irtheta));
            });

            const double tau = 1e-6;
            implicit_loop(
                    advection_field_host,
                    get_const_field(advection_field_coefs_k_alloc_host),
                    get_field(feet_coords_alloc_host),
                    dt / 4.,
                    tau);

            // Evaluate A^n at X^P:
            m_evaluator(
                    ddcHelper::get<X>(advection_field_k_host),
                    get_const_field(feet_coords_alloc_host),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs_k_alloc_host)));
            m_evaluator(
                    ddcHelper::get<Y>(advection_field_k_host),
                    get_const_field(feet_coords_alloc_host),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs_k_alloc_host)));

            // Compute the new advection field (E^n(X^n) + E^n(X^P)) /2:
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                ddcHelper::assign_vector_field_element(
                        advection_field_k_tot_host,
                        irtheta,
                        (advection_field_host(irtheta) + advection_field_k_host(irtheta)) / 2.0);
            });


            // X^P = X^n - dt/2 * ( E^n(X^n) + E^n(X^P) )/2:
            // --- Copy rho^n because it will be modified:
            ddc::parallel_deepcopy(density_predicted, density_host);
            ddcHelper::deepcopy(
                    get_field(advection_field_k_tot_alloc),
                    get_const_field(advection_field_k_tot_host));
            m_advection_solver(
                    density_predicted,
                    get_const_field(advection_field_k_tot_alloc),
                    dt / 2.);

            // --- advect also the feet because it is needed for the next step
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                feet_coords_alloc_host(irtheta) = CoordRTheta(ddc::coordinate(irtheta));
            });
            m_foot_finder(
                    get_field(feet_coords_alloc_host),
                    get_const_field(advection_field_k_tot_host),
                    dt / 2.);


            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            ddc::parallel_deepcopy(density_predicted_host, density_predicted);
            m_builder(density_coef_host, get_const_field(density_predicted_host));
            m_poisson_solver(
                    charge_density,
                    get_field(electrostatic_potential_coef_alloc_host));

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(
                    get_field(electrostatic_potential_coef_alloc_host),
                    advection_field_host);


            // STEP 6: From rho^n and A^P, we compute rho^{n+1}: Vlasov equation
            m_builder(
                    ddcHelper::get<X>(advection_field_coefs_k_alloc_host),
                    ddcHelper::get<X>(get_const_field(advection_field_host)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs_k_alloc_host),
                    ddcHelper::get<Y>(get_const_field(advection_field_host)));


            // initialisation:
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                feet_coords_alloc_host(irtheta) = CoordRTheta(ddc::coordinate(irtheta));
            });

            implicit_loop(
                    advection_field_host,
                    get_const_field(advection_field_coefs_k_alloc_host),
                    get_field(feet_coords_alloc_host),
                    dt / 2.,
                    tau);

            // Evaluate A^P at X^P:
            m_evaluator(
                    ddcHelper::get<X>(advection_field_k_host),
                    get_const_field(feet_coords_alloc_host),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs_k_alloc_host)));
            m_evaluator(
                    ddcHelper::get<Y>(advection_field_k_host),
                    get_const_field(feet_coords_alloc_host),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs_k_alloc_host)));

            // Computed advection field (A^P(X^n) + A^P(X^P)) /2:
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                ddcHelper::assign_vector_field_element(
                        advection_field_k_tot_host,
                        irtheta,
                        (advection_field_host(irtheta) + advection_field_k_host(irtheta)) / 2.);
            });
            // X^k = X^n - dt * ( A^P(X^n) + A^P(X^P) )/2
            ddc::parallel_deepcopy(density, get_const_field(density_host));
            ddcHelper::deepcopy(advection_field_k_tot_alloc, advection_field_k_tot_host);
            m_advection_solver(density, get_const_field(advection_field_k_tot_alloc), dt);
            ddc::parallel_deepcopy(density_host, get_const_field(density));
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        m_builder(density_coef_host, get_const_field(density_host));
        m_poisson_solver(charge_density, get_field(electrical_potential_alloc));
        ddc::parallel_deepcopy(
                get_field(electrical_potential_alloc_host),
                get_const_field(electrical_potential_alloc));

        ddc::PdiEvent("last_iteration")
                .with("iter", steps)
                .with("time", steps * dt)
                .with("density", density_host)
                .with("electrical_potential", electrical_potential_alloc_host);

        end_time = std::chrono::system_clock::now();
        display_time_difference("Iterations time: ", start_time, end_time);


        return density_host;
    }



private:
    double compute_square_polar_distance(CoordRTheta const& coord1, CoordRTheta const& coord2) const
    {
        CoordXY coord_xy1(m_logical_to_physical(coord1));
        CoordXY coord_xy2(m_logical_to_physical(coord2));

        const double x1 = ddc::select<X>(coord_xy1);
        const double y1 = ddc::select<Y>(coord_xy1);
        const double x2 = ddc::select<X>(coord_xy2);
        const double y2 = ddc::select<Y>(coord_xy2);

        return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    }


    void implicit_loop(
            host_t<DVectorFieldRTheta<X, Y>> advection_field_host,
            host_t<ConstVectorSplineCoeffs2D<X, Y>> advection_field_coefs_k_host,
            host_t<FieldRTheta<CoordRTheta>> feet_coords_host,
            double const dt,
            double const tau) const
    {
        IdxRangeRTheta const grid = get_idx_range(advection_field_host);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_k_alloc_host(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_k_tot_alloc_host(grid);
        host_t<FieldMemRTheta<CoordRTheta>> feet_coords_tmp_alloc_host(grid);

        double square_difference_feet = 0.;
        int count = 0;
        const int max_count = 50;
        do {
            count++;

            // Evaluate A at X^{k-1}:
            m_evaluator(
                    ddcHelper::get<X>(advection_field_k_alloc_host),
                    get_const_field(feet_coords_host),
                    ddcHelper::get<X>(advection_field_coefs_k_host));
            m_evaluator(
                    ddcHelper::get<Y>(advection_field_k_alloc_host),
                    get_const_field(feet_coords_host),
                    ddcHelper::get<Y>(advection_field_coefs_k_host));

            // Compute the new advection field A(X^n) + A(X^{k-1}):
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                ddcHelper::assign_vector_field_element(
                        get_field(advection_field_k_tot_alloc_host),
                        irtheta,
                        advection_field_host(irtheta) + advection_field_k_alloc_host(irtheta));
            });

            // X^{k-1} = X^k:
            ddc::parallel_deepcopy(feet_coords_tmp_alloc_host, feet_coords_host);

            // X^k = X^n - dt* X^k:
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                feet_coords_host(irtheta) = CoordRTheta(ddc::coordinate(irtheta));
            });
            m_foot_finder(feet_coords_host, get_const_field(advection_field_k_tot_alloc_host), dt);


            // Convergence test:
            square_difference_feet = 0.;
            ddc::for_each(grid, [&](IdxRTheta const irtheta) {
                double sqr_diff_feet = compute_square_polar_distance(
                        feet_coords_host(irtheta),
                        feet_coords_tmp_alloc_host(irtheta));
                square_difference_feet = square_difference_feet > sqr_diff_feet
                                                 ? square_difference_feet
                                                 : sqr_diff_feet;
            });

        } while ((square_difference_feet > tau * tau) and (count < max_count));
    }
};
