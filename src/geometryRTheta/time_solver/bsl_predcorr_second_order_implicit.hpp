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
    SplinePolarFootFinderType const m_foot_finder;

    PolarSplineFEMPoissonLikeSolver<
            GridR,
            GridTheta,
            PolarBSplinesRTheta,
            SplineRThetaEvaluatorNullBound> const& m_poisson_solver;

    SplineRThetaBuilder const& m_builder;
    SplineRThetaEvaluatorConstBound const& m_evaluator;



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
            SplineRThetaBuilder const& builder,
            PolarSplineFEMPoissonLikeSolver<
                    GridR,
                    GridTheta,
                    PolarBSplinesRTheta,
                    SplineRThetaEvaluatorNullBound> const& poisson_solver,
            SplineRThetaEvaluatorConstBound const& advection_evaluator)
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

        PolarSplineMemRTheta electrostatic_potential_coef_alloc(
                ddc::discrete_space<PolarBSplinesRTheta>().full_domain());

        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        DVectorFieldMemRTheta<X, Y> advection_field_alloc(grid);
        DVectorFieldMemRTheta<X, Y> advection_field_k_alloc(grid);
        DVectorFieldMemRTheta<X, Y> advection_field_k_tot_alloc(grid);
        VectorSplineCoeffsMem2D<X, Y> advection_field_coefs_k_alloc(
                get_spline_idx_range(m_builder));
        FieldMemRTheta<CoordRTheta> feet_coords_alloc(grid);
        DFieldMemRTheta density_predicted_alloc(grid);
        Spline2DMem density_coef_alloc(get_spline_idx_range(m_builder));

        auto advection_field_alloc_host = ddcHelper::create_mirror_view_and_copy(
                Kokkos::DefaultHostExecutionSpace(),
                get_field(advection_field_alloc));
        auto density_alloc
                = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), density_host);
        auto electrostatic_potential_coef_alloc_host
                = ddc::create_mirror_view(get_field(electrostatic_potential_coef_alloc));

        DVectorFieldRTheta<X, Y> advection_field(advection_field_alloc);
        DVectorFieldRTheta<X, Y> advection_field_k(advection_field_k_alloc);
        DVectorFieldRTheta<X, Y> advection_field_k_tot(advection_field_k_tot_alloc);
        host_t<DVectorFieldRTheta<X, Y>> advection_field_host(advection_field_alloc_host);
        VectorSplineCoeffs2D<X, Y> advection_field_coefs_k(advection_field_coefs_k_alloc);
        DFieldRTheta density_predicted = get_field(density_predicted_alloc);
        DFieldRTheta density = get_field(density_alloc);
        Spline2D density_coef(density_coef_alloc);
        FieldRTheta<CoordRTheta> feet_coords(feet_coords_alloc);

        PolarSplineRTheta electrostatic_potential_coef(electrostatic_potential_coef_alloc);
        host_t<PolarSplineRTheta> electrostatic_potential_coef_host(
                electrostatic_potential_coef_alloc_host);

        // Operators
        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<
                Kokkos::DefaultExecutionSpace,
                Kokkos::DefaultExecutionSpace::memory_space,
                PolarBSplinesRTheta,
                ddc::NullExtrapolationRule>
                polar_spline_evaluator(extrapolation_rule);

        AdvectionFieldFinder advection_field_computer(m_logical_to_physical);

        PoissonLikeRHSFunction const charge_density(get_const_field(density_coef), m_evaluator);

        ddc::parallel_deepcopy(density, get_const_field(density_host));

        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            m_builder(density_coef, get_const_field(density));
            m_poisson_solver(charge_density, electrostatic_potential_coef);

            polar_spline_evaluator(
                    get_field(electrical_potential_alloc),
                    get_const_field(electrostatic_potential_coef));

            ddc::parallel_deepcopy(
                    get_field(electrical_potential_alloc_host),
                    get_const_field(electrical_potential_alloc));
            ddc::parallel_deepcopy(density_host, get_const_field(density));

            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .with("time", iter * dt)
                    .with("density", density_host)
                    .with("electrical_potential", electrical_potential_alloc_host);

            ddc::parallel_deepcopy(
                    electrostatic_potential_coef_host,
                    get_const_field(electrostatic_potential_coef));

            // STEP 2: From phi^n, we compute A^n:
            advection_field_computer(electrostatic_potential_coef_host, advection_field_host);

            ddcHelper::deepcopy(advection_field, get_const_field(advection_field_host));

            // STEP 3: From rho^n and A^n, we compute rho^P: Vlasov equation
            m_builder(
                    ddcHelper::get<X>(advection_field_coefs_k),
                    ddcHelper::get<X>(get_const_field(advection_field)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs_k),
                    ddcHelper::get<Y>(get_const_field(advection_field)));

            // initialisation:
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        feet_coords(irtheta) = ddc::coordinate(irtheta);
                    });

            const double tau = 1e-6;
            implicit_loop(
                    advection_field,
                    get_const_field(advection_field_coefs_k),
                    feet_coords,
                    dt / 4.,
                    tau);

            // Evaluate A^n at X^P:
            m_evaluator(
                    ddcHelper::get<X>(advection_field_k),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs_k)));
            m_evaluator(
                    ddcHelper::get<Y>(advection_field_k),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs_k)));

            // Compute the new advection field (E^n(X^n) + E^n(X^P)) /2:
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        ddcHelper::assign_vector_field_element(
                                advection_field_k_tot,
                                irtheta,
                                (advection_field(irtheta) + advection_field_k(irtheta)) / 2.0);
                    });


            // X^P = X^n - dt/2 * ( E^n(X^n) + E^n(X^P) )/2:
            // --- Copy rho^n because it will be modified:
            ddc::parallel_deepcopy(density_predicted, density);
            m_advection_solver(density_predicted, get_const_field(advection_field_k), dt / 2.);

            // --- advect also the feet because it is needed for the next step
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        feet_coords(irtheta) = ddc::coordinate(irtheta);
                    });
            m_foot_finder(feet_coords, get_const_field(advection_field_k_tot), dt / 2.);


            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            m_builder(density_coef, get_const_field(density_predicted));
            m_poisson_solver(charge_density, electrostatic_potential_coef);

            ddc::parallel_deepcopy(
                    electrostatic_potential_coef_host,
                    get_const_field(electrostatic_potential_coef));

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(electrostatic_potential_coef_host, advection_field_host);

            ddcHelper::deepcopy(advection_field, get_const_field(advection_field_host));


            // STEP 6: From rho^n and A^P, we compute rho^{n+1}: Vlasov equation
            m_builder(
                    ddcHelper::get<X>(advection_field_coefs_k),
                    ddcHelper::get<X>(get_const_field(advection_field)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs_k),
                    ddcHelper::get<Y>(get_const_field(advection_field)));


            // initialisation:
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        feet_coords(irtheta) = ddc::coordinate(irtheta);
                    });

            implicit_loop(
                    advection_field,
                    get_const_field(advection_field_coefs_k),
                    get_field(feet_coords_alloc),
                    dt / 2.,
                    tau);

            // Evaluate A^P at X^P:
            m_evaluator(
                    ddcHelper::get<X>(advection_field_k),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs_k)));
            m_evaluator(
                    ddcHelper::get<Y>(advection_field_k),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs_k)));

            // Computed advection field (A^P(X^n) + A^P(X^P)) /2:
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        ddcHelper::assign_vector_field_element(
                                advection_field_k_tot,
                                irtheta,
                                (advection_field(irtheta) + advection_field_k(irtheta)) / 2.);
                    });
            // X^k = X^n - dt * ( A^P(X^n) + A^P(X^P) )/2
            m_advection_solver(density, get_const_field(advection_field_k_tot), dt);
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        m_builder(density_coef, get_const_field(density));
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

        ddc::parallel_deepcopy(density_host, get_const_field(density));
        return density_host;
    }


    /**
     * @brief The implicit loop which calculates the feet of the charateristics.
     *
     * This function should be private but cannot be as it contains Kokkos lambda
     * functions.
     *
     * @param[in] advection_field The current advection field.
     * @param[in] advection_field_coefs_k The coefficients of the spline representation of the current advection field.
     * @param[inout] feet_coords The feet of the characteristics.
     * @param[in] dt The time step.
     * @param[in] tau The stopping criteria. The convergence is reached when the old and new
     *                feet of the characteristics are less than tau apart.
     */
    void implicit_loop(
            DVectorConstFieldRTheta<X, Y> advection_field,
            ConstVectorSplineCoeffs2D<X, Y> advection_field_coefs_k,
            FieldRTheta<CoordRTheta> feet_coords,
            double const dt,
            double const tau) const
    {
        IdxRangeRTheta const grid = get_idx_range(advection_field);
        DVectorFieldMemRTheta<X, Y> advection_field_k_alloc(grid);
        DVectorFieldMemRTheta<X, Y> advection_field_k_tot_alloc(grid);
        FieldMemRTheta<CoordRTheta> feet_coords_tmp_alloc(grid);

        DVectorFieldRTheta<X, Y> advection_field_k(advection_field_k_alloc);
        DVectorFieldRTheta<X, Y> advection_field_k_tot(advection_field_k_tot_alloc);
        FieldRTheta<CoordRTheta> feet_coords_tmp(feet_coords_tmp_alloc);

        double square_difference_feet = 0.;
        int count = 0;
        const int max_count = 50;
        do {
            count++;

            // Evaluate A at X^{k-1}:
            m_evaluator(
                    ddcHelper::get<X>(advection_field_k),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(advection_field_coefs_k));
            m_evaluator(
                    ddcHelper::get<Y>(advection_field_k),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(advection_field_coefs_k));

            // Compute the new advection field A(X^n) + A(X^{k-1}):
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        ddcHelper::assign_vector_field_element(
                                advection_field_k_tot,
                                irtheta,
                                advection_field(irtheta) + advection_field_k(irtheta));
                    });

            // X^{k-1} = X^k:
            ddc::parallel_deepcopy(feet_coords_tmp, feet_coords);

            // X^k = X^n - dt* X^k:
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        feet_coords(irtheta) = ddc::coordinate(irtheta);
                    });
            m_foot_finder(feet_coords, get_const_field(advection_field_k_tot), dt);


            // Convergence test:
            LogicalToPhysicalMapping logical_to_physical_proxy = m_logical_to_physical;
            square_difference_feet = ddc::parallel_transform_reduce(
                    Kokkos::DefaultExecutionSpace(),
                    grid,
                    0.0,
                    ddc::reducer::max<double>(),
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        DVector<X, Y> distance(
                                logical_to_physical_proxy(feet_coords(irtheta))
                                - logical_to_physical_proxy(feet_coords_tmp(irtheta)));
                        return scalar_product(distance, distance);
                    });

        } while ((square_difference_feet > tau * tau) and (count < max_count));
    }
};
