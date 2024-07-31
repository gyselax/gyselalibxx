// SPDX-License-Identifier: MIT

#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <geometry.hpp>

#include "advection_domain.hpp"
#include "advection_field_rp.hpp"
#include "bsl_advection_rp.hpp"
#include "geometry.hpp"
#include "itimesolver.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "spline_foot_finder.hpp"
#include "spline_interpolator_2d_rp.hpp"



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
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes implicitly @f$\rho^P@f$ with a BslAdvectionRP on @f$ \frac{dt}{4} @f$:
 *      - the characteristic feet @f$X^P@f$ is such that @f$X^P = X^k@f$ with @f$X^k@f$ the result of the implicit method:
 *          - @f$ X^k = X^n - \frac{dt}{4} \partial_t X^k@f$.
 *
 * Secondly, it corrects:
 * - 4. From @f$\rho^P@f$, it computes @f$\phi^P@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 5. From @f$\phi^P@f$, it computes @f$A^P@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$ A^{P} @f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionRP on @f$ \frac{dt}{2} @f$.
 *      - the characteristic feet @f$X^C@f$ is such that @f$X^C = X^k@f$ with @f$X^k@f$ the result of the implicit method:
 *          - @f$\partial_t X^k = A^P(X^n) + A^P(X^{k-1}) @f$,
 *          - @f$ X^k = X^n - \frac{dt}{2} \partial_t X^k @f$.
 *
 *
 * @tparam Mapping
 *      A Curvilinear2DToCartesian class or one of its child classes.
 * @tparam AdvectionDomain
 *      An AdvectionDomain class.
 * @tparam FootFinder
 *      A IFootFinder class.
 *
 */
template <class Mapping, class AdvectionDomain>
class BslImplicitPredCorrRP : public ITimeSolverRP
{
private:
    using EulerMethod = Euler<FieldRP<CoordRP>, VectorDFieldRP<RDimX, RDimY>>;

    Mapping const& m_mapping;

    BslAdvectionRP<SplineFootFinder<EulerMethod, AdvectionDomain>, Mapping> const&
            m_advection_solver;

    EulerMethod const m_euler;
    SplineFootFinder<EulerMethod, AdvectionDomain> const m_foot_finder;

    PolarSplineFEMPoissonLikeSolver const& m_poisson_solver;

    SplineRPBuilder const& m_builder;
    SplineRPEvaluatorConstBound const& m_evaluator;



public:
    /**
     * @brief Instantiate a BslImplicitPredCorrRP.
     *
     * @param[in] advection_domain
     *      An AdvectionDomain object which gives the information
     *      in which domain we advect.
     * @param[in] mapping
     *      The mapping function from the logical domain to the
     *      physical domain.
     * @param[in] advection_solver
     *      The advection operator with an Euler method.
     * @param[in] grid
     *      The domain on which the functions are defined.
     * @param[in] builder
     *      A spline builder to get the spline representation of the
     *      advection field and the rhs.
     * @param[in] rhs_evaluator
     *      The evaluator of B-splines for the rhs.
     * @param[in] poisson_solver
     *      The PDE solver which computes the electrical
     *      potential.
     * @param[in] advection_evaluator
     *      An evaluator of B-splines for the spline advection field.
     */
    BslImplicitPredCorrRP(
            AdvectionDomain const& advection_domain,
            Mapping const& mapping,
            BslAdvectionRP<SplineFootFinder<EulerMethod, AdvectionDomain>, Mapping> const&
                    advection_solver,
            IDomainRP const& grid,
            SplineRPBuilder const& builder,
            SplineRPEvaluatorNullBound const& rhs_evaluator,
            PolarSplineFEMPoissonLikeSolver const& poisson_solver,
            SplineRPEvaluatorConstBound const& advection_evaluator)
        : m_mapping(mapping)
        , m_advection_solver(advection_solver)
        , m_euler(grid)
        , m_foot_finder(m_euler, advection_domain, builder, advection_evaluator)
        , m_poisson_solver(poisson_solver)
        , m_builder(builder)
        , m_evaluator(advection_evaluator)

    {
    }


    ~BslImplicitPredCorrRP() {};



    DSpanRP operator()(DSpanRP allfdistribu, double const dt, int const steps) const
    {
        std::chrono::time_point<std::chrono::system_clock> start_time
                = std::chrono::system_clock::now();
        std::chrono::time_point<std::chrono::system_clock> end_time;

        // Grid. ------------------------------------------------------------------------------------------
        IDomainRP const grid(allfdistribu.domain<IDimR, IDimP>());

        FieldRP<CoordRP> coords(grid);
        ddc::for_each(grid, [&](IndexRP const irp) { coords(irp) = ddc::coordinate(irp); });

        BSDomainR radial_bsplines(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                ddc::DiscreteVector<BSplinesR> {PolarBSplinesRP::continuity + 1}));
        BSDomainP polar_domain(ddc::discrete_space<BSplinesP>().full_domain());

        // --- Electrostatic potential (phi). -------------------------------------------------------------
        DFieldRP electrical_potential(grid);

        SplinePolar electrostatic_potential_coef(
                PolarBSplinesRP::singular_idx_range<PolarBSplinesRP>(),
                BSDomainRP(radial_bsplines, polar_domain));

        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<PolarBSplinesRP, ddc::NullExtrapolationRule> polar_spline_evaluator(
                extrapolation_rule);

        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        VectorDFieldRP<RDimX, RDimY> electric_field(grid);
        VectorDFieldRP<RDimX, RDimY> advection_field(grid);

        AdvectionFieldFinder advection_field_computer(m_mapping);



        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            Spline2D allfdistribu_coef(m_builder.spline_domain());
            m_builder(allfdistribu_coef.span_view(), allfdistribu.span_cview());
            PoissonLikeRHSFunction const
                    charge_density_coord_1(allfdistribu_coef.span_cview(), m_evaluator);
            m_poisson_solver(charge_density_coord_1, electrostatic_potential_coef);

            polar_spline_evaluator(
                    electrical_potential.span_view(),
                    coords.span_cview(),
                    electrostatic_potential_coef);

            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .and_with("time", iter * dt)
                    .and_with("density", allfdistribu)
                    .and_with("electrical_potential", electrical_potential);


            // STEP 2: From phi^n, we compute A^n:
            advection_field_computer(electrostatic_potential_coef, advection_field);


            // STEP 3: From rho^n and A^n, we compute rho^P: Vlasov equation
            VectorDFieldRP<RDimX, RDimY> advection_field_k(grid);
            VectorDFieldRP<RDimX, RDimY> advection_field_k_tot(grid);

            VectorSpline2D<RDimX, RDimY> advection_field_coefs_k(m_builder.spline_domain());
            m_builder(
                    ddcHelper::get<RDimX>(advection_field_coefs_k),
                    ddcHelper::get<RDimX>(advection_field.span_cview()));
            m_builder(
                    ddcHelper::get<RDimY>(advection_field_coefs_k),
                    ddcHelper::get<RDimY>(advection_field.span_cview()));

            FieldRP<CoordRP> feet_coords(grid);
            FieldRP<CoordRP> feet_coords_tmp(grid);


            // initialisation:
            ddc::for_each(grid, [&](IndexRP const irp) {
                feet_coords(irp) = CoordRP(ddc::coordinate(irp));
            });

            const double tau = 1e-6;
            implicit_loop(advection_field, advection_field_coefs_k, feet_coords, dt / 4., tau);

            // Evaluate A^n at X^P:
            m_evaluator(
                    ddcHelper::get<RDimX>(advection_field_k).span_view(),
                    feet_coords.span_cview(),
                    ddcHelper::get<RDimX>(advection_field_coefs_k.span_cview()));
            m_evaluator(
                    ddcHelper::get<RDimY>(advection_field_k).span_view(),
                    feet_coords.span_cview(),
                    ddcHelper::get<RDimY>(advection_field_coefs_k.span_cview()));

            // Compute the new advection field (E^n(X^n) + E^n(X^P)) /2:
            ddc::for_each(grid, [&](IndexRP const irp) {
                ddcHelper::get<RDimX>(advection_field_k_tot)(irp)
                        = (ddcHelper::get<RDimX>(advection_field)(irp)
                           + ddcHelper::get<RDimX>(advection_field_k)(irp))
                          / 2.;
                ddcHelper::get<RDimY>(advection_field_k_tot)(irp)
                        = (ddcHelper::get<RDimY>(advection_field)(irp)
                           + ddcHelper::get<RDimY>(advection_field_k)(irp))
                          / 2.;
            });


            // X^P = X^n - dt/2 * ( E^n(X^n) + E^n(X^P) )/2:
            // --- Copy phi^n because it will be modified:
            DFieldRP allfdistribu_predicted(grid);
            ddc::parallel_deepcopy(allfdistribu_predicted, allfdistribu);
            m_advection_solver(
                    allfdistribu_predicted.span_view(),
                    advection_field_k_tot.span_cview(),
                    dt / 2.);

            // --- advect also the feet because it is needed for the next step
            ddc::for_each(grid, [&](IndexRP const irp) {
                feet_coords(irp) = CoordRP(ddc::coordinate(irp));
            });
            m_foot_finder(feet_coords.span_view(), advection_field_k_tot.span_cview(), dt / 2.);


            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            m_builder(allfdistribu_coef.span_view(), allfdistribu.span_cview());
            PoissonLikeRHSFunction const
                    charge_density_coord_4(allfdistribu_coef.span_cview(), m_evaluator);
            m_poisson_solver(charge_density_coord_4, electrostatic_potential_coef);

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(electrostatic_potential_coef, advection_field);


            // STEP 6: From rho^n and A^P, we compute rho^{n+1}: Vlasov equation
            m_builder(
                    ddcHelper::get<RDimX>(advection_field_coefs_k),
                    ddcHelper::get<RDimX>(advection_field.span_cview()));
            m_builder(
                    ddcHelper::get<RDimY>(advection_field_coefs_k),
                    ddcHelper::get<RDimY>(advection_field.span_cview()));


            // initialisation:
            ddc::for_each(grid, [&](IndexRP const irp) {
                feet_coords(irp) = CoordRP(ddc::coordinate(irp));
            });

            implicit_loop(advection_field, advection_field_coefs_k, feet_coords, dt / 2., tau);

            // Evaluate A^P at X^P:
            m_evaluator(
                    ddcHelper::get<RDimX>(advection_field_k).span_view(),
                    feet_coords.span_cview(),
                    ddcHelper::get<RDimX>(advection_field_coefs_k.span_cview()));
            m_evaluator(
                    ddcHelper::get<RDimY>(advection_field_k).span_view(),
                    feet_coords.span_cview(),
                    ddcHelper::get<RDimY>(advection_field_coefs_k.span_cview()));

            // Computed advection field (A^P(X^n) + A^P(X^P)) /2:
            ddc::for_each(grid, [&](IndexRP const irp) {
                ddcHelper::get<RDimX>(advection_field_k_tot)(irp)
                        = (ddcHelper::get<RDimX>(advection_field)(irp)
                           + ddcHelper::get<RDimX>(advection_field_k)(irp))
                          / 2.;
                ddcHelper::get<RDimY>(advection_field_k_tot)(irp)
                        = (ddcHelper::get<RDimY>(advection_field)(irp)
                           + ddcHelper::get<RDimY>(advection_field_k)(irp))
                          / 2.;
            });
            // X^k = X^n - dt * ( A^P(X^n) + A^P(X^P) )/2
            m_advection_solver(allfdistribu, advection_field_k_tot.span_cview(), dt);
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        Spline2D allfdistribu_coef(m_builder.spline_domain());
        m_builder(allfdistribu_coef.span_view(), allfdistribu.span_cview());
        PoissonLikeRHSFunction const
                charge_density_coord(allfdistribu_coef.span_cview(), m_evaluator);
        m_poisson_solver(charge_density_coord, coords, electrical_potential);

        ddc::PdiEvent("last_iteration")
                .with("iter", steps)
                .and_with("time", steps * dt)
                .and_with("density", allfdistribu)
                .and_with("electrical_potential", electrical_potential);

        end_time = std::chrono::system_clock::now();
        display_time_difference("Iterations time: ", start_time, end_time);



        return allfdistribu;
    }



private:
    double compute_square_polar_distance(CoordRP const& coord1, CoordRP const& coord2) const
    {
        CoordXY coord_xy1(m_mapping(coord1));
        CoordXY coord_xy2(m_mapping(coord2));

        const double x1 = ddc::select<RDimX>(coord_xy1);
        const double y1 = ddc::select<RDimY>(coord_xy1);
        const double x2 = ddc::select<RDimX>(coord_xy2);
        const double y2 = ddc::select<RDimY>(coord_xy2);

        return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    }


    void implicit_loop(
            VectorDSpanRP<RDimX, RDimY> advection_field,
            VectorSpline2DView<RDimX, RDimY> advection_field_coefs_k,
            SpanRP<CoordRP> feet_coords,
            double const dt,
            double const tau) const
    {
        IDomainRP const grid = advection_field.domain();
        VectorDFieldRP<RDimX, RDimY> advection_field_k(grid);
        VectorDFieldRP<RDimX, RDimY> advection_field_k_tot(grid);
        FieldRP<CoordRP> feet_coords_tmp(grid);

        double square_difference_feet = 0.;
        int count = 0;
        const int max_count = 50;
        do {
            count++;

            // Evaluate A at X^{k-1}:
            m_evaluator(
                    ddcHelper::get<RDimX>(advection_field_k).span_view(),
                    feet_coords.span_cview(),
                    ddcHelper::get<RDimX>(advection_field_coefs_k));
            m_evaluator(
                    ddcHelper::get<RDimY>(advection_field_k).span_view(),
                    feet_coords.span_cview(),
                    ddcHelper::get<RDimY>(advection_field_coefs_k));

            // Compute the new advection field A(X^n) + A(X^{k-1}):
            ddc::for_each(grid, [&](IndexRP const irp) {
                ddcHelper::get<RDimX>(advection_field_k_tot)(irp)
                        = ddcHelper::get<RDimX>(advection_field)(irp)
                          + ddcHelper::get<RDimX>(advection_field_k)(irp);
                ddcHelper::get<RDimY>(advection_field_k_tot)(irp)
                        = ddcHelper::get<RDimY>(advection_field)(irp)
                          + ddcHelper::get<RDimY>(advection_field_k)(irp);
            });

            // X^{k-1} = X^k:
            ddc::parallel_deepcopy(feet_coords_tmp, feet_coords);

            // X^k = X^n - dt* X^k:
            ddc::for_each(grid, [&](IndexRP const irp) {
                feet_coords(irp) = CoordRP(ddc::coordinate(irp));
            });
            m_foot_finder(feet_coords, advection_field_k_tot.span_cview(), dt);


            // Convergence test:
            square_difference_feet = 0.;
            ddc::for_each(grid, [&](IndexRP const irp) {
                double sqr_diff_feet
                        = compute_square_polar_distance(feet_coords(irp), feet_coords_tmp(irp));
                square_difference_feet = square_difference_feet > sqr_diff_feet
                                                 ? square_difference_feet
                                                 : sqr_diff_feet;
            });

        } while ((square_difference_feet > tau * tau) and (count < max_count));
    }
};
