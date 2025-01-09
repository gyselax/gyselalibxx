// SPDX-License-Identifier: MIT

#pragma once
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "advection_domain.hpp"
#include "advection_field_rp.hpp"
#include "bsl_advection_rp.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
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
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes implicitly @f$\rho^P@f$ with a BslAdvectionRTheta on @f$ \frac{dt}{4} @f$:
 *      - the characteristic feet @f$X^P@f$ is such that @f$X^P = X^k@f$ with @f$X^k@f$ the result of the implicit method:
 *          - @f$ X^k = X^n - \frac{dt}{4} \partial_t X^k@f$.
 * 
 * Secondly, it corrects: 
 * - 4. From @f$\rho^P@f$, it computes @f$\phi^P@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 5. From @f$\phi^P@f$, it computes @f$A^P@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$ A^{P} @f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionRTheta on @f$ \frac{dt}{2} @f$.
 *      - the characteristic feet @f$X^C@f$ is such that @f$X^C = X^k@f$ with @f$X^k@f$ the result of the implicit method:
 *          - @f$\partial_t X^k = A^P(X^n) + A^P(X^{k-1}) @f$,
 *
 *
 * @tparam Mapping
 *      A class describing a mapping from curvilinear coordinates to cartesian coordinates.
 * @tparam AdvectionDomain
 *      An AdvectionDomain class.
 * @tparam FootFinder
 *      A IFootFinder class.
 *
 */
template <class Mapping, class AdvectionDomain>
class BslImplicitPredCorrRTheta : public ITimeSolverRTheta
{
private:
    using EulerMethod
            = Euler<host_t<FieldMemRTheta<CoordRTheta>>,
                    host_t<DVectorFieldMemRTheta<X, Y>>,
                    Kokkos::DefaultHostExecutionSpace>;

    using SplineFootFinderType = SplineFootFinder<EulerMethod, AdvectionDomain, Mapping>;

    Mapping const& m_mapping;

    BslAdvectionRTheta<SplineFootFinderType, Mapping> const& m_advection_solver;

    EulerMethod const m_euler;
    SplineFootFinderType const m_foot_finder;

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
     * @param[in] advection_domain
     *      An AdvectionDomain object which gives the information
     *      in which index range we advect.
     * @param[in] mapping
     *      The mapping function from the logical index range to the
     *      physical index range.
     * @param[in] advection_solver
     *      The advection operator with an Euler method.
     * @param[in] grid
     *      The index range on which the functions are defined.
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
    BslImplicitPredCorrRTheta(
            AdvectionDomain const& advection_domain,
            Mapping const& mapping,
            BslAdvectionRTheta<SplineFootFinderType, Mapping> const& advection_solver,
            IdxRangeRTheta const& grid,
            SplineRThetaBuilder_host const& builder,
            SplineRThetaEvaluatorNullBound_host const& rhs_evaluator,
            PolarSplineFEMPoissonLikeSolver<
                    GridR,
                    GridTheta,
                    PolarBSplinesRTheta,
                    SplineRThetaEvaluatorNullBound> const& poisson_solver,
            SplineRThetaEvaluatorConstBound_host const& advection_evaluator)
        : m_mapping(mapping)
        , m_advection_solver(advection_solver)
        , m_euler(grid)
        , m_foot_finder(m_euler, advection_domain, mapping, builder, advection_evaluator)
        , m_poisson_solver(poisson_solver)
        , m_builder(builder)
        , m_evaluator(advection_evaluator)

    {
    }



    host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> allfdistribu,
            double const dt,
            int const steps) const final
    {
        std::chrono::time_point<std::chrono::system_clock> start_time;
        std::chrono::time_point<std::chrono::system_clock> end_time;

        // Grid. ------------------------------------------------------------------------------------------
        IdxRangeRTheta const grid(get_idx_range<GridR, GridTheta>(allfdistribu));

        host_t<FieldMemRTheta<CoordRTheta>> coords(grid);
        ddc::for_each(grid, [&](IdxRTheta const irp) { coords(irp) = ddc::coordinate(irp); });

        IdxRangeBSR radial_bsplines(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                IdxStep<BSplinesR> {PolarBSplinesRTheta::continuity + 1}));
        IdxRangeBSTheta polar_idx_range(ddc::discrete_space<BSplinesTheta>().full_domain());

        // --- Electrostatic potential (phi). -------------------------------------------------------------
        DFieldMemRTheta electrical_potential(grid);
        host_t<DFieldMemRTheta> electrical_potential_host(grid);

        host_t<SplinePolar> electrostatic_potential_coef(
                PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>(),
                IdxRangeBSRTheta(radial_bsplines, polar_idx_range));

        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule, Kokkos::HostSpace>
                polar_spline_evaluator(extrapolation_rule);

        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        host_t<DVectorFieldMemRTheta<X, Y>> electric_field_alloc(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_alloc(grid);
        host_t<DVectorFieldRTheta<X, Y>> electric_field(electric_field_alloc);
        host_t<DVectorFieldRTheta<X, Y>> advection_field(advection_field_alloc);

        AdvectionFieldFinder advection_field_computer(m_mapping);

        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            host_t<Spline2DMem> allfdistribu_coef(get_spline_idx_range(m_builder));
            m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
            PoissonLikeRHSFunction const
                    charge_density_coord_1(get_const_field(allfdistribu_coef), m_evaluator);
            m_poisson_solver(charge_density_coord_1, electrostatic_potential_coef);

            polar_spline_evaluator(
                    get_field(electrical_potential_host),
                    get_const_field(coords),
                    electrostatic_potential_coef);

            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .and_with("time", iter * dt)
                    .and_with("density", allfdistribu)
                    .and_with("electrical_potential", electrical_potential_host);


            // STEP 2: From phi^n, we compute A^n:
            advection_field_computer(electrostatic_potential_coef, advection_field);


            // STEP 3: From rho^n and A^n, we compute rho^P: Vlasov equation
            host_t<DVectorFieldMemRTheta<X, Y>> advection_field_k(grid);
            host_t<DVectorFieldMemRTheta<X, Y>> advection_field_k_tot(grid);

            host_t<VectorSplineCoeffsMem2D<X, Y>> advection_field_coefs_k(
                    get_spline_idx_range(m_builder));
            m_builder(
                    ddcHelper::get<X>(advection_field_coefs_k),
                    ddcHelper::get<X>(get_const_field(advection_field)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs_k),
                    ddcHelper::get<Y>(get_const_field(advection_field)));

            host_t<FieldMemRTheta<CoordRTheta>> feet_coords(grid);
            host_t<FieldMemRTheta<CoordRTheta>> feet_coords_tmp(grid);


            // initialisation:
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                feet_coords(irp) = CoordRTheta(ddc::coordinate(irp));
            });

            const double tau = 1e-6;
            implicit_loop(
                    advection_field,
                    get_const_field(advection_field_coefs_k),
                    get_field(feet_coords),
                    dt / 4.,
                    tau);

            // Evaluate A^n at X^P:
            m_evaluator(
                    get_field(ddcHelper::get<X>(advection_field_k)),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs_k)));
            m_evaluator(
                    get_field(ddcHelper::get<Y>(advection_field_k)),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs_k)));

            // Compute the new advection field (E^n(X^n) + E^n(X^P)) /2:
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                ddcHelper::get<X>(advection_field_k_tot)(irp)
                        = (ddcHelper::get<X>(advection_field)(irp)
                           + ddcHelper::get<X>(advection_field_k)(irp))
                          / 2.;
                ddcHelper::get<Y>(advection_field_k_tot)(irp)
                        = (ddcHelper::get<Y>(advection_field)(irp)
                           + ddcHelper::get<Y>(advection_field_k)(irp))
                          / 2.;
            });


            // X^P = X^n - dt/2 * ( E^n(X^n) + E^n(X^P) )/2:
            // --- Copy phi^n because it will be modified:
            host_t<DFieldMemRTheta> allfdistribu_predicted(grid);
            ddc::parallel_deepcopy(allfdistribu_predicted, allfdistribu);
            m_advection_solver(
                    get_field(allfdistribu_predicted),
                    get_const_field(advection_field_k_tot),
                    dt / 2.);

            // --- advect also the feet because it is needed for the next step
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                feet_coords(irp) = CoordRTheta(ddc::coordinate(irp));
            });
            m_foot_finder(get_field(feet_coords), get_const_field(advection_field_k_tot), dt / 2.);


            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
            PoissonLikeRHSFunction const
                    charge_density_coord_4(get_const_field(allfdistribu_coef), m_evaluator);
            m_poisson_solver(charge_density_coord_4, electrostatic_potential_coef);

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(electrostatic_potential_coef, advection_field);


            // STEP 6: From rho^n and A^P, we compute rho^{n+1}: Vlasov equation
            m_builder(
                    ddcHelper::get<X>(advection_field_coefs_k),
                    ddcHelper::get<X>(get_const_field(advection_field)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs_k),
                    ddcHelper::get<Y>(get_const_field(advection_field)));


            // initialisation:
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                feet_coords(irp) = CoordRTheta(ddc::coordinate(irp));
            });

            implicit_loop(
                    advection_field,
                    get_const_field(advection_field_coefs_k),
                    get_field(feet_coords),
                    dt / 2.,
                    tau);

            // Evaluate A^P at X^P:
            m_evaluator(
                    get_field(ddcHelper::get<X>(advection_field_k)),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs_k)));
            m_evaluator(
                    get_field(ddcHelper::get<Y>(advection_field_k)),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs_k)));

            // Computed advection field (A^P(X^n) + A^P(X^P)) /2:
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                ddcHelper::get<X>(advection_field_k_tot)(irp)
                        = (ddcHelper::get<X>(advection_field)(irp)
                           + ddcHelper::get<X>(advection_field_k)(irp))
                          / 2.;
                ddcHelper::get<Y>(advection_field_k_tot)(irp)
                        = (ddcHelper::get<Y>(advection_field)(irp)
                           + ddcHelper::get<Y>(advection_field_k)(irp))
                          / 2.;
            });
            // X^k = X^n - dt * ( A^P(X^n) + A^P(X^P) )/2
            m_advection_solver(allfdistribu, get_const_field(advection_field_k_tot), dt);
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        host_t<Spline2DMem> allfdistribu_coef(get_spline_idx_range(m_builder));
        m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
        PoissonLikeRHSFunction const
                charge_density_coord(get_const_field(allfdistribu_coef), m_evaluator);
        m_poisson_solver(charge_density_coord, get_field(electrical_potential));
        ddc::parallel_deepcopy(electrical_potential_host, electrical_potential);

        ddc::PdiEvent("last_iteration")
                .with("iter", steps)
                .and_with("time", steps * dt)
                .and_with("density", allfdistribu)
                .and_with("electrical_potential", electrical_potential_host);

        end_time = std::chrono::system_clock::now();
        display_time_difference("Iterations time: ", start_time, end_time);



        return allfdistribu;
    }



private:
    double compute_square_polar_distance(CoordRTheta const& coord1, CoordRTheta const& coord2) const
    {
        CoordXY coord_xy1(m_mapping(coord1));
        CoordXY coord_xy2(m_mapping(coord2));

        const double x1 = ddc::select<X>(coord_xy1);
        const double y1 = ddc::select<Y>(coord_xy1);
        const double x2 = ddc::select<X>(coord_xy2);
        const double y2 = ddc::select<Y>(coord_xy2);

        return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    }


    void implicit_loop(
            host_t<DVectorFieldRTheta<X, Y>> advection_field,
            host_t<ConstVectorSplineCoeffs2D<X, Y>> advection_field_coefs_k,
            host_t<FieldRTheta<CoordRTheta>> feet_coords,
            double const dt,
            double const tau) const
    {
        IdxRangeRTheta const grid = get_idx_range(advection_field);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_k(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_k_tot(grid);
        host_t<FieldMemRTheta<CoordRTheta>> feet_coords_tmp(grid);

        double square_difference_feet = 0.;
        int count = 0;
        const int max_count = 50;
        do {
            count++;

            // Evaluate A at X^{k-1}:
            m_evaluator(
                    get_field(ddcHelper::get<X>(advection_field_k)),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(advection_field_coefs_k));
            m_evaluator(
                    get_field(ddcHelper::get<Y>(advection_field_k)),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(advection_field_coefs_k));

            // Compute the new advection field A(X^n) + A(X^{k-1}):
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                ddcHelper::get<X>(advection_field_k_tot)(irp)
                        = ddcHelper::get<X>(advection_field)(irp)
                          + ddcHelper::get<X>(advection_field_k)(irp);
                ddcHelper::get<Y>(advection_field_k_tot)(irp)
                        = ddcHelper::get<Y>(advection_field)(irp)
                          + ddcHelper::get<Y>(advection_field_k)(irp);
            });

            // X^{k-1} = X^k:
            ddc::parallel_deepcopy(feet_coords_tmp, feet_coords);

            // X^k = X^n - dt* X^k:
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                feet_coords(irp) = CoordRTheta(ddc::coordinate(irp));
            });
            m_foot_finder(feet_coords, get_const_field(advection_field_k_tot), dt);


            // Convergence test:
            square_difference_feet = 0.;
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                double sqr_diff_feet
                        = compute_square_polar_distance(feet_coords(irp), feet_coords_tmp(irp));
                square_difference_feet = square_difference_feet > sqr_diff_feet
                                                 ? square_difference_feet
                                                 : sqr_diff_feet;
            });

        } while ((square_difference_feet > tau * tau) and (count < max_count));
    }
};
