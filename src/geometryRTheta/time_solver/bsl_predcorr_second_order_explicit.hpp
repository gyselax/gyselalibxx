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
#include "euler.hpp"
#include "geometry.hpp"
#include "itimesolver.hpp"
#include "poisson_rhs_function.hpp"
#include "polarpoissonsolver.hpp"
#include "spline_foot_finder.hpp"
#include "spline_interpolator_2d_rp.hpp"



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
 * - 1. From @f$\rho^n@f$, it computes @f$\phi^n@f$ with a PolarSplineFEMPoissonSolver;
 * - 2. From @f$\phi^n@f$, it computes @f$A^n@f$ with a AdvectionFieldFinder;
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes @f$\rho^P@f$ with a BslAdvectionRP on @f$ dt @f$;
 *
 * We write @f$X^P@f$ the characteristic feet such that @f$\partial_t X^P = A^n(X^n)@f$.
 *
 * Secondly, it corrects:
 * - 4. From @f$\rho^P@f$, it computes @f$\phi^P@f$ with a PolarSplineFEMPoissonSolver;
 * - 5. From @f$\phi^P@f$, it computes @f$A^P@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$\frac{A^{P}(X^n) + A^n(X^P)}{2} @f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionRP on @f$ dt @f$.
 *
 * (With @f$X^C@f$ the characteristic feet such that @f$\partial_t X^C = \frac{A^{P}(X^n) + A^n(X^P)}{2} @f$.)
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
class BslExplicitPredCorrRP : public ITimeSolverRP
{
private:
    using EulerMethod = Euler<FieldRP<CoordRP>, VectorDFieldRP<RDimX, RDimY>>;


    Mapping const& m_mapping;

    BslAdvectionRP<SplineFootFinder<EulerMethod, AdvectionDomain>, Mapping> const&
            m_advection_solver;

    EulerMethod const m_euler;
    SplineFootFinder<EulerMethod, AdvectionDomain> const m_find_feet;

    PolarSplineFEMPoissonSolver const& m_poisson_solver;

    SplineRPBuilder const& m_builder;
    SplineRPEvaluatorConstBound const& m_evaluator;



public:
    /**
     * @brief Instantiate a BslExplicitPredCorrRP.
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
     *      advection field and the RHS.
     * @param[in] rhs_evaluator
     *      The evaluator of B-splines for the RHS.
     * @param[in] poisson_solver
     *      The Poisson solver which computes the electrical
     *      potential.
     * @param[in] advection_evaluator
     *      An evaluator of B-splines for the spline advection field.
     */
    BslExplicitPredCorrRP(
            AdvectionDomain const& advection_domain,
            Mapping const& mapping,
            BslAdvectionRP<SplineFootFinder<EulerMethod, AdvectionDomain>, Mapping>&
                    advection_solver,
            IDomainRP const& grid,
            SplineRPBuilder const& builder,
            SplineRPEvaluatorNullBound const& rhs_evaluator,
            PolarSplineFEMPoissonSolver const& poisson_solver,
            SplineRPEvaluatorConstBound const& advection_evaluator)
        : m_mapping(mapping)
        , m_advection_solver(advection_solver)
        , m_euler(grid)
        , m_find_feet(m_euler, advection_domain, builder, advection_evaluator)
        , m_poisson_solver(poisson_solver)
        , m_builder(builder)
        , m_evaluator(advection_evaluator)

    {
    }


    ~BslExplicitPredCorrRP() {};



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
                PolarBSplinesRP::singular_domain(),
                BSDomainRP(radial_bsplines, polar_domain));

        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<PolarBSplinesRP, ddc::NullExtrapolationRule> polar_spline_evaluator(
                extrapolation_rule);

        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        VectorDFieldRP<RDimX, RDimY> electric_field(grid);
        VectorDFieldRP<RDimX, RDimY> electric_field_predicted(grid);
        VectorDFieldRP<RDimX, RDimY> advection_field(grid);
        VectorDFieldRP<RDimX, RDimY> advection_field_predicted(grid);

        AdvectionFieldFinder advection_field_computer(m_mapping);



        // --- Parameter for linearisation of advection field: --------------------------------------------
        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            double const time = iter * dt;
            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            Spline2D allfdistribu_coef(m_builder.spline_domain());
            m_builder(allfdistribu_coef.span_view(), allfdistribu.span_cview());
            PoissonRHSFunction const
                    charge_density_coord_1(allfdistribu_coef.span_cview(), m_evaluator);
            m_poisson_solver(charge_density_coord_1, electrostatic_potential_coef);

            polar_spline_evaluator(
                    electrical_potential.span_view(),
                    coords.span_cview(),
                    electrostatic_potential_coef);

            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .and_with("time", time)
                    .and_with("density", allfdistribu)
                    .and_with("electrical_potential", electrical_potential);

            // STEP 2: From phi^n, we compute A^n:
            advection_field_computer(electrostatic_potential_coef, advection_field);


            // STEP 3: From rho^n and A^n, we compute rho^P: Vlasov equation
            // --- Copy rho^n because it will be modified:
            DFieldRP allfdistribu_predicted(grid);
            ddc::parallel_deepcopy(allfdistribu_predicted.span_view(), allfdistribu);
            m_advection_solver(allfdistribu_predicted.span_view(), advection_field.span_view(), dt);

            // --- advect also the feet because it is needed for the next step
            FieldRP<CoordRP> feet_coords(grid);
            ddc::for_each(grid, [&](IndexRP const irp) {
                feet_coords(irp) = CoordRP(ddc::coordinate(irp));
            });
            m_find_feet(feet_coords.span_view(), advection_field.span_view(), dt);


            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            m_builder(allfdistribu_coef.span_view(), allfdistribu.span_cview());
            PoissonRHSFunction const
                    charge_density_coord_4(allfdistribu_coef.span_cview(), m_evaluator);
            m_poisson_solver(charge_density_coord_4, electrostatic_potential_coef);

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(electrostatic_potential_coef, advection_field_predicted);


            // ---  we evaluate the advection field A^n at the characteristic feet X^P
            VectorDFieldRP<RDimX, RDimY> advection_field_evaluated(grid);
            VectorSpline2D<RDimX, RDimY> advection_field_coefs(m_builder.spline_domain());

            m_builder(
                    ddcHelper::get<RDimX>(advection_field_coefs),
                    ddcHelper::get<RDimX>(advection_field.span_cview()));
            m_builder(
                    ddcHelper::get<RDimY>(advection_field_coefs),
                    ddcHelper::get<RDimY>(advection_field.span_cview()));

            m_evaluator(
                    ddcHelper::get<RDimX>(advection_field_evaluated).span_view(),
                    feet_coords.span_cview(),
                    ddcHelper::get<RDimX>(advection_field_coefs.span_cview()));
            m_evaluator(
                    ddcHelper::get<RDimY>(advection_field_evaluated).span_view(),
                    feet_coords.span_cview(),
                    ddcHelper::get<RDimY>(advection_field_coefs.span_cview()));


            // STEP 6: From rho^n and (A^n(X^P) + A^P(X^n))/2, we compute rho^{n+1}: Vlasov equation
            ddc::for_each(grid, [&](IndexRP const irp) {
                ddcHelper::get<RDimX>(advection_field)(irp)
                        = (ddcHelper::get<RDimX>(advection_field_evaluated)(irp)
                           + ddcHelper::get<RDimX>(advection_field_predicted)(irp))
                          / 2.;
                ddcHelper::get<RDimY>(advection_field)(irp)
                        = (ddcHelper::get<RDimY>(advection_field_evaluated)(irp)
                           + ddcHelper::get<RDimY>(advection_field_predicted)(irp))
                          / 2.;
            });


            m_advection_solver(allfdistribu, advection_field.span_view(), dt);
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        Spline2D allfdistribu_coef(m_builder.spline_domain());
        m_builder(allfdistribu_coef.span_view(), allfdistribu.span_cview());
        PoissonRHSFunction const charge_density_coord(allfdistribu_coef.span_cview(), m_evaluator);
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
};
