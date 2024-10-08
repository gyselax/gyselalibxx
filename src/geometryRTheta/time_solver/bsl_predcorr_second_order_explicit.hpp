// SPDX-License-Identifier: MIT

#pragma once
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include "advection_domain.hpp"
#include "advection_field_rp.hpp"
#include "bsl_advection_rp.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "itimesolver.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
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
 * - 1. From @f$\rho^n@f$, it computes @f$\phi^n@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 2. From @f$\phi^n@f$, it computes @f$A^n@f$ with a AdvectionFieldFinder;
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes @f$\rho^P@f$ with a BslAdvectionRTheta on @f$ dt @f$;
 *
 * We write @f$X^P@f$ the characteristic feet such that @f$\partial_t X^P = A^n(X^n)@f$.
 *
 * Secondly, it corrects:
 * - 4. From @f$\rho^P@f$, it computes @f$\phi^P@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 5. From @f$\phi^P@f$, it computes @f$A^P@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$\frac{A^{P}(X^n) + A^n(X^P)}{2} @f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionRTheta on @f$ dt @f$.
 *
 * (With @f$X^C@f$ the characteristic feet such that @f$\partial_t X^C = \frac{A^{P}(X^n) + A^n(X^P)}{2} @f$.)
 *
 * @tparam Mapping
 *      A Curvilinear2DToCartesian class or one of its child classes.
 * @tparam IdxRangeAdvection
 *      An IdxRangeAdvection class.
 * @tparam FootFinder
 *      A IFootFinder class.
 *
 */
template <class Mapping, class IdxRangeAdvection>
class BslExplicitPredCorrRTheta : public ITimeSolverRTheta
{
private:
    using EulerMethod
            = Euler<host_t<FieldMemRTheta<CoordRTheta>>, host_t<DVectorFieldMemRTheta<X, Y>>>;


    Mapping const& m_mapping;

    BslAdvectionRTheta<SplineFootFinder<EulerMethod, IdxRangeAdvection>, Mapping> const&
            m_advection_solver;

    EulerMethod const m_euler;
    SplineFootFinder<EulerMethod, IdxRangeAdvection> const m_find_feet;

    PolarSplineFEMPoissonLikeSolver const& m_poisson_solver;

    SplineRThetaBuilder const& m_builder;
    SplineRThetaEvaluatorConstBound const& m_evaluator;



public:
    /**
     * @brief Instantiate a BslExplicitPredCorrRTheta.
     *
     * @param[in] advection_idx_range
     *      An IdxRangeAdvection object which gives the information
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
     *      advection field and the RHS.
     * @param[in] rhs_evaluator
     *      The evaluator of B-splines for the RHS.
     * @param[in] poisson_solver
     *      The PDE solver which computes the electrical
     *      potential.
     * @param[in] advection_evaluator
     *      An evaluator of B-splines for the spline advection field.
     */
    BslExplicitPredCorrRTheta(
            IdxRangeAdvection const& advection_idx_range,
            Mapping const& mapping,
            BslAdvectionRTheta<SplineFootFinder<EulerMethod, IdxRangeAdvection>, Mapping>&
                    advection_solver,
            IdxRangeRTheta const& grid,
            SplineRThetaBuilder const& builder,
            SplineRThetaEvaluatorNullBound const& rhs_evaluator,
            PolarSplineFEMPoissonLikeSolver const& poisson_solver,
            SplineRThetaEvaluatorConstBound const& advection_evaluator)
        : m_mapping(mapping)
        , m_advection_solver(advection_solver)
        , m_euler(grid)
        , m_find_feet(m_euler, advection_idx_range, builder, advection_evaluator)
        , m_poisson_solver(poisson_solver)
        , m_builder(builder)
        , m_evaluator(advection_evaluator)

    {
    }


    ~BslExplicitPredCorrRTheta() {};



    host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> allfdistribu,
            double const dt,
            int const steps) const
    {
        std::chrono::time_point<std::chrono::system_clock> start_time
                = std::chrono::system_clock::now();
        std::chrono::time_point<std::chrono::system_clock> end_time;

        // Grid. ------------------------------------------------------------------------------------------
        IdxRangeRTheta const grid(get_idx_range<GridR, GridTheta>(allfdistribu));

        host_t<FieldMemRTheta<CoordRTheta>> coords(grid);
        ddc::for_each(grid, [&](IdxRTheta const irp) { coords(irp) = ddc::coordinate(irp); });

        IdxRangeBSR radial_bsplines(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                IdxStep<BSplinesR> {PolarBSplinesRTheta::continuity + 1}));
        IdxRangeBSTheta polar_idx_range(ddc::discrete_space<BSplinesTheta>().full_domain());

        // --- Electrostatic potential (phi). -------------------------------------------------------------
        host_t<DFieldMemRTheta> electrical_potential(grid);

        SplinePolar electrostatic_potential_coef(
                PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>(),
                IdxRangeBSRTheta(radial_bsplines, polar_idx_range));

        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule>
                polar_spline_evaluator(extrapolation_rule);

        // --- For the computation of advection field from the electrostatic potential (phi): -------------
        host_t<DVectorFieldMemRTheta<X, Y>> electric_field(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> electric_field_predicted(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field(grid);
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_predicted(grid);

        AdvectionFieldFinder advection_field_computer(m_mapping);



        // --- Parameter for linearisation of advection field: --------------------------------------------
        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            double const time = iter * dt;
            // STEP 1: From rho^n, we compute phi^n: Poisson equation
            host_t<Spline2D> allfdistribu_coef(get_spline_idx_range(m_builder));
            m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
            PoissonLikeRHSFunction const
                    charge_density_coord_1(get_const_field(allfdistribu_coef), m_evaluator);
            m_poisson_solver(charge_density_coord_1, electrostatic_potential_coef);

            polar_spline_evaluator(
                    get_field(electrical_potential),
                    get_const_field(coords),
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
            host_t<DFieldMemRTheta> allfdistribu_predicted(grid);
            ddc::parallel_deepcopy(get_field(allfdistribu_predicted), allfdistribu);
            m_advection_solver(get_field(allfdistribu_predicted), get_field(advection_field), dt);

            // --- advect also the feet because it is needed for the next step
            host_t<FieldMemRTheta<CoordRTheta>> feet_coords(grid);
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                feet_coords(irp) = CoordRTheta(ddc::coordinate(irp));
            });
            m_find_feet(get_field(feet_coords), get_field(advection_field), dt);


            // STEP 4: From rho^P, we compute phi^P: Poisson equation
            m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
            PoissonLikeRHSFunction const
                    charge_density_coord_4(get_const_field(allfdistribu_coef), m_evaluator);
            m_poisson_solver(charge_density_coord_4, electrostatic_potential_coef);

            // STEP 5: From phi^P, we compute A^P:
            advection_field_computer(electrostatic_potential_coef, advection_field_predicted);


            // ---  we evaluate the advection field A^n at the characteristic feet X^Theta
            host_t<DVectorFieldMemRTheta<X, Y>> advection_field_evaluated(grid);
            host_t<VectorSplineCoeffsMem2D<X, Y>> advection_field_coefs(
                    get_spline_idx_range(m_builder));

            m_builder(
                    ddcHelper::get<X>(advection_field_coefs),
                    ddcHelper::get<X>(get_const_field(advection_field)));
            m_builder(
                    ddcHelper::get<Y>(advection_field_coefs),
                    ddcHelper::get<Y>(get_const_field(advection_field)));

            m_evaluator(
                    get_field(ddcHelper::get<X>(advection_field_evaluated)),
                    get_const_field(feet_coords),
                    ddcHelper::get<X>(get_const_field(advection_field_coefs)));
            m_evaluator(
                    get_field(ddcHelper::get<Y>(advection_field_evaluated)),
                    get_const_field(feet_coords),
                    ddcHelper::get<Y>(get_const_field(advection_field_coefs)));


            // STEP 6: From rho^n and (A^n(X^P) + A^P(X^n))/2, we compute rho^{n+1}: Vlasov equation
            ddc::for_each(grid, [&](IdxRTheta const irp) {
                ddcHelper::get<X>(advection_field)(irp)
                        = (ddcHelper::get<X>(advection_field_evaluated)(irp)
                           + ddcHelper::get<X>(advection_field_predicted)(irp))
                          / 2.;
                ddcHelper::get<Y>(advection_field)(irp)
                        = (ddcHelper::get<Y>(advection_field_evaluated)(irp)
                           + ddcHelper::get<Y>(advection_field_predicted)(irp))
                          / 2.;
            });


            m_advection_solver(allfdistribu, get_field(advection_field), dt);
        }

        // STEP 1: From rho^n, we compute phi^n: Poisson equation
        host_t<Spline2D> allfdistribu_coef(get_spline_idx_range(m_builder));
        m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
        PoissonLikeRHSFunction const
                charge_density_coord(get_const_field(allfdistribu_coef), m_evaluator);
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
