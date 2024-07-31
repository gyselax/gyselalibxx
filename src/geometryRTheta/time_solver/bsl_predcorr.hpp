// SPDX-License-Identifier: MIT

#pragma once
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <utils_tools.hpp>

#include "advection_domain.hpp"
#include "advection_field_rp.hpp"
#include "bsl_advection_rp.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "ifoot_finder.hpp"
#include "itimesolver.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "rk2.hpp"
#include "spline_interpolator_2d_rp.hpp"

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
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes @f$\rho^{n+1/2}@f$ with a BslAdvectionRTheta on @f$\frac{dt}{2}@f$;
 *
 * Secondly, it advects on a full time step:
 * - 4. From @f$\rho^{n+1/2}@f$, it computes @f$\phi^{n+1/2}@f$ with a PolarSplineFEMPoissonLikeSolver;
 * - 5. From @f$\phi^{n+1/2}@f$, it computes @f$A^{n+1/2}@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$A^{n+1/2}@f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionRTheta on @f$dt@f$.
 *
 * @tparam Mapping
 *      A Curvilinear2DToCartesian class or one of its child classes.
 * @tparam FootFinder
 *      A IFootFinder class.
 *
 */
template <class Mapping, class FootFinder>
class BslPredCorrRTheta : public ITimeSolverRTheta
{
private:
    Mapping const& m_mapping;

    BslAdvectionRTheta<FootFinder, Mapping> const& m_advection_solver;

    PolarSplineFEMPoissonLikeSolver const& m_poisson_solver;

    SplineRThetaBuilder const& m_builder;
    SplineRThetaEvaluatorNullBound const& m_spline_evaluator;


public:
    /**
     * @brief Instantiate a BslPredCorrRTheta.
     *
     * @param[in] mapping
     *      The mapping function from the logical index range to the
     *      physical index range.
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
            BslAdvectionRTheta<FootFinder, Mapping> const& advection_solver,
            SplineRThetaBuilder const& builder,
            SplineRThetaEvaluatorNullBound const& rhs_evaluator,
            PolarSplineFEMPoissonLikeSolver const& poisson_solver)
        : m_mapping(mapping)
        , m_advection_solver(advection_solver)
        , m_poisson_solver(poisson_solver)
        , m_builder(builder)
        , m_spline_evaluator(rhs_evaluator)
    {
    }


    ~BslPredCorrRTheta() {};


    DFieldRTheta operator()(DFieldRTheta allfdistribu, double const dt, int const steps) const
    {
        std::chrono::time_point<std::chrono::system_clock> start_time
                = std::chrono::system_clock::now();
        std::chrono::time_point<std::chrono::system_clock> end_time;


        // Grid. ------------------------------------------------------------------------------------------
        IdxRangeRTheta grid(get_idx_range<GridR, GridTheta>(allfdistribu));
        FieldMemRTheta<CoordRTheta> coords(grid);
        ddc::for_each(grid, [&](IdxRTheta const irp) { coords(irp) = ddc::coordinate(irp); });
        AdvectionFieldFinder advection_field_computer(m_mapping);

        BSIdxRangeR radial_bsplines(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                IdxStep<BSplinesR> {PolarBSplinesRTheta::continuity + 1}));
        BSIdxRangeTheta polar_idx_range(ddc::discrete_space<BSplinesTheta>().full_domain());

        SplinePolar electrostatic_potential_coef(
                PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>(),
                BSIdxRangeRTheta(radial_bsplines, polar_idx_range));
        ddc::NullExtrapolationRule extrapolation_rule;
        PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule>
                polar_spline_evaluator(extrapolation_rule);


        DFieldMemRTheta electrical_potential0(grid);

        Spline2D allfdistribu_coef(get_spline_idx_range(m_builder));
        m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
        PoissonLikeRHSFunction const
                charge_density_coord(get_const_field(allfdistribu_coef), m_spline_evaluator);
        m_poisson_solver(charge_density_coord, coords, electrical_potential0);

        ddc::PdiEvent("iteration")
                .with("iter", 0)
                .and_with("time", 0)
                .and_with("density", allfdistribu)
                .and_with("electrical_potential", electrical_potential0);


        std::function<void(DVectorFieldRTheta<X, Y>, DConstFieldRTheta)> define_advection_field
                = [&](DVectorFieldRTheta<X, Y> advection_field, DConstFieldRTheta allfdistribu) {
                      // --- compute electrostatic potential:
                      Spline2D allfdistribu_coef(get_spline_idx_range(m_builder));
                      m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
                      PoissonLikeRHSFunction const charge_density_coord(
                              get_const_field(allfdistribu_coef),
                              m_spline_evaluator);
                      m_poisson_solver(charge_density_coord, electrostatic_potential_coef);

                      // --- compute advection field:
                      advection_field_computer(electrostatic_potential_coef, advection_field);
                  };

        std::function<void(DFieldRTheta, DConstVectorFieldRTheta<X, Y>, double)> advect_allfdistribu
                = [&](DFieldRTheta allfdistribu,
                      DConstVectorFieldRTheta<X, Y> advection_field,
                      double dt) { m_advection_solver(allfdistribu, advection_field, dt); };

        RK2<DFieldMemRTheta, DVectorFieldMemRTheta<X, Y>> time_stepper(grid);

        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            time_stepper
                    .update(Kokkos::DefaultHostExecutionSpace(),
                            allfdistribu,
                            dt,
                            define_advection_field,
                            advect_allfdistribu);

            DFieldMemRTheta electrical_potential(grid);
            Spline2D allfdistribu_coef(get_spline_idx_range(m_builder));
            m_builder(get_field(allfdistribu_coef), get_const_field(allfdistribu));
            PoissonLikeRHSFunction const
                    charge_density_coord(get_const_field(allfdistribu_coef), m_spline_evaluator);
            m_poisson_solver(charge_density_coord, coords, electrical_potential);

            ddc::PdiEvent("iteration")
                    .with("iter", iter + 1)
                    .and_with("time", iter * dt)
                    .and_with("density", allfdistribu)
                    .and_with("electrical_potential", electrical_potential);
        }
        end_time = std::chrono::system_clock::now();


        display_time_difference("Iterations time: ", start_time, end_time);


        return allfdistribu;
    }
};
