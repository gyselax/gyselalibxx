// SPDX-License-Identifier: MIT

#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <utils_tools.hpp>

#include "sll/spline_evaluator_2d.hpp"

#include "advection_domain.hpp"
#include "advection_field_rp.hpp"
#include "bsl_advection_rp.hpp"
#include "geometry.hpp"
#include "ifoot_finder.hpp"
#include "itimesolver.hpp"
#include "poisson_rhs_function.hpp"
#include "polarpoissonsolver.hpp"
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
 * - 1. From @f$\rho^n@f$, it computes @f$\phi^n@f$ with a PolarSplineFEMPoissonSolver;
 * - 2. From @f$\phi^n@f$, it computes @f$A^n@f$ with a AdvectionFieldFinder;
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes @f$\rho^{n+1/2}@f$ with a BslAdvectionRP on @f$\frac{dt}{2}@f$;
 *
 * Secondly, it advects on a full time step:
 * - 4. From @f$\rho^{n+1/2}@f$, it computes @f$\phi^{n+1/2}@f$ with a PolarSplineFEMPoissonSolver;
 * - 5. From @f$\phi^{n+1/2}@f$, it computes @f$A^{n+1/2}@f$ with a AdvectionFieldFinder;
 * - 6. From @f$\rho^n@f$ and @f$A^{n+1/2}@f$, it computes @f$\rho^{n+1}@f$ with a BslAdvectionRP on @f$dt@f$.
 *
 * @tparam Mapping
 *      A Curvilinear2DToCartesian class or one of its child classes.
 * @tparam AdvectionDomain
 *      An AdvectionDomain class.
 * @tparam FootFinder
 *      A IFootFinder class.
 *
 */
template <class Mapping, class AdvectionDomain, class FootFinder>
class BslPredCorrRP : public ITimeSolverRP
{
private:
    using Evaluator = SplineEvaluator2D<BSplinesR, BSplinesP>;
    using Builder = SplineBuilder2D<SplineRBuilder, SplinePBuilder>;


    AdvectionDomain const& m_advection_domain;

    Mapping const& m_mapping;

    BslAdvectionRP<FootFinder, Mapping> const& m_advection_solver;

    PolarSplineFEMPoissonSolver const& m_poisson_solver;

    SplineRPBuilder const& m_builder;
    SplineRPEvaluator const& m_spline_evaluator;


public:
    /**
     * @brief Instantiate a BslPredCorrRP.
     *
     * @param[in] advection_domain
     *      An AdvectionDomain object which gives the information
     *      in which domain we advect.
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
     *      The Poisson solver which computes the electrical
     *      potential.
     */
    BslPredCorrRP(
            AdvectionDomain const& advection_domain,
            Mapping const& mapping,
            BslAdvectionRP<FootFinder, Mapping> const& advection_solver,
            SplineRPBuilder const& builder,
            SplineRPEvaluator const& rhs_evaluator,
            PolarSplineFEMPoissonSolver const& poisson_solver)
        : m_advection_domain(advection_domain)
        , m_mapping(mapping)
        , m_advection_solver(advection_solver)
        , m_poisson_solver(poisson_solver)
        , m_builder(builder)
        , m_spline_evaluator(rhs_evaluator)
    {
    }


    ~BslPredCorrRP() {};


    DSpanRP operator()(DSpanRP allfdistribu, double const dt, int const steps) const
    {
        std::chrono::time_point<std::chrono::system_clock> start_time
                = std::chrono::system_clock::now();
        std::chrono::time_point<std::chrono::system_clock> end_time;


        // Grid. ------------------------------------------------------------------------------------------
        IDomainRP grid(allfdistribu.domain<IDimR, IDimP>());
        FieldRP<CoordRP> coords(grid);
        ddc::for_each(grid, [&](IndexRP const irp) { coords(irp) = ddc::coordinate(irp); });
        AdvectionFieldFinder advection_field_computer(m_mapping);

        BSDomainR radial_bsplines(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                ddc::DiscreteVector<BSplinesR> {PolarBSplinesRP::continuity + 1}));
        BSDomainP polar_domain(ddc::discrete_space<BSplinesP>().full_domain());

        SplinePolar electrostatic_potential_coef(
                PolarBSplinesRP::singular_domain(),
                BSDomainRP(radial_bsplines, polar_domain));
        PolarSplineEvaluator<PolarBSplinesRP> polar_spline_evaluator(
                g_polar_null_boundary_2d<PolarBSplinesRP>);


        DFieldRP electrical_potential0(grid);

        Spline2D allfdistribu_coef(m_builder.spline_domain());
        m_builder(allfdistribu_coef, allfdistribu);
        PoissonRHSFunction const charge_density_coord(allfdistribu_coef, m_spline_evaluator);
        m_poisson_solver(charge_density_coord, coords, electrical_potential0);

        ddc::PdiEvent("iteration")
                .with("iter", 0)
                .and_with("time", 0)
                .and_with("density", allfdistribu)
                .and_with("electrical_potential", electrical_potential0);


        std::function<void(VectorDSpanRP<RDimX, RDimY>, DViewRP)> define_advection_field
                = [&](VectorDSpanRP<RDimX, RDimY> advection_field, DViewRP allfdistribu) {
                      // --- compute electrostatic potential:
                      Spline2D allfdistribu_coef(m_builder.spline_domain());
                      m_builder(allfdistribu_coef, allfdistribu);
                      PoissonRHSFunction const
                              charge_density_coord(allfdistribu_coef, m_spline_evaluator);
                      m_poisson_solver(charge_density_coord, electrostatic_potential_coef);

                      // --- compute advection field:
                      advection_field_computer(electrostatic_potential_coef, advection_field);
                  };

        std::function<void(DSpanRP, VectorDViewRP<RDimX, RDimY>, double)> advect_allfdistribu =
                [&](DSpanRP allfdistribu, VectorDViewRP<RDimX, RDimY> advection_field, double dt) {
                    m_advection_solver(allfdistribu, advection_field, dt);
                };

        RK2<DFieldRP, VectorDFieldRP<RDimX, RDimY>> time_stepper(grid);

        start_time = std::chrono::system_clock::now();
        for (int iter(0); iter < steps; ++iter) {
            time_stepper
                    .update(Kokkos::Serial(),
                            allfdistribu,
                            dt,
                            define_advection_field,
                            advect_allfdistribu);

            DFieldRP electrical_potential(grid);
            Spline2D allfdistribu_coef(m_builder.spline_domain());
            m_builder(allfdistribu_coef, allfdistribu);
            PoissonRHSFunction const charge_density_coord(allfdistribu_coef, m_spline_evaluator);
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
