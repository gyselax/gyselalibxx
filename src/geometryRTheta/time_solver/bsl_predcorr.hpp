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
#include "bsl_advection_rp.hpp"
#include "geometry.hpp"
#include "ifoot_finder.hpp"
#include "itimesolver.hpp"
#include "poisson_rhs_function.hpp"
#include "polarpoissonsolver.hpp"
#include "rk2.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "vlasovpoissonsolver.hpp"


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
 * - 1. & 2. From @f$\rho^n@f$, it computes @f$\phi^n@f$ and  @f$E^n@f$ with a VlasovPoissonSolver;
 * - 3. From @f$\rho^n@f$ and @f$A^n@f$, it computes @f$\rho^{n+1/2}@f$ with a BslAdvectionRP on @f$\frac{dt}{2}@f$;
 *
 * Secondly, it advects on a full time step:
 * - 4. & 5. From @f$\rho^{n+1/2}@f$, it computes @f$\phi^{n+1/2}@f$ and @f$E^{n+1/2}@f$ with a VlasovPoissonSolver;
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

    BslAdvectionRP<FootFinder> const& m_advection_solver;

    VlasovPoissonSolver<Mapping> const m_vlasov_poisson_solver;


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
            BslAdvectionRP<FootFinder> const& advection_solver,
            SplineRPBuilder const& builder,
            SplineRPEvaluator const& rhs_evaluator,
            PolarSplineFEMPoissonSolver const& poisson_solver)
        : m_advection_domain(advection_domain)
        , m_mapping(mapping)
        , m_advection_solver(advection_solver)
        , m_vlasov_poisson_solver(mapping, builder, rhs_evaluator, poisson_solver)
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

        DFieldRP electrical_potential0(grid);
        VectorDFieldRP<RDimX, RDimY> electric_field0(grid);
        m_vlasov_poisson_solver(electrical_potential0, electric_field0, allfdistribu);
        ddc::PdiEvent("iteration")
                .with("iter", 0)
                .and_with("time", 0)
                .and_with("density", allfdistribu)
                .and_with("electrical_potential", electrical_potential0);


        std::function<void(VectorDSpanRP<RDimX, RDimY>, DViewRP)> define_advection_field
                = [&](VectorDSpanRP<RDimX, RDimY> advection_field, DViewRP allfdistribu) {
                      // --- compute electric field:
                      DFieldRP electrical_potential(grid);
                      VectorDFieldRP<RDimX, RDimY> electric_field(grid);
                      m_vlasov_poisson_solver(electrical_potential, electric_field, allfdistribu);

                      // --- compute advection field from electric field:
                      ddc::for_each(advection_field.domain(), [&](IndexRP const idx) {
                          ddcHelper::get<RDimX>(advection_field)(idx)
                                  = -ddcHelper::get<RDimY>(electric_field)(idx);
                      });
                      ddc::deepcopy(
                              ddcHelper::get<RDimY>(advection_field),
                              ddcHelper::get<RDimX>(electric_field));
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
            VectorDFieldRP<RDimX, RDimY> electric_field(grid);
            m_vlasov_poisson_solver(electrical_potential, electric_field, allfdistribu);

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
