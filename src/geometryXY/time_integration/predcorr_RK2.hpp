// SPDX-License-Identifier: MIT

#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include "bsl_advection_1d.hpp"
#include "directional_tag.hpp"
#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "rk2.hpp"
#include "utils_tools.hpp"
#include "vector_field.hpp"
#include "vector_field_span.hpp"

/**
 * @brief Predictor-corrector based on RK2 for the guiding-center model. 
 * 
 * It solves in time the following guiding-center equations system:
 *
 * - @f$  -\Delta \phi = f @f$,
 * - @f$ E = - \nabla \phi @f$,
 * - @f$ \partial_t f - E_y \partial_x f + E_x \partial_y f= 0@f$.
 *
 * This method is mainly a Runge-Kutta 2 method:
 *
 * for @f$ n \geq 0 @f$,
 *
 * First, it advects on a half time step:
 * - 1./2. From @f$f^n@f$, it computes @f$E^n@f$ with a FFTPoissonSolver;
 * - 3. From @f$f^n@f$ and @f$E^n@f$, it computes @f$f^{n+1/2}@f$ with a BslAdvection1D on @f$\frac{dt}{2}@f$;
 *
 * Secondly, it advects on a full time step:
 * - 4./5. From @f$f^{n+1/2}@f$, it computes @f$E^{n+1/2}@f$ with a FFTPoissonSolver;
 * - 6. From @f$f^n@f$ and @f$E^{n+1/2}@f$, it computes @f$f^{n+1}@f$ with a BslAdvectionRP on @f$dt@f$.
 * 
 * @tparam PoissonSolver Type of the Poisson solver applied in the method. 
 * @tparam AdvectionX Type of the 1D advection operator applied to advect along X. 
 * @tparam AdvectionY Type of the 1D advection operator applied to advect along Y. 
 */
template <class PoissonSolver, class AdvectionX, class AdvectionY>
class PredCorrRK2XY
{
private:
    PoissonSolver const& m_poisson_solver;

    AdvectionX const& m_advection_x;
    AdvectionY const& m_advection_y;


public:
    /**
     * @brief Instantiate the predictor-corrector.
     * @param poisson_solver Poisson solver also computing the electric field.  
     * @param advection_x 1D advection operator along @f$ x @f$ direction.
     * @param advection_y 1D advection operator along @f$ y @f$ direction.
     */
    PredCorrRK2XY(
            PoissonSolver const& poisson_solver,
            AdvectionX const& advection_x,
            AdvectionY const& advection_y)
        : m_poisson_solver(poisson_solver)
        , m_advection_x(advection_x)
        , m_advection_y(advection_y) {};

    ~PredCorrRK2XY() = default;

    /**
     * @brief Apply the predictor-corrector method on several time steps. 
     *  
     *  Along the simulation, the data are saved in an output folder. 
     * 
     * @param allfdistribu Initial function  @f$f (0, x, y)@f$.
     * @param dt Time step. 
     * @param nbiter Number of time steps. 
     */
    void operator()(DSpanXY allfdistribu, double const dt, int const nbiter)
    {
        // Domain
        IDomainXY const meshXY = allfdistribu.domain();

        // Output of the Poisson solver
        DFieldXY electrostatic_potential_alloc(meshXY);
        DSpanXY electrostatic_potential = electrostatic_potential_alloc.span_view();

        VectorFieldXY_XY electric_field_alloc(meshXY);
        VectorSpanXY_XY electric_field = electric_field_alloc.span_view();

        // Definition of the RK2
        RK2<DFieldXY, VectorFieldXY_XY> predictor_corrector(meshXY);

        // Computation of the advection field: Poisson equation ---
        std::function<void(VectorSpanXY_XY, DViewXY)> define_electric_field
                = [&](VectorSpanXY_XY electric_field, DViewXY allfdistribu_view) {
                      IDomainXY xy_dom(allfdistribu_view.domain<IDimX, IDimY>());

                      // --- compute electrostatic potential and electric field:
                      DFieldXY electrostatic_potential_alloc(xy_dom);
                      DSpanXY electrostatic_potential = electrostatic_potential_alloc.span_view();

                      /*
                        The applied Poisson solver needs a span type for allfdistribu. 
                        The time stepper uses a view type. So we create a temporary  
                        allfdistribu_alloc chunk containing the values of the view 
                        allfdistribu to solve the type conflict in the Poisson solver. 
                    */
                      DFieldXY allfdistribu_alloc(xy_dom);
                      DSpanXY allfdistribu_span = allfdistribu_alloc.span_view();
                      ddc::parallel_deepcopy(
                              Kokkos::DefaultExecutionSpace(),
                              allfdistribu_span,
                              allfdistribu_view);


                      m_poisson_solver(electrostatic_potential, electric_field, allfdistribu_span);
                  };

        // Advection operator ---
        std::function<void(DSpanXY, VectorViewXY_XY, double)> advect_allfdistribu
                = [&](DSpanXY allfdistribu, VectorViewXY_XY electric_field, double dt) {
                      DViewXY electric_field_x(ddcHelper::get<RDimX>(electric_field));
                      DViewXY electric_field_y(ddcHelper::get<RDimY>(electric_field));

                      // --- compute advection field:
                      IDomainXY dom = electric_field.domain();
                      DFieldXY advection_field_x_alloc(dom);
                      DFieldXY advection_field_y_alloc(dom);
                      DSpanXY advection_field_x = advection_field_x_alloc.span_view();
                      DSpanXY advection_field_y = advection_field_y_alloc.span_view();
                      ddc::parallel_for_each(
                              Kokkos::DefaultExecutionSpace(),
                              meshXY,
                              KOKKOS_LAMBDA(IndexXY const i_xy) {
                                  advection_field_x(i_xy) = -electric_field_y(i_xy);
                                  advection_field_y(i_xy) = electric_field_x(i_xy);
                              });


                      // --- Strang splitting for the advection
                      m_advection_x(allfdistribu, advection_field_x, dt / 2);
                      m_advection_y(allfdistribu, advection_field_y, dt);
                      m_advection_x(allfdistribu, advection_field_x, dt / 2);
                  };



        // Iteration on the number of steps.
        for (int iter(1); iter < nbiter + 1; ++iter) {
            predictor_corrector
                    .update(Kokkos::DefaultExecutionSpace(),
                            allfdistribu,
                            dt,
                            define_electric_field,
                            advect_allfdistribu);

            // Save the data ---
            m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);
            auto allfdistribu_host = ddc::create_mirror_and_copy(allfdistribu);
            auto electrostatic_potential_host
                    = ddc::create_mirror_and_copy(electrostatic_potential);
            auto electric_field_x_host
                    = ddc::create_mirror_and_copy(ddcHelper::get<RDimX>(electric_field));
            auto electric_field_y_host
                    = ddc::create_mirror_and_copy(ddcHelper::get<RDimY>(electric_field));
            ddc::PdiEvent("iteration")
                    .with("iter", iter)
                    .and_with("time_saved", iter * dt)
                    .and_with("fdistribu", allfdistribu_host)
                    .and_with("electrostatic_potential", electrostatic_potential_host)
                    .and_with("electric_field_x", electric_field_x_host)
                    .and_with("electric_field_y", electric_field_y_host);
        }
    };
};