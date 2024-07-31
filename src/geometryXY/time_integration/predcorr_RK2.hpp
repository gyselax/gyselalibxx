// SPDX-License-Identifier: MIT

#pragma once
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include "bsl_advection_1d.hpp"
#include "ddc_aliases.hpp"
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
    void operator()(DFieldXY allfdistribu, double const dt, int const nbiter)
    {
        // Domain
        IdxRangeXY const meshXY = get_idx_range(allfdistribu);

        // Output of the Poisson solver
        DFieldMemXY electrostatic_potential_alloc(meshXY);
        DFieldXY electrostatic_potential = get_field(electrostatic_potential_alloc);

        VectorFieldXY_XY electric_field_alloc(meshXY);
        VectorSpanXY_XY electric_field = get_field(electric_field_alloc);

        // Definition of the RK2
        RK2<DFieldMemXY, VectorFieldXY_XY> predictor_corrector(meshXY);

        // Computation of the advection field: Poisson equation ---
        std::function<void(VectorSpanXY_XY, DConstFieldXY)> define_electric_field
                = [&](VectorSpanXY_XY electric_field, DConstFieldXY allfdistribu_view) {
                      IdxRangeXY xy_dom(get_idx_range<GridX, GridY>(allfdistribu_view));

                      // --- compute electrostatic potential and electric field:
                      DFieldMemXY electrostatic_potential_alloc(xy_dom);
                      DFieldXY electrostatic_potential = get_field(electrostatic_potential_alloc);

                      /*
                        The applied Poisson solver needs a span type for allfdistribu. 
                        The time stepper uses a view type. So we create a temporary  
                        allfdistribu_alloc chunk containing the values of the view 
                        allfdistribu to solve the type conflict in the Poisson solver. 
                    */
                      DFieldMemXY allfdistribu_alloc(xy_dom);
                      DFieldXY allfdistribu_span = get_field(allfdistribu_alloc);
                      ddc::parallel_deepcopy(
                              Kokkos::DefaultExecutionSpace(),
                              allfdistribu_span,
                              allfdistribu_view);


                      m_poisson_solver(electrostatic_potential, electric_field, allfdistribu_span);
                  };

        // Advection operator ---
        std::function<void(DFieldXY, VectorViewXY_XY, double)> advect_allfdistribu
                = [&](DFieldXY allfdistribu, VectorViewXY_XY electric_field, double dt) {
                      DConstFieldXY electric_field_x(ddcHelper::get<X>(electric_field));
                      DConstFieldXY electric_field_y(ddcHelper::get<Y>(electric_field));

                      // --- compute advection field:
                      IdxRangeXY dom = get_idx_range(electric_field);
                      DFieldMemXY advection_field_x_alloc(dom);
                      DFieldMemXY advection_field_y_alloc(dom);
                      DFieldXY advection_field_x = get_field(advection_field_x_alloc);
                      DFieldXY advection_field_y = get_field(advection_field_y_alloc);
                      ddc::parallel_for_each(
                              Kokkos::DefaultExecutionSpace(),
                              meshXY,
                              KOKKOS_LAMBDA(IdxXY const i_xy) {
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
                    = ddc::create_mirror_and_copy(ddcHelper::get<X>(electric_field));
            auto electric_field_y_host
                    = ddc::create_mirror_and_copy(ddcHelper::get<Y>(electric_field));
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