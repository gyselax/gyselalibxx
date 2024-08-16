// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "utils_tools.hpp"


/**
 * @brief Equilibrium solution of a Vlasov-Poissson equations system
 * in polar coordinates.
 *
 * @tparam Mapping
 *      A Curvilinear2DToCartesian mapping function class.
 */
template <class Mapping>
class VortexMergerEquilibria
{
private:
    Mapping const& m_mapping;
    IdxRangeRTheta const& m_grid;
    SplineRThetaBuilder const& m_builder;
    SplineRThetaEvaluatorNullBound const& m_evaluator;
    PolarSplineFEMPoissonLikeSolver const& m_poisson_solver;

public:
    /**
     * @brief Instantiate a VortexMergerEquilibria.
     *
     * @param[in] mapping
     *      The mapping function from the logical index range
     *      to the physical index range.
     * @param[in] grid
     *      The index range where the equilibrium is defined.
     * @param[in] builder
     *      A spline builder to get the spline representation
     *      of the RHS of the PDE.
     * @param[in] evaluator
     *      The evaluator of B-splines for the RHS of the
     *      PDE.
     * @param[in] poisson_solver
     *      The PDE solver which computes the electrical potential.
     */
    VortexMergerEquilibria(
            Mapping const& mapping,
            IdxRangeRTheta const& grid,
            SplineRThetaBuilder const& builder,
            SplineRThetaEvaluatorNullBound const& evaluator,
            PolarSplineFEMPoissonLikeSolver const& poisson_solver)
        : m_mapping(mapping)
        , m_grid(grid)
        , m_builder(builder)
        , m_evaluator(evaluator)
        , m_poisson_solver(poisson_solver) {};


    ~VortexMergerEquilibria() {};

    /**
     * @brief Get an equilibrium.
     *
     * The equilibrium is determined by the eigenvalue problem.
     * For the given initial data @f$ (\sigma^0, \phi^0) @f$,
     *
     * 1. compute @f$ \rho^i = \sigma^{i-1} f(\phi^{i-1})@f$;
     * 2. compute @f$ \phi_*^i@f$ with @f$ -\nabla\cdot\nabla\phi_*^i = \rho^i@f$;
     * 3. for @f$ \phi_{\text{max}}@f$ given, compute @f$ c^i@f$ with @f$ c^i = \phi_{\text{max}} / \Vert \phi_*^i \Vert_{\mathcal{L}^\infty} @f$ ;
     * 4. compute @f$ (\sigma^i, \phi^i) = c^i (\sigma^{i-1}, \phi_*^i) @f$.
     *
     * We iterate until @f$ |\sigma^i - \sigma^{i-1}|\leq \tau @f$.
     *
     * For the vortex merger simulation, @f$ f(\phi) = \phi^2@f$.
     *
     * The algorithm is also detailed in Edoardo Zoni's article
     * (https://doi.org/10.1016/j.jcp.2019.108889).
     *
     * @param[in] sigma
     *      The @f$\sigma @f$ parameter.
     * @param[in] phi_eq
     *      The equilibrium electrical potential @f$ \phi @f$.
     * @param[out] rho_eq
     *      The equilibrium density  @f$ \rho @f$.
     * @param[in] function
     *      The function @f$ f @f$.
     * @param[in] phi_max
     *      The maximal value of the electrical potential @f$ \phi @f$.
     * @param[in] tau
     *      The @f$ \tau @f$ parameter.
     * @param[in] count_max
     *       The maximal number of iteration of the
     *       implicit loop.
     */
    void find_equilibrium(
            DFieldRTheta sigma,
            DFieldRTheta phi_eq,
            DFieldRTheta rho_eq,
            std::function<double(double const)> const& function,
            double const phi_max, // ToDo: ADD CASE WHERE RHO_MAX IS GIVEN.
            double const tau,
            int count_max = 25) const
    {
        DFieldMemRTheta phi_star(m_grid);
        DFieldMemRTheta ci(m_grid);

        auto idx_range_bsplinesRTheta = get_spline_idx_range(m_builder);
        Spline2D rho_coef(idx_range_bsplinesRTheta);

        FieldMemRTheta<CoordRTheta> coords(m_grid);
        ddc::for_each(m_grid, [&](IdxRTheta const irp) { coords(irp) = ddc::coordinate(irp); });

        double difference_sigma(0.);
        int count = 0;

        do {
            count += 1;
            // STEP 1: compute rho^i
            ddc::for_each(m_grid, [&](IdxRTheta const irp) {
                rho_eq(irp) = sigma(irp) * function(phi_eq(irp));
            });


            // STEP 2: compute phi_star^i with PDE solver
            m_builder(get_field(rho_coef), get_const_field(rho_eq));
            PoissonLikeRHSFunction poisson_rhs(rho_coef, m_evaluator);
            m_poisson_solver(poisson_rhs, get_const_field(coords), get_field(phi_star));

            // STEP 3: compute c^i
            // If phi_max is given:
            double norm_Linf_phi_star(0.);
            ddc::for_each(m_grid, [&](IdxRTheta const irp) {
                double const abs_phi_star = fabs(phi_star(irp));
                norm_Linf_phi_star
                        = norm_Linf_phi_star > abs_phi_star ? norm_Linf_phi_star : abs_phi_star;
            });

            ddc::for_each(m_grid, [&](IdxRTheta const irp) {
                ci(irp) = phi_max / norm_Linf_phi_star;
            });


            // STEP 4: update sigma and phi
            difference_sigma = 0.;
            ddc::for_each(m_grid, [&](IdxRTheta const irp) {
                double const abs_diff_sigma = fabs(sigma(irp) - ci(irp) * sigma(irp));
                difference_sigma
                        = difference_sigma > abs_diff_sigma ? difference_sigma : abs_diff_sigma;

                sigma(irp) = ci(irp) * sigma(irp);
                phi_eq(irp) = ci(irp) * phi_star(irp);
            });

        } while ((difference_sigma > tau) and (count < count_max));


        // STEP 1: compute rho^i
        ddc::for_each(m_grid, [&](IdxRTheta const irp) {
            rho_eq(irp) = sigma(irp) * function(phi_eq(irp));
        });

        // Unify at the center point:
        auto r_idx_range = get_idx_range<GridR>(rho_eq);
        auto theta_idx_range = get_idx_range<GridTheta>(rho_eq);
        if (std::fabs(ddc::coordinate(r_idx_range.front())) < 1e-15) {
            ddc::for_each(theta_idx_range, [&](const IdxTheta ip) {
                rho_eq(r_idx_range.front(), ip)
                        = rho_eq(r_idx_range.front(), theta_idx_range.front());
            });
        }
    };


    /**
     * @brief Set an equilibrium.
     *
     *
     * @param[out] rho_eq
     *      The equilibrium density @f$ \rho @f$.     
     * @param[in] function
     *      The function @f$ f @f$ used to compute the equilibrium.
     * @param[in] phi_max
     *      The maximal value of the electrical potential @f$ \phi @f$.
     * @param[in] tau
      *      The @f$ \tau @f$ parameter.
     *
     */
    void set_equilibrium(
            DFieldRTheta rho_eq,
            std::function<double(double const)> function,
            double const phi_max,
            double const tau)
    {
        IdxRangeRTheta grid = get_idx_range<GridR, GridTheta>(rho_eq);

        // Equilibrium:
        DFieldMemRTheta sigma_0(grid);
        DFieldMemRTheta phi_eq(grid);
        const double sig = 0.3;
        ddc::for_each(grid, [&](IdxRTheta const irp) {
            const CoordRTheta coord_rp(ddc::coordinate(irp));
            const CoordXY coord_xy(m_mapping(coord_rp));
            const double x = ddc::get<X>(coord_xy);
            const double y = ddc::get<Y>(coord_xy);
            sigma_0(irp) = sig;
            phi_eq(irp) = std::exp(-(x * x + y * y) / (2 * sig * sig));
        });
        find_equilibrium(
                get_field(sigma_0),
                get_field(phi_eq),
                get_field(rho_eq),
                function,
                phi_max,
                tau);
    };
};
