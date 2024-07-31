#pragma once

#include <functional>

#include <ddc/ddc.hpp>

#include <geometry.hpp>

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
    IDomainRP const& m_grid;
    SplineRPBuilder const& m_builder;
    SplineRPEvaluatorNullBound const& m_evaluator;
    PolarSplineFEMPoissonLikeSolver const& m_poisson_solver;

public:
    /**
     * @brief Instantiate a VortexMergerEquilibria.
     *
     * @param[in] mapping
     *      The mapping function from the logical domain
     *      to the physical domain.
     * @param[in] grid
     *      The domain where the equilibrium is defined.
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
            IDomainRP const& grid,
            SplineRPBuilder const& builder,
            SplineRPEvaluatorNullBound const& evaluator,
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
            DSpanRP sigma,
            DSpanRP phi_eq,
            DSpanRP rho_eq,
            std::function<double(double const)> const& function,
            double const phi_max, // ToDo: ADD CASE WHERE RHO_MAX IS GIVEN.
            double const tau,
            int count_max = 25) const
    {
        DFieldRP phi_star(m_grid);
        DFieldRP ci(m_grid);

        auto dom_bsplinesRP = m_builder.spline_domain();
        Spline2D rho_coef(dom_bsplinesRP);

        FieldRP<CoordRP> coords(m_grid);
        ddc::for_each(m_grid, [&](IndexRP const irp) { coords(irp) = ddc::coordinate(irp); });

        double difference_sigma(0.);
        int count = 0;

        do {
            count += 1;
            // STEP 1: compute rho^i
            ddc::for_each(m_grid, [&](IndexRP const irp) {
                rho_eq(irp) = sigma(irp) * function(phi_eq(irp));
            });


            // STEP 2: compute phi_star^i with PDE solver
            m_builder(rho_coef.span_view(), rho_eq.span_cview());
            PoissonLikeRHSFunction poisson_rhs(rho_coef, m_evaluator);
            m_poisson_solver(poisson_rhs, coords.span_cview(), phi_star.span_view());

            // STEP 3: compute c^i
            // If phi_max is given:
            double norm_Linf_phi_star(0.);
            ddc::for_each(m_grid, [&](IndexRP const irp) {
                double const abs_phi_star = fabs(phi_star(irp));
                norm_Linf_phi_star
                        = norm_Linf_phi_star > abs_phi_star ? norm_Linf_phi_star : abs_phi_star;
            });

            ddc::for_each(m_grid, [&](IndexRP const irp) {
                ci(irp) = phi_max / norm_Linf_phi_star;
            });


            // STEP 4: update sigma and phi
            difference_sigma = 0.;
            ddc::for_each(m_grid, [&](IndexRP const irp) {
                double const abs_diff_sigma = fabs(sigma(irp) - ci(irp) * sigma(irp));
                difference_sigma
                        = difference_sigma > abs_diff_sigma ? difference_sigma : abs_diff_sigma;

                sigma(irp) = ci(irp) * sigma(irp);
                phi_eq(irp) = ci(irp) * phi_star(irp);
            });

        } while ((difference_sigma > tau) and (count < count_max));


        // STEP 1: compute rho^i
        ddc::for_each(m_grid, [&](IndexRP const irp) {
            rho_eq(irp) = sigma(irp) * function(phi_eq(irp));
        });

        // Unify at the center point:
        auto r_domain = ddc::get_domain<IDimR>(rho_eq);
        auto theta_domain = ddc::get_domain<IDimP>(rho_eq);
        if (std::fabs(ddc::coordinate(r_domain.front())) < 1e-15) {
            ddc::for_each(theta_domain, [&](const IndexP ip) {
                rho_eq(r_domain.front(), ip) = rho_eq(r_domain.front(), theta_domain.front());
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
            DSpanRP rho_eq,
            std::function<double(double const)> function,
            double const phi_max,
            double const tau)
    {
        IDomainRP grid = ddc::get_domain<IDimR, IDimP>(rho_eq);

        // Equilibrium:
        DFieldRP sigma_0(grid);
        DFieldRP phi_eq(grid);
        const double sig = 0.3;
        ddc::for_each(grid, [&](IndexRP const irp) {
            const CoordRP coord_rp(ddc::coordinate(irp));
            const CoordXY coord_xy(m_mapping(coord_rp));
            const double x = ddc::get<RDimX>(coord_xy);
            const double y = ddc::get<RDimY>(coord_xy);
            sigma_0(irp) = sig;
            phi_eq(irp) = std::exp(-(x * x + y * y) / (2 * sig * sig));
        });
        find_equilibrium(
                sigma_0.span_view(),
                phi_eq.span_view(),
                rho_eq.span_view(),
                function,
                phi_max,
                tau);
    };
};
