#pragma once


#include <ddc/ddc.hpp>

#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>

#include <directional_tag.hpp>
#include <geometry.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

#include "ipoissonsolver.hpp"
#include "poisson_rhs_function.hpp"
#include "polarpoissonsolver.hpp"



/**
 * @brief Solve the Poisson equation and return the electric
 * field for the coupled Vlasov equation.
 *
 * The Vlasov-Poisson equations are given by
 *
 * - (1) @f$ \partial_t \rho - E_y \partial_x \rho + E_x \partial_y\rho = 0 @f$,
 *
 * - (2) @f$ - \Delta \phi = \rho  @f$,
 *
 * - (3) and @f$ E = -\nabla \phi  @f$.
 *
 * The functions are defined on a logical domain, and the mapping from the logical
 * domain to the physical domain is written @f$\mathcal{F}@f$.
 *
 * We here focus on equations (2) and (3). We first compute @f$ \phi @f$
 * on B-splines with the given Poisson solver. Then in the VlasovPoissonSolver::operator()
 * we compute the electric field thanks to (3) using the B-splines coefficients.
 * Depending on the given mapping, the computation at the center point is not
 * always well-defined so we linearize around the center point as explained
 * in Edoardo Zoni's article (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * With (3), we compute the electric field
 *
 * - @f$ E(r, \theta) = -\nabla_{x,y} \phi(r, \theta) @f$,
 *
 * with @f$\nabla_{x,y} \phi(r, \theta) @f$ from
 * @f$ \nabla_{x,y} \phi(r, \theta)
 * = (J_{\mathcal{F}}^{-1})^T \nabla_{r, \theta} \phi(r, \theta)  @f$,
 *
 * and
 * @f$ \nabla \phi = \sum_{ij} \partial_{x_i} \phi g^{i,j} e_j @f$
 * with @f$ g^{i,j} @f$ the coefficients of the inverse metric tensor
 * and @f$ e_j @f$ the unnormalized local covariant base:
 *
 *
 * - @f$ \nabla \phi(r,\theta) \cdot e_r = \partial_r \phi(r,\theta) g^{1,1}  + \partial_\theta \phi(r,\theta)g^{2,1} @f$,
 *
 * - @f$ \nabla \phi(r,\theta) \cdot e_\theta = \partial_r \phi(r,\theta) g^{1,2}  + \partial_\theta \phi(r,\theta)g^{2,2} @f$.
 *
 *
 *
 * For @f$ r < \varepsilon @f$,
 * @f$ E(r, \theta) = \left( 1 - \frac{r}{\varepsilon} \right)  E(0, \theta)
 * + \frac{r}{\varepsilon} E(\varepsilon, \theta) @f$,
 *
 * with @f$ E(0, \theta) @f$ computed thanks to
 *
 *  - @f$ \partial_r \phi (0, \theta_1) = \left[\partial_r x  \partial_x \phi
 * + \partial_r y  \partial_y \phi \right](0, \theta_1) @f$,
 *
 *  - @f$ \partial_r \phi (0, \theta_2) = \left[\partial_r x  \partial_x \phi
 * + \partial_r y  \partial_y \phi \right] (0, \theta_2) @f$,
 *
 * where @f$ \theta_1 @f$ and @f$ \theta_2 @f$  correspond to
 * linearly independent directions.
 *
 *
 * The equation (1) is computed thanks to advection operator.
 *
 *
 * @tparam Mapping A Curvilinear2DToCartesian class.
 *
 *
 * @see PolarSplineFEMPoissonSolver
 *
 */
template <class Mapping>
class VlasovPoissonSolver : public IPoissonSolver
{
public:
    /**
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

private:
    Mapping const& m_mapping;

    SplineRPBuilder const& m_rhs_builder;
    SplineRPEvaluator const& m_rhs_evaluator;

    PolarSplineEvaluator<PolarBSplinesRP> const m_polar_spline_evaluator;

    PolarSplineFEMPoissonSolver const& m_poisson_solver;

    double const m_epsilon;


    static constexpr int n_overlap_cells = PolarBSplinesRP::continuity + 1;

public:
    /**
	 * @brief Instantiate a VlasovPoissonSolver.
	 *
	 * @param[in] mapping
	 *      The mapping @f$ \mathcal{F} @f$ from the logical domain to the physical domain.
	 * @param[in] builder
	 *      The splines builder for the rhs.
	 * @param[in] evaluator
	 *      The splines evaluator for the rhs.
	 * @param[in] poisson_solver
	 *      The Poisson solver used to solve (2).
	 * @param[in] epsilon
	 *      The parameter @f$ \varepsilon @f$ for the linearization of the
	 *      electric field.
	 */
    VlasovPoissonSolver(
            Mapping const& mapping,
            SplineRPBuilder const& builder,
            SplineRPEvaluator const& evaluator,
            PolarSplineFEMPoissonSolver const& poisson_solver,
            double const epsilon = 1e-12)
        : m_mapping(mapping)
        , m_rhs_builder(builder)
        , m_rhs_evaluator(evaluator)
        , m_polar_spline_evaluator(g_polar_null_boundary_2d<PolarBSplinesRP>)
        , m_poisson_solver(poisson_solver)
        , m_epsilon(epsilon) {};

    ~VlasovPoissonSolver() {};


    /**
	 * @brief Compute the electric field from the Poisson equation.
	 *
	 * @param[out] electrostatic_potential
	 *      The solution @f$\phi@f$ of the Poisson equation (2).
	 * @param[out] electric_field
     *      The electric field @f$E = -\nabla \phi@f$.
     * @param[in] allfdistribu
     *      The rhs @f$\rho@f$ of the Poisson equation (2).
	 */
    void operator()(
            DSpanRP electrostatic_potential,
            VectorFieldSpan<double, IDomainRP, NDTag<RDimX, RDimY>> electric_field,
            DViewRP allfdistribu) const
    {
        auto const grid = ddc::get_domain<IDimR, IDimP>(allfdistribu);

        FieldRP<CoordRP> coords(grid);
        ddc::for_each(grid, [&](IndexRP const irp) { coords(irp) = ddc::coordinate(irp); });

        auto const dom_bsplinesRP = m_rhs_builder.spline_domain();

        Spline2D allfdistribu_coef(dom_bsplinesRP);

        BSDomainR radial_bsplines(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                ddc::DiscreteVector<BSplinesR> {n_overlap_cells}));
        BSDomainP polar_domain(ddc::discrete_space<BSplinesP>().full_domain());
        SplinePolar electrostatic_potential_coef(
                PolarBSplinesRP::singular_domain(),
                BSDomainRP(radial_bsplines, polar_domain));


        m_rhs_builder(allfdistribu_coef, allfdistribu);
        PoissonRHSFunction const charge_density_coord(allfdistribu_coef, m_rhs_evaluator);

        m_poisson_solver(charge_density_coord, electrostatic_potential_coef);


        m_polar_spline_evaluator(
                electrostatic_potential.span_view(),
                coords.span_cview(),
                electrostatic_potential_coef);


        // Computation of the advection_field:
        DFieldRP deriv_r_phi(grid);
        DFieldRP deriv_p_phi(grid);


        m_polar_spline_evaluator.deriv_dim_1(
                deriv_r_phi.span_view(),
                coords.span_cview(),
                electrostatic_potential_coef);
        m_polar_spline_evaluator.deriv_dim_2(
                deriv_p_phi.span_view(),
                coords.span_cview(),
                electrostatic_potential_coef);

        ddc::for_each(grid, [&](IndexRP const irp) {
            double const r = ddc::coordinate(ddc::select<IDimR>(irp));
            double const th = ddc::coordinate(ddc::select<IDimP>(irp));

            if (r > m_epsilon) {
                CoordRP const coord_rp(r, th);

                Matrix_2x2 J; // Jacobian matrix
                m_mapping.jacobian_matrix(coord_rp, J);
                Matrix_2x2 inv_G; // Inverse metric tensor
                m_mapping.inverse_metric_tensor(coord_rp, inv_G);

                // E = -grad phi
                double const electric_field_r
                        = -deriv_r_phi(irp) * inv_G[0][0] - deriv_p_phi(irp) * inv_G[1][0];
                double const electric_field_th
                        = -deriv_r_phi(irp) * inv_G[0][1] - deriv_p_phi(irp) * inv_G[1][1];

                // compute the electric field in the physical domain (Cartesian domain)
                ddcHelper::get<RDimX>(electric_field)(irp)
                        = electric_field_r * J[0][0] + electric_field_th * J[0][1];
                ddcHelper::get<RDimY>(electric_field)(irp)
                        = electric_field_r * J[1][0] + electric_field_th * J[1][1];

            } else {
                // Linearisation of the electric field
                double const th1 = M_PI / 4.;
                double const th2 = -M_PI / 4. + 2 * M_PI;

                // --- Value at r = 0:
                CoordRP const coord_1_0(0, th1);
                CoordRP const coord_2_0(0, th2);

                double const dr_x_1 = m_mapping.jacobian_11(coord_1_0); // dr_x (0, th1)
                double const dr_y_1 = m_mapping.jacobian_21(coord_1_0); // dr_y (0, th1)

                double const dr_x_2 = m_mapping.jacobian_11(coord_2_0); // dr_x (0, th2)
                double const dr_y_2 = m_mapping.jacobian_21(coord_2_0); // dr_y (0, th2)

                double deriv_r_phi_1
                        = m_polar_spline_evaluator
                                  .deriv_dim_1(coord_1_0, electrostatic_potential_coef);
                double deriv_r_phi_2
                        = m_polar_spline_evaluator
                                  .deriv_dim_1(coord_2_0, electrostatic_potential_coef);

                double const deriv_y_phi_0 = 1 / (dr_y_2 - dr_y_1 * dr_x_2 / dr_x_1)
                                             * (deriv_r_phi_2 - deriv_r_phi_1 * dr_x_2 / dr_x_1);

                double const deriv_x_phi_0 = 1 / dr_x_1 * (deriv_r_phi_1 - deriv_y_phi_0 * dr_y_1);

                // E = -grad phi
                double const electric_field_x_0 = -deriv_x_phi_0;
                double const electric_field_y_0 = -deriv_y_phi_0;



                // --- Value at r = m_epsilon:
                CoordRP const coord_rp_epsilon(m_epsilon, th);
                Matrix_2x2 J_eps; // Jacobian matrix
                m_mapping.jacobian_matrix(coord_rp_epsilon, J_eps);
                Matrix_2x2 inv_G_eps; // Inverse metric tensor
                m_mapping.inverse_metric_tensor(coord_rp_epsilon, inv_G_eps);

                double deriv_r_phi_epsilon
                        = m_polar_spline_evaluator
                                  .deriv_dim_1(coord_rp_epsilon, electrostatic_potential_coef);
                double deriv_p_phi_epsilon
                        = m_polar_spline_evaluator
                                  .deriv_dim_2(coord_rp_epsilon, electrostatic_potential_coef);


                // E = -grad phi
                double const electric_field_r = -deriv_r_phi_epsilon * inv_G_eps[0][0]
                                                - deriv_p_phi_epsilon * inv_G_eps[1][0];
                double const electric_field_th = -deriv_r_phi_epsilon * inv_G_eps[0][1]
                                                 - deriv_p_phi_epsilon * inv_G_eps[1][1];

                // compute the electric field in the physical domain (cartesian domain)
                double const electric_field_x_epsilon
                        = electric_field_r * J_eps[0][0] + electric_field_th * J_eps[0][1];
                double const electric_field_y_epsilon
                        = electric_field_r * J_eps[1][0] + electric_field_th * J_eps[1][1];


                // --- Linearisation:
                ddcHelper::get<RDimX>(electric_field)(irp)
                        = electric_field_x_0 * (1 - r / m_epsilon)
                          + electric_field_x_epsilon * r / m_epsilon;
                ddcHelper::get<RDimY>(electric_field)(irp)
                        = electric_field_y_0 * (1 - r / m_epsilon)
                          + electric_field_y_epsilon * r / m_epsilon;
            }
        });
    }
};
