// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include <sll/mapping/metric_tensor.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "directional_tag.hpp"
#include "geometry.hpp"
#include "iqnsolver.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"



/**
 * @brief Solve the Poisson-like equation and return the electric
 * field for the coupled Vlasov equation.
 *
 * The Vlasov-Poisson equations are given by
 *
 * - (1) @f$ \partial_t \rho - E_y \partial_x \rho + E_x \partial_y\rho = 0 @f$,
 *
 * - (2) @f$ - L \phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho  @f$,
 *
 * - (3) and @f$ E = -\nabla \phi  @f$.
 *
 * The functions are defined on a logical index range, and the mapping from the logical
 * index range to the physical index range is written @f$\mathcal{F}@f$.
 *
 * We here focus on equation (3). The @f$ \phi @f$ is already computed 
 * on B-splines with the given Poisson solver. Then in the AdvectionFieldRTheta::operator()
 * we compute the advection field (@f$A = E \wedge e_z@f$) thanks to (3) using the B-splines coefficients.
 * Depending on the given mapping, the computation at the center point is not
 * always well-defined so we linearize around the center point as explained
 * in Edoardo Zoni's article (https://doi.org/10.1016/j.jcp.2019.108889).
 * 
 * The advection field can be computed along the logical index range axis or the physical index range
 * axis. 
 * 
 * 1- In the first case, we compute the electric field thanks to (3) and 
 * - @f$ \nabla_{x,y} \phi(r, \theta) = (J^{-1})^{T} [\partial_r \phi, \partial_\theta \phi]^T @f$,
 * - @f$ E(r, \theta) = -\nabla_{x,y} \phi(r, \theta) @f$,
 * 
 * For @f$ r < \varepsilon @f$, @f$(J^{-1})^{T}@f$ is  ill-defined so we linearize 
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
 * Then the advection field along the physical index range axis 
 * is given by @f$A = E \wedge e_z@f$.
 * 
 * 
 * 2- In the second case, the advection field along the logical index range axis
 * is computed with 
 * - @f$ \nabla \phi = \sum_{i,j} \partial_{x_i} f g^{ij} \sqrt{g_{jj}} \hat{e}_j@f$, 
 * - with @f$g^{ij}@f$, the coefficients of the inverse metric tensor,
 * - @f$g_{jj}@f$, the coefficients of the metric tensor,
 * - @f$\hat{e}_j@f$, the normalized covariants vectors.
 * 
 * Then, we compute @f$ E = -\nabla \phi  @f$ and @f$A = E \wedge e_z@f$.
 * 
 *
 * The equation (1) is solved thanks to advection operator (IAdvectionRTheta).
 *
 *
 * @tparam Mapping
 *      A class describing a mapping from curvilinear coordinates to cartesian coordinates.
 *
 *
 * @see PolarSplineFEMPoissonLikeSolver
 *
 */
template <class Mapping>
class AdvectionFieldFinder
{
public:
    /**
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

private:
    Mapping const& m_mapping;

    PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule, Kokkos::HostSpace> const
            m_polar_spline_evaluator;

    SplineRThetaEvaluatorNullBound const m_spline_evaluator;

    double const m_epsilon;

    static constexpr int n_overlap_cells = PolarBSplinesRTheta::continuity + 1;

public:
    /**
     * @brief Instantiate a AdvectionFieldRTheta .
     *
     * @param[in] mapping
     *      The mapping @f$ \mathcal{F} @f$ from the logical index range to the physical index range.
     * @param[in] epsilon
     *      The parameter @f$ \varepsilon @f$ for the linearization of the
     *      electric field.
     */
    explicit AdvectionFieldFinder(Mapping const& mapping, double const epsilon = 1e-12)
        : m_mapping(mapping)
        , m_polar_spline_evaluator(ddc::NullExtrapolationRule())
        , m_spline_evaluator {ddc::NullExtrapolationRule(), ddc::NullExtrapolationRule(), ddc::PeriodicExtrapolationRule<Theta>(), ddc::PeriodicExtrapolationRule<Theta>()}
        , m_epsilon(epsilon) {};

    // -------------------------------------------------------------------------------------------
    // COMPUTE ADVECTION FIELD IN XY:                                                            |
    // Advection field along the physical directions.                                            |
    // -------------------------------------------------------------------------------------------


    /**
     * @brief Compute the advection field from a Field of @f$\phi@f$ values.
     *
     * @param[in] electrostatic_potential
     *      The values of the solution @f$\phi@f$ of the Poisson-like equation (2).
     * @param[out] advection_field_xy
     *      The advection field on the physical axis. 
     */
    void operator()(
            host_t<DFieldRTheta> electrostatic_potential,
            host_t<DVectorFieldRTheta<X, Y>> advection_field_xy) const
    {
        IdxRangeRTheta const grid = get_idx_range(advection_field_xy);

        // Compute the spline representation of the electrostatic potential
        SplineRThetaBuilder const builder(grid);
        IdxRangeBSRTheta const idx_range_bsplinesRTheta = get_spline_idx_range(builder);
        host_t<Spline2DMem> electrostatic_potential_coef(idx_range_bsplinesRTheta);
        builder(get_field(electrostatic_potential_coef), get_const_field(electrostatic_potential));

        (*this)(get_field(electrostatic_potential_coef), advection_field_xy);
    }



    /**
     * @brief Compute the advection field from a spline representation of @f$\phi@f$ solution.
     * The B-splines basis used is the cross-product of two 1D B-splines basis. 
     *
     * @param[in] electrostatic_potential_coef
     *      The spline representation of the solution @f$\phi@f$ of the Poisson-like equation (2).
     * @param[out] advection_field_xy
     *      The advection field on the physical axis. 
     */
    void operator()(
            host_t<Spline2D> electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<X, Y>> advection_field_xy) const
    {
        compute_advection_field_XY(
                m_spline_evaluator,
                electrostatic_potential_coef,
                advection_field_xy);
    }


    /**
     * @brief Compute the advection field from the Poisson-like equation solution.
     * The B-splines basis used is the polar B-splines (PolarSpline). 
     *
     * @param[in] electrostatic_potential_coef
     *      The polar spline representation of the solution @f$\phi@f$ of the Poisson-like equation (2).
     * @param[out] advection_field_xy
     *      The advection field on the physical axis. 
     */
    void operator()(
            host_t<SplinePolar>& electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<X, Y>> advection_field_xy) const
    {
        compute_advection_field_XY(
                m_polar_spline_evaluator,
                electrostatic_potential_coef,
                advection_field_xy);
    }


private:
    /**
     * @brief Compute the advection field along the physical axis.
     *
     * @param[in] evaluator 
     *      The spline evaluator used to evaluated electrostatic_potential_coef.
     * @param[in] electrostatic_potential_coef
     *      The spline representation of the solution @f$\phi@f$ of the Poisson-like equation (2).
     * @param[out] advection_field_xy
     *      The advection field on the physical axis. 
     */
    template <class SplineType, class Evaluator>
    void compute_advection_field_XY(
            Evaluator evaluator,
            SplineType& electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<X, Y>> advection_field_xy) const
    {
        static_assert(
                (std::is_same_v<
                         Evaluator,
                         SplineRThetaEvaluatorNullBound> && std::is_same_v<SplineType, host_t<Spline2D>>)
                || (std::is_same_v<
                            Evaluator,
                            PolarSplineEvaluator<
                                    PolarBSplinesRTheta,
                                    ddc::NullExtrapolationRule,
                                    Kokkos::HostSpace>> && std::is_same_v<SplineType, host_t<SplinePolar>>));

        IdxRangeRTheta const grid = get_idx_range(advection_field_xy);
        host_t<DVectorFieldMemRTheta<X, Y>> electric_field(grid);

        host_t<FieldMemRTheta<CoordRTheta>> coords(grid);
        ddc::for_each(grid, [&](IdxRTheta const irp) { coords(irp) = ddc::coordinate(irp); });

        // > computation of the phi derivatives
        host_t<DFieldMemRTheta> deriv_r_phi(grid);
        host_t<DFieldMemRTheta> deriv_p_phi(grid);

        evaluator.deriv_dim_1(
                get_field(deriv_r_phi),
                get_const_field(coords),
                get_const_field(electrostatic_potential_coef));
        evaluator.deriv_dim_2(
                get_field(deriv_p_phi),
                get_const_field(coords),
                get_const_field(electrostatic_potential_coef));

        InverseJacobianMatrix<Mapping, CoordRTheta> inv_jacobian_matrix(m_mapping);

        // > computation of the electric field
        ddc::for_each(grid, [&](IdxRTheta const irp) {
            double const r = ddc::coordinate(ddc::select<GridR>(irp));
            double const th = ddc::coordinate(ddc::select<GridTheta>(irp));

            if (r > m_epsilon) {
                CoordRTheta const coord_rp(r, th);

                Matrix_2x2 inv_J = inv_jacobian_matrix(coord_rp);

                // Gradiant of phi in the physical index range (Cartesian index range)
                double const deriv_x_phi
                        = deriv_r_phi(irp) * inv_J[0][0] + deriv_p_phi(irp) * inv_J[1][0];
                double const deriv_y_phi
                        = deriv_r_phi(irp) * inv_J[0][1] + deriv_p_phi(irp) * inv_J[1][1];

                // E = -grad phi
                ddcHelper::get<X>(electric_field)(irp) = -deriv_x_phi;
                ddcHelper::get<Y>(electric_field)(irp) = -deriv_y_phi;

            } else {
                // Linearisation of the electric field
                double const th1 = M_PI / 4.;
                double const th2 = -M_PI / 4. + 2 * M_PI;

                // --- Value at r = 0:
                CoordRTheta const coord_1_0(0, th1);
                CoordRTheta const coord_2_0(0, th2);

                double const dr_x_1 = m_mapping.jacobian_11(coord_1_0); // dr_x (0, th1)
                double const dr_y_1 = m_mapping.jacobian_21(coord_1_0); // dr_y (0, th1)

                double const dr_x_2 = m_mapping.jacobian_11(coord_2_0); // dr_x (0, th2)
                double const dr_y_2 = m_mapping.jacobian_21(coord_2_0); // dr_y (0, th2)

                double deriv_r_phi_1 = evaluator.deriv_dim_1(
                        coord_1_0,
                        get_const_field(electrostatic_potential_coef));
                double deriv_r_phi_2 = evaluator.deriv_dim_1(
                        coord_2_0,
                        get_const_field(electrostatic_potential_coef));

                double const determinant = dr_x_1 * dr_y_2 - dr_x_2 * dr_y_1;

                double const deriv_x_phi_0
                        = (dr_y_2 * deriv_r_phi_1 - dr_y_1 * deriv_r_phi_2) / determinant;
                double const deriv_y_phi_0
                        = (-dr_x_2 * deriv_r_phi_1 + dr_x_1 * deriv_r_phi_2) / determinant;

                // E = -grad phi
                double const electric_field_x_0 = -deriv_x_phi_0;
                double const electric_field_y_0 = -deriv_y_phi_0;



                // --- Value at r = m_epsilon:
                CoordRTheta const coord_rp_epsilon(m_epsilon, th);

                Matrix_2x2 inv_J_eps = inv_jacobian_matrix(coord_rp_epsilon);

                double const deriv_r_phi_epsilon = evaluator.deriv_dim_1(
                        coord_rp_epsilon,
                        get_const_field(electrostatic_potential_coef));
                double const deriv_p_phi_epsilon = evaluator.deriv_dim_2(
                        coord_rp_epsilon,
                        get_const_field(electrostatic_potential_coef));

                // Gradiant of phi in the physical index range (Cartesian index range)
                double const deriv_x_phi_epsilon = deriv_r_phi_epsilon * inv_J_eps[0][0]
                                                   + deriv_p_phi_epsilon * inv_J_eps[1][0];
                double const deriv_y_phi_epsilon = deriv_r_phi_epsilon * inv_J_eps[0][1]
                                                   + deriv_p_phi_epsilon * inv_J_eps[1][1];

                // E = -grad phi
                double const electric_field_x_epsilon = -deriv_x_phi_epsilon;
                double const electric_field_y_epsilon = -deriv_y_phi_epsilon;


                // --- Linearisation:
                ddcHelper::get<X>(electric_field)(irp) = electric_field_x_0 * (1 - r / m_epsilon)
                                                         + electric_field_x_epsilon * r / m_epsilon;
                ddcHelper::get<Y>(electric_field)(irp) = electric_field_y_0 * (1 - r / m_epsilon)
                                                         + electric_field_y_epsilon * r / m_epsilon;
            }

            // > computation of the advection field
            ddcHelper::get<X>(advection_field_xy)(irp) = -ddcHelper::get<Y>(electric_field)(irp);
            ddcHelper::get<Y>(advection_field_xy)(irp) = ddcHelper::get<X>(electric_field)(irp);
        });
    }



public:
    // -------------------------------------------------------------------------------------------
    // COMPUTE ADVECTION FIELD IN RTheta:                                                            |
    // Advection field along the logical directions.                                             |
    // -------------------------------------------------------------------------------------------


    /**
     * @brief Compute the advection field from a Field of @f$\phi@f$ values.
     *
     * @param[in] electrostatic_potential
     *      The values of the solution @f$\phi@f$ of the Poisson-like equation (2).
     * @param[out] advection_field_rp
     *      The advection field on the logical axis. 
     * @param[out] advection_field_xy_center
     *      The advection field on the physical axis at the O-point. 
     */
    void operator()(
            host_t<DFieldRTheta> electrostatic_potential,
            host_t<DVectorFieldRTheta<R, Theta>> advection_field_rp,
            CoordXY& advection_field_xy_center) const
    {
        IdxRangeRTheta const grid = get_idx_range(electrostatic_potential);

        // Compute the spline representation of the electrostatic potential
        SplineRThetaBuilder const builder(grid);
        IdxRangeBSRTheta const idx_range_bsplinesRTheta = get_spline_idx_range(builder);
        host_t<Spline2DMem> electrostatic_potential_coef(idx_range_bsplinesRTheta);
        builder(get_field(electrostatic_potential_coef), get_const_field(electrostatic_potential));

        (*this)(get_field(electrostatic_potential_coef),
                advection_field_rp,
                advection_field_xy_center);
    }



    /**
     * @brief Compute the advection field from a spline representation of @f$\phi@f$.
     * The B-splines basis used is the cross-product of two 1D B-splines basis. 
     *
     * @param[in] electrostatic_potential_coef
     *      The spline representation of the solution @f$\phi@f$ of the Poisson-like equation (2).
     * @param[out] advection_field_rp
     *      The advection field on the logical axis. 
     * @param[out] advection_field_xy_center
     *      The advection field on the physical axis at the O-point.  
     */
    void operator()(
            host_t<Spline2D> electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<R, Theta>> advection_field_rp,
            CoordXY& advection_field_xy_center) const
    {
        compute_advection_field_RTheta(
                m_spline_evaluator,
                electrostatic_potential_coef,
                advection_field_rp,
                advection_field_xy_center);
    }


    /**
     * @brief Compute the advection field from the Poisson-like equation.
     * The B-splines basis used is the polar B-splines (PolarSpline). 
     *
     * @param[in] electrostatic_potential_coef
     *      The polar spline representation of the solution @f$\phi@f$ of the Poisson-like equation (2).
     * @param[out] advection_field_rp
     *      The advection field on the logical axis. 
     * @param[out] advection_field_xy_center
     *      The advection field on the physical axis at the O-point. 
     */
    void operator()(
            host_t<SplinePolar>& electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<R, Theta>> advection_field_rp,
            CoordXY& advection_field_xy_center) const
    {
        compute_advection_field_RTheta(
                m_polar_spline_evaluator,
                electrostatic_potential_coef,
                advection_field_rp,
                advection_field_xy_center);
    }



private:
    /**
     * @brief Compute the advection field along the logical axis.
     *
     * @param[in] evaluator 
     *      The spline evaluator used to evaluated electrostatic_potential_coef.
     * @param[in] electrostatic_potential_coef
     *      The spline representation of the solution @f$\phi@f$ of the Poisson-like equation (2).
     * @param[out] advection_field_rp
     *      The advection field on the logical axis on an index range without O-point. 
     * @param[out] advection_field_xy_center
     *      The advection field on the physical axis at the O-point. 
     */
    template <class SplineType, class Evaluator>
    void compute_advection_field_RTheta(
            Evaluator evaluator,
            SplineType& electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<R, Theta>> advection_field_rp,
            CoordXY& advection_field_xy_center) const
    {
        static_assert(
                (std::is_same_v<
                         Evaluator,
                         SplineRThetaEvaluatorNullBound> && std::is_same_v<SplineType, host_t<Spline2D>>)
                || (std::is_same_v<
                            Evaluator,
                            PolarSplineEvaluator<
                                    PolarBSplinesRTheta,
                                    ddc::NullExtrapolationRule,
                                    Kokkos::HostSpace>> && std::is_same_v<SplineType, host_t<SplinePolar>>));

        IdxRangeRTheta const grid_without_Opoint = get_idx_range(advection_field_rp);

        host_t<FieldMemRTheta<CoordRTheta>> coords(grid_without_Opoint);
        ddc::for_each(grid_without_Opoint, [&](IdxRTheta const irp) {
            coords(irp) = ddc::coordinate(irp);
        });

        // > computation of the phi derivatives
        host_t<DFieldMemRTheta> deriv_r_phi(grid_without_Opoint);
        host_t<DFieldMemRTheta> deriv_p_phi(grid_without_Opoint);

        evaluator.deriv_dim_1(
                get_field(deriv_r_phi),
                get_const_field(coords),
                get_const_field(electrostatic_potential_coef));
        evaluator.deriv_dim_2(
                get_field(deriv_p_phi),
                get_const_field(coords),
                get_const_field(electrostatic_potential_coef));

        MetricTensor<Mapping, CoordRTheta> metric_tensor(m_mapping);

        // > computation of the advection field
        ddc::for_each(grid_without_Opoint, [&](IdxRTheta const irp) {
            CoordRTheta const coord_rp(ddc::coordinate(irp));

            Matrix_2x2 J; // Jacobian matrix
            m_mapping.jacobian_matrix(coord_rp, J);
            Matrix_2x2 inv_G; // Inverse metric tensor
            metric_tensor.inverse(inv_G, coord_rp);
            Matrix_2x2 G; // Metric tensor
            metric_tensor(G, coord_rp);

            // E = -grad phi
            double const electric_field_r
                    = (-deriv_r_phi(irp) * inv_G[0][0] - deriv_p_phi(irp) * inv_G[1][0])
                      * std::sqrt(G[0][0]);
            double const electric_field_p
                    = (-deriv_r_phi(irp) * inv_G[0][1] - deriv_p_phi(irp) * inv_G[1][1])
                      * std::sqrt(G[1][1]);

            // A = E \wedge e_z
            ddcHelper::get<R>(advection_field_rp)(irp) = -electric_field_p;
            ddcHelper::get<Theta>(advection_field_rp)(irp) = electric_field_r;
        });

        // SPECIAL TREATMENT FOR THE O-POINT =====================================================
        // Linearisation of the electric field
        double const th1 = M_PI / 4.;
        double const th2 = -M_PI / 4. + 2 * M_PI;

        // --- Value at r = 0:
        CoordRTheta const coord_1_0(0, th1);
        CoordRTheta const coord_2_0(0, th2);

        double const dr_x_1 = m_mapping.jacobian_11(coord_1_0); // dr_x (0, th1)
        double const dr_y_1 = m_mapping.jacobian_21(coord_1_0); // dr_y (0, th1)

        double const dr_x_2 = m_mapping.jacobian_11(coord_2_0); // dr_x (0, th2)
        double const dr_y_2 = m_mapping.jacobian_21(coord_2_0); // dr_y (0, th2)

        double const deriv_r_phi_1
                = evaluator.deriv_dim_1(coord_1_0, get_const_field(electrostatic_potential_coef));
        double const deriv_r_phi_2
                = evaluator.deriv_dim_1(coord_2_0, get_const_field(electrostatic_potential_coef));

        double const determinant = dr_x_1 * dr_y_2 - dr_x_2 * dr_y_1;

        double const deriv_x_phi_0
                = (dr_y_2 * deriv_r_phi_1 - dr_y_1 * deriv_r_phi_2) / determinant;
        double const deriv_y_phi_0
                = (-dr_x_2 * deriv_r_phi_1 + dr_x_1 * deriv_r_phi_2) / determinant;

        advection_field_xy_center = CoordXY(deriv_y_phi_0, -deriv_x_phi_0);
    }
};
