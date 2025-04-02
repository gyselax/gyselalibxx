

# File advection\_field\_rtheta.hpp

[**File List**](files.md) **>** [**advection\_field**](dir_93bd9c5dcb1561501879bb6d3649fe64.md) **>** [**advection\_field\_rtheta.hpp**](advection__field__rtheta_8hpp.md)

[Go to the documentation of this file](advection__field__rtheta_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "iqnsolver.hpp"
#include "metric_tensor_evaluator.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polar_spline.hpp"
#include "polar_spline_evaluator.hpp"
#include "polarpoissonlikesolver.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"



template <class Mapping>
class AdvectionFieldFinder
{
public:
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

private:
    Mapping const& m_mapping;

    PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule> const
            m_polar_spline_evaluator;

    SplineRThetaEvaluatorNullBound_host const m_spline_evaluator;

    double const m_epsilon;

    static constexpr int n_overlap_cells = PolarBSplinesRTheta::continuity + 1;

public:
    explicit AdvectionFieldFinder(Mapping const& mapping, double const epsilon = 1e-12)
        : m_mapping(mapping)
        , m_polar_spline_evaluator(ddc::NullExtrapolationRule())
        , m_spline_evaluator {ddc::NullExtrapolationRule(), ddc::NullExtrapolationRule(), ddc::PeriodicExtrapolationRule<Theta>(), ddc::PeriodicExtrapolationRule<Theta>()}
        , m_epsilon(epsilon) {};

    // -------------------------------------------------------------------------------------------
    // COMPUTE ADVECTION FIELD IN XY:                                                            |
    // Advection field along the physical directions.                                            |
    // -------------------------------------------------------------------------------------------


    void operator()(
            host_t<DFieldRTheta> electrostatic_potential,
            host_t<DVectorFieldRTheta<X, Y>> advection_field_xy) const
    {
        IdxRangeRTheta const grid = get_idx_range(advection_field_xy);

        // Compute the spline representation of the electrostatic potential
        SplineRThetaBuilder_host const builder(grid);
        IdxRangeBSRTheta const idx_range_bsplinesRTheta = get_spline_idx_range(builder);
        host_t<Spline2DMem> electrostatic_potential_coef(idx_range_bsplinesRTheta);
        builder(get_field(electrostatic_potential_coef), get_const_field(electrostatic_potential));

        (*this)(get_field(electrostatic_potential_coef), advection_field_xy);
    }



    void operator()(
            host_t<Spline2D> electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<X, Y>> advection_field_xy) const
    {
        compute_advection_field_XY(
                m_spline_evaluator,
                electrostatic_potential_coef,
                advection_field_xy);
    }


    void operator()(
            host_t<PolarSplineMemRTheta>& electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<X, Y>> advection_field_xy) const
    {
        compute_advection_field_XY(
                m_polar_spline_evaluator,
                electrostatic_potential_coef,
                advection_field_xy);
    }


private:
    template <class SplineType, class Evaluator>
    void compute_advection_field_XY(
            Evaluator evaluator,
            SplineType& electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<X, Y>> advection_field_xy) const
    {
        static_assert(
                (std::is_same_v<
                         Evaluator,
                         SplineRThetaEvaluatorNullBound_host> && std::is_same_v<SplineType, host_t<Spline2D>>)
                || (std::is_same_v<
                            Evaluator,
                            PolarSplineEvaluator<
                                    PolarBSplinesRTheta,
                                    ddc::NullExtrapolationRule>> && std::is_same_v<SplineType, host_t<PolarSplineMemRTheta>>));

        IdxRangeRTheta const grid = get_idx_range(advection_field_xy);
        host_t<DVectorFieldMemRTheta<X, Y>> electric_field(grid);

        host_t<FieldMemRTheta<CoordRTheta>> coords(grid);
        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            coords(irtheta) = ddc::coordinate(irtheta);
        });

        // > computation of the phi derivatives
        host_t<DVectorFieldMemRTheta<R_cov, Theta_cov>> deriv_phi(grid);

        evaluator.deriv_dim_1(
                ddcHelper::get<R_cov>(deriv_phi),
                get_const_field(coords),
                get_const_field(electrostatic_potential_coef));
        evaluator.deriv_dim_2(
                ddcHelper::get<Theta_cov>(deriv_phi),
                get_const_field(coords),
                get_const_field(electrostatic_potential_coef));

        InverseJacobianMatrix<Mapping, CoordRTheta> inv_jacobian_matrix(m_mapping);

        // > computation of the electric field
        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            double const r = ddc::coordinate(ddc::select<GridR>(irtheta));
            double const th = ddc::coordinate(ddc::select<GridTheta>(irtheta));

            if (r > m_epsilon) {
                CoordRTheta const coord_rtheta(r, th);

                DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X, Y>> inv_J
                        = inv_jacobian_matrix(coord_rtheta);

                // Gradient of phi in the physical index range (Cartesian index range)
                // grad_{x,y} phi = J^{-T} grad_{r,theta} phi
                DVector<X, Y> grad_phi
                        = tensor_mul(index<'j', 'i'>(inv_J), index<'j'>(deriv_phi(irtheta)));

                // E = -grad phi
                ddcHelper::get<X>(electric_field)(irtheta) = -ddcHelper::get<X>(grad_phi);
                ddcHelper::get<Y>(electric_field)(irtheta) = -ddcHelper::get<Y>(grad_phi);

            } else {
                // Linearisation of the electric field
                double const th1 = M_PI / 4.;
                double const th2 = -M_PI / 4. + 2 * M_PI;

                // --- Value at r = 0:
                CoordRTheta const coord_1_0(0, th1);
                CoordRTheta const coord_2_0(0, th2);

                double const dr_x_1 = m_mapping.template jacobian_component<X, R_cov>(
                        coord_1_0); // dr_x (0, th1)
                double const dr_y_1 = m_mapping.template jacobian_component<Y, R_cov>(
                        coord_1_0); // dr_y (0, th1)

                double const dr_x_2 = m_mapping.template jacobian_component<X, R_cov>(
                        coord_2_0); // dr_x (0, th2)
                double const dr_y_2 = m_mapping.template jacobian_component<Y, R_cov>(
                        coord_2_0); // dr_y (0, th2)

                double deriv_r_phi_1 = evaluator.deriv_dim_1(
                        coord_1_0,
                        get_const_field(electrostatic_potential_coef));
                double deriv_r_phi_2 = evaluator.deriv_dim_1(
                        coord_2_0,
                        get_const_field(electrostatic_potential_coef));

                double const determinant = dr_x_1 * dr_y_2 - dr_x_2 * dr_y_1;

                DVector<X, Y> deriv_phi_0(
                        (dr_y_2 * deriv_r_phi_1 - dr_y_1 * deriv_r_phi_2) / determinant,
                        (-dr_x_2 * deriv_r_phi_1 + dr_x_1 * deriv_r_phi_2) / determinant);

                // E = -grad phi
                DVector<X, Y> electric_field_0 = -deriv_phi_0;

                // --- Value at r = m_epsilon:
                CoordRTheta const coord_rtheta_epsilon(m_epsilon, th);

                Tensor inv_J_eps = inv_jacobian_matrix(coord_rtheta_epsilon);

                DVector<R_cov, Theta_cov> deriv_phi_epsilon(
                        evaluator.deriv_dim_1(
                                coord_rtheta_epsilon,
                                get_const_field(electrostatic_potential_coef)),
                        evaluator.deriv_dim_2(
                                coord_rtheta_epsilon,
                                get_const_field(electrostatic_potential_coef)));

                // Gradient of phi in the physical domain (Cartesian domain)
                // (dx phi, dy phi) = J^{-T} (dr phi, dtheta phi)
                // E = -grad phi
                DVector<X, Y> electric_field_epsilon
                        = -tensor_mul(index<'j', 'i'>(inv_J_eps), index<'j'>(deriv_phi_epsilon));

                DVector<X, Y> E = electric_field_0 * (1 - r / m_epsilon)
                                  + electric_field_epsilon * r / m_epsilon;

                // --- Linearisation:
                ddcHelper::get<X>(electric_field)(irtheta) = ddcHelper::get<X>(E);
                ddcHelper::get<Y>(electric_field)(irtheta) = ddcHelper::get<Y>(E);
            }

            // > computation of the advection field
            ddcHelper::get<X>(advection_field_xy)(irtheta)
                    = ddcHelper::get<Y>(electric_field)(irtheta);
            ddcHelper::get<Y>(advection_field_xy)(irtheta)
                    = -ddcHelper::get<X>(electric_field)(irtheta);
        });
    }



public:
    // -------------------------------------------------------------------------------------------
    // COMPUTE ADVECTION FIELD IN RTheta:                                                        |
    // Advection field along the logical directions.                                             |
    // -------------------------------------------------------------------------------------------


    void operator()(
            host_t<DFieldRTheta> electrostatic_potential,
            host_t<DVectorFieldRTheta<R, Theta>> advection_field_rtheta,
            DVector<X, Y>& advection_field_xy_centre) const
    {
        IdxRangeRTheta const grid = get_idx_range(electrostatic_potential);

        // Compute the spline representation of the electrostatic potential
        SplineRThetaBuilder_host const builder(grid);
        IdxRangeBSRTheta const idx_range_bsplinesRTheta = get_spline_idx_range(builder);
        host_t<Spline2DMem> electrostatic_potential_coef(idx_range_bsplinesRTheta);
        builder(get_field(electrostatic_potential_coef), get_const_field(electrostatic_potential));

        (*this)(get_field(electrostatic_potential_coef),
                advection_field_rtheta,
                advection_field_xy_centre);
    }



    void operator()(
            host_t<Spline2D> electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<R, Theta>> advection_field_rtheta,
            DVector<X, Y>& advection_field_xy_centre) const
    {
        compute_advection_field_RTheta(
                m_spline_evaluator,
                electrostatic_potential_coef,
                advection_field_rtheta,
                advection_field_xy_centre);
    }


    void operator()(
            host_t<PolarSplineMemRTheta>& electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<R, Theta>> advection_field_rtheta,
            DVector<X, Y>& advection_field_xy_centre) const
    {
        compute_advection_field_RTheta(
                m_polar_spline_evaluator,
                electrostatic_potential_coef,
                advection_field_rtheta,
                advection_field_xy_centre);
    }



private:
    template <class SplineType, class Evaluator>
    void compute_advection_field_RTheta(
            Evaluator evaluator,
            SplineType& electrostatic_potential_coef,
            host_t<DVectorFieldRTheta<R, Theta>> advection_field_rtheta,
            DVector<X, Y>& advection_field_xy_centre) const
    {
        static_assert(
                (std::is_same_v<
                         Evaluator,
                         SplineRThetaEvaluatorNullBound_host> && std::is_same_v<SplineType, host_t<Spline2D>>)
                || (std::is_same_v<
                            Evaluator,
                            PolarSplineEvaluator<
                                    PolarBSplinesRTheta,
                                    ddc::NullExtrapolationRule>> && std::is_same_v<SplineType, host_t<PolarSplineMemRTheta>>));

        IdxRangeRTheta const grid_without_Opoint = get_idx_range(advection_field_rtheta);

        host_t<FieldMemRTheta<CoordRTheta>> coords(grid_without_Opoint);
        ddc::for_each(grid_without_Opoint, [&](IdxRTheta const irtheta) {
            coords(irtheta) = ddc::coordinate(irtheta);
        });

        // > computation of the phi derivatives
        host_t<DVectorFieldMemRTheta<R_cov, Theta_cov>> deriv_phi(grid_without_Opoint);

        evaluator.deriv_dim_1(
                ddcHelper::get<R_cov>(deriv_phi),
                get_const_field(coords),
                get_const_field(electrostatic_potential_coef));
        evaluator.deriv_dim_2(
                ddcHelper::get<Theta_cov>(deriv_phi),
                get_const_field(coords),
                get_const_field(electrostatic_potential_coef));

        MetricTensorEvaluator<Mapping, CoordRTheta> metric_tensor(m_mapping);

        // > computation of the advection field
        ddc::for_each(grid_without_Opoint, [&](IdxRTheta const irtheta) {
            CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

            DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<R, Theta>> inv_G
                    = metric_tensor.inverse(coord_rtheta);
            DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> J
                    = m_mapping.jacobian_matrix(coord_rtheta);
            double const jacobian = m_mapping.jacobian(coord_rtheta);

            // E = -grad phi
            DVector<R, Theta> electric_field
                    = -tensor_mul(index<'i', 'j'>(inv_G), index<'j'>(deriv_phi(irtheta)));

            // A (see README for the expression)
            ddcHelper::get<R>(advection_field_rtheta)(irtheta)
                    = (ddcHelper::get<X, R_cov>(J) * ddcHelper::get<X, Theta_cov>(J)
                       + ddcHelper::get<Y, R_cov>(J) * ddcHelper::get<Y, Theta_cov>(J))
                              * ddcHelper::get<R>(electric_field) / jacobian
                      + (ddcHelper::get<Y, Theta_cov>(J) * ddcHelper::get<Y, Theta_cov>(J)
                         + ddcHelper::get<X, Theta_cov>(J) * ddcHelper::get<X, Theta_cov>(J))
                                * ddcHelper::get<Theta>(electric_field) / jacobian;
            ddcHelper::get<Theta>(advection_field_rtheta)(irtheta)
                    = -(ddcHelper::get<X, R_cov>(J) * ddcHelper::get<X, R_cov>(J)
                        + ddcHelper::get<Y, R_cov>(J) * ddcHelper::get<Y, R_cov>(J))
                              * ddcHelper::get<R>(electric_field) / jacobian
                      - (ddcHelper::get<X, R_cov>(J) * ddcHelper::get<X, Theta_cov>(J)
                         + ddcHelper::get<Y, R_cov>(J) * ddcHelper::get<Y, Theta_cov>(J))
                                * ddcHelper::get<Theta>(electric_field) / jacobian;
        });

        // SPECIAL TREATMENT FOR THE O-POINT =====================================================
        // Linearisation of the electric field
        double const th1 = M_PI / 4.;
        double const th2 = -M_PI / 4. + 2 * M_PI;

        // --- Value at r = 0:
        CoordRTheta const coord_1_0(0, th1);
        CoordRTheta const coord_2_0(0, th2);

        double const dr_x_1
                = m_mapping.template jacobian_component<X, R_cov>(coord_1_0); // dr_x (0, th1)
        double const dr_y_1
                = m_mapping.template jacobian_component<Y, R_cov>(coord_1_0); // dr_y (0, th1)

        double const dr_x_2
                = m_mapping.template jacobian_component<X, R_cov>(coord_2_0); // dr_x (0, th2)
        double const dr_y_2
                = m_mapping.template jacobian_component<Y, R_cov>(coord_2_0); // dr_y (0, th2)

        double const deriv_r_phi_1
                = evaluator.deriv_dim_1(coord_1_0, get_const_field(electrostatic_potential_coef));
        double const deriv_r_phi_2
                = evaluator.deriv_dim_1(coord_2_0, get_const_field(electrostatic_potential_coef));

        double const determinant = dr_x_1 * dr_y_2 - dr_x_2 * dr_y_1;

        double const deriv_x_phi_0
                = (dr_y_2 * deriv_r_phi_1 - dr_y_1 * deriv_r_phi_2) / determinant;
        double const deriv_y_phi_0
                = (-dr_x_2 * deriv_r_phi_1 + dr_x_1 * deriv_r_phi_2) / determinant;

        advection_field_xy_centre = DVector<X, Y>(-deriv_y_phi_0, deriv_x_phi_0);
    }
};
```


