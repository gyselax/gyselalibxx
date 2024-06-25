// SPDX-License-Identifier: MIT

#pragma once

#include <sll/gauss_legendre_integration.hpp>
#include <sll/matrix.hpp>

#include <geometry.hpp>

#include "ichargedensitycalculator.hpp"
#include "iqnsolver.hpp"

struct HiddenNUBSplinesX : ddc::NonUniformBSplines<RDimX, BSDegreeX>
{
};
using NUBSplinesX
        = std::conditional_t<ddc::is_uniform_bsplines_v<BSplinesX>, HiddenNUBSplinesX, BSplinesX>;
using NUBSDomainX = ddc::DiscreteDomain<NUBSplinesX>;

using NUBSplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        NUBSplinesX,
        IDimX,
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule,
        IDimX>;

/**
 * @brief An operator which solves the Quasi-Neutrality equation using Finite
 * Elements on a non-periodic domain.
 *
 * An operator which solves the Quasi-Neutrality equation:
 * @f$ - \frac{d^2 \phi}{dx^2} = \rho @f$
 * using Finite Elements on a non-periodic domain.
 * This solver uses spline basis elements and imposes Dirichlet boundary
 * conditions.
 */
class FemNonPeriodicQNSolver : public IQNSolver
{
public:
    /// The discrete dimension of the quadrature points.
    struct QMeshX : ddc::NonUniformPointSampling<RDimX>
    {
    };

private:
    using DQFieldX = device_t<ddc::Chunk<double, ddc::DiscreteDomain<QMeshX>>>;
    using DNUBSSpanX = device_t<ddc::ChunkSpan<double, NUBSDomainX>>;

private:
    // Spline degree in x direction
    static int constexpr s_degree = BSplinesX::degree();

    // Gauss points used for integration computation
    static int constexpr s_npts_gauss = BSplinesX::degree() + 1;

private:
    SplineXBuilder_1d const& m_spline_x_builder;

    SplineXEvaluator_1d m_spline_x_evaluator;

    NUBSplineXEvaluator_1d m_spline_x_nu_evaluator;

    IChargeDensityCalculator const& m_compute_rho;

    // Number of spline basis in x direction
    int m_nbasis;

    // Number of cells in x direction
    int m_ncells;

    DQFieldX m_quad_coef_alloc;

    std::unique_ptr<Matrix> m_fem_matrix;

private:
    static KOKKOS_FUNCTION ddc::Coordinate<RDimX> quad_point_from_coord(
            ddc::Coordinate<RDimX> const& coord)
    {
        return ddc::Coordinate<RDimX>(ddc::get<RDimX>(coord));
    }

    static KOKKOS_FUNCTION ddc::Coordinate<RDimX> coord_from_quad_point(
            ddc::Coordinate<RDimX> const& coord)
    {
        return ddc::Coordinate<RDimX>(ddc::get<RDimX>(coord));
    }

public:
    /**
     * Construct the FemNonPeriodicQNSolver operator.
     *
     * @param spline_x_builder A spline builder which calculates the coefficients of a spline representation.
     * @param spline_x_evaluator A spline evaluator which provides the value of a spline representation from its coefficients.
     * @param compute_rho The operator which calculates the charge density, the right hand side of the equation.
     */
    FemNonPeriodicQNSolver(
            SplineXBuilder_1d const& spline_x_builder,
            SplineXEvaluator_1d const& spline_x_evaluator,
            IChargeDensityCalculator const& compute_rho);

    /**
     * The operator which solves the equation using the method described by the class.
     *
     * @param[out] electrostatic_potential The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field The electric field, the derivative of the electrostatic potential.
     * @param[in] allfdistribu The distribution function.
     */
    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;

private:
    void build_matrix();

public:
    /**
     * [SHOULD BE PRIVATE (Kokkos limitation)]
     *
     * @param[out] phi_spline_coef
     * @param[in] rho_spline_coef
     */
    void solve_matrix_system(DNUBSSpanX phi_spline_coef, DBSViewX rho_spline_coef) const;
};
