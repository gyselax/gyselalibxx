// SPDX-License-Identifier: MIT

#pragma once

#include <sll/gauss_legendre_integration.hpp>
#include <sll/matrix.hpp>

#include <geometry.hpp>

#include "ichargedensitycalculator.hpp"
#include "ipoissonsolver.hpp"

using NUBSplinesX = ddc::NonUniformBSplines<RDimX, BSDegreeX>;
using NUBSDomainX = ddc::DiscreteDomain<NUBSplinesX>;
using UBSplinesX = ddc::UniformBSplines<RDimX, BSDegreeX>;

using NUBSplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        NUBSplinesX,
        IDimX,
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule,
        IDimX>;

/**
 * @brief An operator which solves the Poisson equation using Finite
 * Elements on a non-periodic domain.
 *
 * An operator which solves the Poisson equation:
 * @f$ - \frac{d^2 \phi}{dx^2} = \rho @f$
 * using Finite Elements on a non-periodic domain.
 * This solver uses spline basis elements and imposes Dirichlet boundary
 * conditions.
 */
class FemNonPeriodicPoissonSolver : public IPoissonSolver
{
public:
    /**
     * A tag to represent the dimension where the quadrature points
     * are defined.
     */
    struct QDimX
    {
    };

private:
    using QMeshX = ddc::NonUniformPointSampling<QDimX>;

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

    ddc::Chunk<double, ddc::DiscreteDomain<QMeshX>> m_quad_coef;

    std::unique_ptr<Matrix> m_fem_matrix;

private:
    static ddc::Coordinate<QDimX> quad_point_from_coord(ddc::Coordinate<RDimX> const& coord)
    {
        return ddc::Coordinate<QDimX>(ddc::get<RDimX>(coord));
    }

    static ddc::Coordinate<RDimX> coord_from_quad_point(ddc::Coordinate<QDimX> const& coord)
    {
        return ddc::Coordinate<RDimX>(ddc::get<QDimX>(coord));
    }

public:
    /**
     * Construct the FemNonPeriodicPoissonSolver operator.
     *
     * @param spline_x_builder A spline builder which calculates the coefficients of a spline representation.
     * @param spline_x_evaluator A spline evaluator which provides the value of a spline representation from its coefficients.
     * @param compute_rho The operator which calculates the charge density, the right hand side of the equation.
     */
    FemNonPeriodicPoissonSolver(
            SplineXBuilder_1d const& spline_x_builder,
            SplineXEvaluator_1d const& spline_x_evaluator,
            IChargeDensityCalculator const& compute_rho);

    /**
     * The operator which solves the equation using the method described by the class.
     *
     * @param[out] electrostatic_potential_device The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field_device The electric field, the derivative of the electrostatic potential.
     * @param[in] allfdistribu_device The distribution function.
     */
    void operator()(
            device_t<DSpanX> electrostatic_potential_device,
            device_t<DSpanX> electric_field_device,
            device_t<DViewSpXVx> allfdistribu_device) const override;

private:
    void build_matrix();

    void solve_matrix_system(
            ddc::ChunkSpan<double, NUBSDomainX> phi_spline_coef,
            ddc::ChunkSpan<double, BSDomainX> rho_spline_coef) const;
};
