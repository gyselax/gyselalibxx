// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "ichargedensitycalculator.hpp"
#include "ipoissonsolver.hpp"

/**
 * @brief An operator which solves the Poisson equation using Finite
 * Elements on a periodic domain.
 *
 * An operator which solves the Poisson equation:
 * @f$ - \frac{d^2 \phi}{dx^2} = \rho @f$
 * using Finite Elements on a periodic domain.
 * This solver uses spline basis elements and imposes Dirichlet boundary
 * conditions.
 */
class FemPeriodicPoissonSolver : public IPoissonSolver
{
public:
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
    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

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
     * Construct the FemPeriodicPoissonSolver operator.
     *
     * @param spline_x_builder A spline builder which calculates the coefficients of a spline representation.
     * @param spline_x_evaluator A spline evaluator which provides the value of a spline representation from its coefficients.
     * @param compute_rho The operator which calculates the charge density, the right hand side of the equation.
     */
    FemPeriodicPoissonSolver(
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
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
            ddc::ChunkSpan<double, BSDomainX> phi_spline_coef,
            ddc::ChunkSpan<double, BSDomainX> rho_spline_coef) const;
};
