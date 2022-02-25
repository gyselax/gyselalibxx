// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "chargedensitycalculator.hpp"
#include "ipoissonsolver.hpp"

class SpeciesInformation;

class FemPeriodicPoissonSolver : public IPoissonSolver
{
    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

    ChargeDensityCalculator compute_rho;

private:
    // Spline degree in x direction
    static int constexpr m_degree = BSplinesX::degree();

    // Number of spline basis in x direction
    int const m_nbasis = discretization<BSplinesX>().nbasis();

    // Number of cells in x direction
    int const m_ncells = discretization<BSplinesX>().ncells();

    // Gauss points used for integration computation
    static int constexpr m_npts_gauss = m_degree + 1;
    std::vector<double> m_eval_pts_ptr;
    std::vector<double> m_quad_coef_ptr;
    DSpan1D m_eval_pts;
    DSpan1D m_quad_coef;

    std::unique_ptr<Matrix> m_fem_matrix;

    void build_matrix();

    void solve_matrix_system(
            ChunkSpan<double, BSDomainX> phi_spline_coef,
            ChunkSpan<double, BSDomainX> rho_spline_coef) const;

public:
    FemPeriodicPoissonSolver(
            SpeciesInformation const& species_info,
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;
};
