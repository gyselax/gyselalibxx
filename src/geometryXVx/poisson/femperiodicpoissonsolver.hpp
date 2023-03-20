// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "chargedensitycalculator.hpp"
#include "ipoissonsolver.hpp"

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

    ChargeDensityCalculator m_compute_rho;

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
    FemPeriodicPoissonSolver(
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;

private:
    void build_matrix();

    void solve_matrix_system(
            ddc::ChunkSpan<double, BSDomainX> phi_spline_coef,
            ddc::ChunkSpan<double, BSDomainX> rho_spline_coef) const;
};
