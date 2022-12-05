// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "chargedensitycalculator.hpp"
#include "ipoissonsolver.hpp"

using NUBSplinesX = NonUniformBSplines<RDimX, BSDegreeX>;
using NUBSDomainX = DiscreteDomain<NUBSplinesX>;
using UBSplinesX = UniformBSplines<RDimX, BSDegreeX>;

class FemNonPeriodicPoissonSolver : public IPoissonSolver
{
public:
    struct QDimX
    {
    };

private:
    using QMeshX = NonUniformPointSampling<QDimX>;

private:
    // Spline degree in x direction
    static int constexpr s_degree = BSplinesX::degree();

    // Gauss points used for integration computation
    static int constexpr s_npts_gauss = BSplinesX::degree() + 1;

private:
    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

    SplineEvaluator<NUBSplinesX> m_spline_x_nu_evaluator;

    ChargeDensityCalculator m_compute_rho;

    // Number of spline basis in x direction
    int m_nbasis;

    // Number of cells in x direction
    int m_ncells;

    Chunk<double, DiscreteDomain<QMeshX>> m_quad_coef;

    std::unique_ptr<Matrix> m_fem_matrix;

private:
    static Coordinate<QDimX> quad_point_from_coord(Coordinate<RDimX> const& coord)
    {
        return Coordinate<QDimX>(get<RDimX>(coord));
    }

    static Coordinate<RDimX> coord_from_quad_point(Coordinate<QDimX> const& coord)
    {
        return Coordinate<RDimX>(get<QDimX>(coord));
    }

public:
    FemNonPeriodicPoissonSolver(
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;

private:
    void build_matrix();

    void solve_matrix_system(
            ChunkSpan<double, NUBSDomainX> phi_spline_coef,
            ChunkSpan<double, BSDomainX> rho_spline_coef) const;
};
