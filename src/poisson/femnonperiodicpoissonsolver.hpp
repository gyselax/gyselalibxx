// SPDX-License-Identifier: MIT

#pragma once

#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

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

    static inline Coordinate<QDimX> quad_point_from_coord(Coordinate<RDimX> const& coord)
    {
        return Coordinate<QDimX>(get<RDimX>(coord));
    }
    static inline Coordinate<RDimX> coord_from_quad_point(Coordinate<QDimX> const& coord)
    {
        return Coordinate<RDimX>(get<QDimX>(coord));
    }

private:
    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

    SplineEvaluator<NUBSplinesX> m_spline_x_nu_evaluator;

    ChargeDensityCalculator compute_rho;

private:
    // Spline degree in x direction
    static int constexpr m_degree = BSplinesX::degree();

    // Number of spline basis in x direction
    int const m_nbasis = discrete_space<BSplinesX>().nbasis();

    // Number of cells in x direction
    int const m_ncells = discrete_space<BSplinesX>().ncells();

    // Gauss points used for integration computation
    static int constexpr m_npts_gauss = m_degree + 1;
    Chunk<double, DiscreteDomain<QMeshX>> m_quad_coef;

    std::unique_ptr<Matrix> m_fem_matrix;

    void build_matrix();

    void solve_matrix_system(
            ChunkSpan<double, NUBSDomainX> phi_spline_coef,
            ChunkSpan<double, BSDomainX> rho_spline_coef) const;

public:
    FemNonPeriodicPoissonSolver(
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;
};
