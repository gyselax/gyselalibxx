// SPDX-License-Identifier: MIT

#include <cassert>

#include <ddc/ddc.hpp>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/matrix.hpp>

#include <geometry.hpp>
#include <species_info.hpp>

#include "femperiodicpoissonsolver.hpp"

FemPeriodicPoissonSolver::FemPeriodicPoissonSolver(
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , compute_rho(spline_vx_builder, spline_vx_evaluator)
    , m_eval_pts_ptr(m_npts_gauss * m_ncells)
    , m_quad_coef_ptr(m_npts_gauss * m_ncells)
    , m_eval_pts(m_eval_pts_ptr.data(), m_npts_gauss * m_ncells)
    , m_quad_coef(m_quad_coef_ptr.data(), m_npts_gauss * m_ncells)
{
    assert(SplineXBuilder::bsplines_type::is_periodic());

    std::vector<double> knots_ptr(m_ncells + 1);
    DSpan1D knots(knots_ptr.data(), m_ncells + 1);

    for (size_t i(0); i < m_ncells + 1; ++i) {
        knots(i) = discrete_space<BSplinesX>().get_knot(i);
    }

    // Calculate the integration coefficients
    GaussLegendre const gl(m_npts_gauss);
    gl.compute_points_and_weights_on_mesh(m_eval_pts, m_quad_coef, knots);

    // Build the finite elements matrix
    build_matrix();
}


//===========================================================================
// Construct the finite element matrix with Bsplines as test functions
//===========================================================================
void FemPeriodicPoissonSolver::build_matrix()
{
    int constexpr n_lower_diags = m_degree + 1;
    int const matrix_size = m_nbasis + 1;
    bool const positive_definite_symmetric = false;

    // Matrix with block is used instead of periodic to contain the
    // Dirichlet boundary conditions
    m_fem_matrix = Matrix::make_new_block_with_banded_region(
            matrix_size,
            n_lower_diags,
            n_lower_diags,
            positive_definite_symmetric,
            n_lower_diags);

    // Fill the banded part of the matrix
    double derivs_ptr[m_degree + 1];
    DSpan1D derivs(derivs_ptr, m_degree + 1);
    for (int i = 0; i < m_eval_pts.size(); ++i) {
        auto jmin = discrete_space<BSplinesX>().eval_deriv(derivs, m_eval_pts[i]);
        for (int j = 0; j < m_degree + 1; ++j) {
            for (int k = 0; k < m_degree + 1; ++k) {
                int const j_idx = (j + jmin.uid()) % m_nbasis;
                int const k_idx = (k + jmin.uid()) % m_nbasis;
                double a_jk = m_fem_matrix->get_element(j_idx, k_idx);
                // Update element
                a_jk = a_jk + derivs(j) * derivs(k) * m_quad_coef[i];

                m_fem_matrix->set_element(j_idx, k_idx, a_jk);
            }
        }
    }

    // Impose the boundary conditions
    const DiscreteDomain bspline_full_domain = discrete_space<BSplinesX>().full_domain();
    const DiscreteDomain bspline_dom
            = bspline_full_domain.take_first(DiscreteVector<BSplinesX>(m_nbasis));

    Chunk<double, DiscreteDomain<BSplinesX>> int_vals(bspline_dom);
    discrete_space<BSplinesX>().integrals(int_vals);

    for (auto ix : bspline_dom) {
        const int i = ix.uid();
        m_fem_matrix->set_element(m_nbasis, i, int_vals(ix));
        m_fem_matrix->set_element(i, m_nbasis, int_vals(ix));
    }

    // Factorize the matrix ready to call solve
    m_fem_matrix->factorize();
}


//===========================================================================
//               Solve the Poisson equation
//---------------------------------------------------------------------------
void FemPeriodicPoissonSolver::solve_matrix_system(
        ChunkSpan<double, BSDomainX> phi_spline_coef,
        ChunkSpan<double, BSDomainX> rho_spline_coef) const
{
    int const nb_eval_pts = m_eval_pts.size();
    double values_ptr[m_degree + 1];
    DSpan1D values(values_ptr, m_degree + 1);

    int const rhs_size = m_nbasis + 1;
    DSpan1D phi_rhs(phi_spline_coef.data(), rhs_size);
    for (int i(0); i < rhs_size; ++i) {
        phi_rhs(i) = 0.0;
    }

    // Fill phi_rhs(i) with \int rho(x) b_i(x) dx
    // Rk: phi_rhs no longer contains spline coefficients, but is the
    //     RHS of the matrix equation
    for (int i = 0; i < nb_eval_pts; ++i) {
        auto jmin = discrete_space<BSplinesX>().eval_basis(values, m_eval_pts[i]);
        double const rho_val = m_spline_x_evaluator(m_eval_pts[i], rho_spline_coef);
        for (int j = 0; j < m_degree + 1; ++j) {
            int const j_idx = (jmin.uid() + j) % m_nbasis;
            phi_rhs(j_idx) = phi_rhs(j_idx) + rho_val * values(j) * m_quad_coef[i];
        }
    }

    // Solve the matrix equation to find the spline coefficients of phi
    m_fem_matrix->solve_inplace(phi_rhs);

    // Copy the first d coefficients into the last d coefficients
    // These coefficients refer to the same BSplines which cross the boundaries
    for (int i = 0; i < m_degree; i++) {
        phi_spline_coef(DiscreteElement<BSplinesX>(m_nbasis + i)) = phi_rhs(i);
    }
}



//===========================================================================
// Resolution of Poisson equation through the finite difference method (FDM)
//  -d^2Phi\dx^2 = rho
//
// THOMAS ALGORITHM FOR PERIODIC TRIDIAGONAL SYSTEMS
//
//       |a1 c1   ...     b1  | |Phi1  |  |rho1  |
//       |b2 a2 c2 ...    0   | |Phi2  |  |rho2  |
//       |0 b3 a3 c3 ...  0   | |Phi3  |  |rho3  |
//       |    ...             |*|...   |= |...   |
//       |     ...        cN-1| |PhiN-1|  |rhoN-1|
//       |cN    ...  bN   aN  | |PhiN  |  |rhoN  |
//
// with Lagrangian multipliers
//----------------------------------------------------------------------------
void FemPeriodicPoissonSolver::operator()(
        DSpanX electrostatic_potential,
        DSpanX electric_field,
        DViewSpXVx allfdistribu) const
{
    assert(electrostatic_potential.domain() == get_domain<IDimX>(allfdistribu));
    IDomainX dom_x = electrostatic_potential.domain();

    // Compute the RHS of the Poisson equation
    Chunk<double, IDomainX> rho(dom_x);
    compute_rho(rho, allfdistribu);

    //
    Chunk<double, BSDomainX> rho_spline_coef(m_spline_x_builder.spline_domain());
    m_spline_x_builder(rho_spline_coef, rho);
    Chunk<double, BSDomainX> phi_spline_coef(m_spline_x_builder.spline_domain());
    solve_matrix_system(phi_spline_coef, rho_spline_coef);

    //
    for_each(dom_x, [&](IndexX const ix) {
        electrostatic_potential(ix) = m_spline_x_evaluator(coordinate(ix), phi_spline_coef);
        electric_field(ix) = -m_spline_x_evaluator.deriv(coordinate(ix), phi_spline_coef);
    });
}
