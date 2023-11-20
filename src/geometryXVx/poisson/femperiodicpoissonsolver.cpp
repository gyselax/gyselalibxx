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
    , m_compute_rho(spline_vx_builder, spline_vx_evaluator)
    , m_nbasis(ddc::discrete_space<BSplinesX>().nbasis())
    , m_ncells(ddc::discrete_space<BSplinesX>().ncells())
    , m_quad_coef(ddc::DiscreteDomain<QMeshX>(
              ddc::DiscreteElement<QMeshX>(0),
              ddc::DiscreteVector<QMeshX>(s_npts_gauss * m_ncells)))
{
    static_assert(SplineXBuilder::bsplines_type::is_periodic());

    BSDomainX const
            domain(ddc::DiscreteElement<BSplinesX>(0),
                   ddc::DiscreteVector<BSplinesX>(m_ncells + 1));
    ddc::Chunk<ddc::Coordinate<QDimX>, BSDomainX> knots(domain);

    for (ddc::DiscreteElement<BSplinesX> const i : domain) {
        knots(i) = quad_point_from_coord(ddc::discrete_space<BSplinesX>().get_knot(i.uid()));
    }

    // Calculate the integration coefficients
    GaussLegendre<QDimX> const gl(s_npts_gauss);
    std::vector<ddc::Coordinate<QDimX>> eval_pts_data(m_quad_coef.domain().size());
    ddc::ChunkSpan<ddc::Coordinate<QDimX>, ddc::DiscreteDomain<QMeshX>> const
            eval_pts(eval_pts_data.data(), m_quad_coef.domain());
    gl.compute_points_and_weights_on_mesh(eval_pts, m_quad_coef.span_view(), knots.span_cview());

    ddc::init_discrete_space<QMeshX>(eval_pts_data);

    // Build the finite elements matrix
    build_matrix();
}


//===========================================================================
// Construct the finite element matrix with Bsplines as test functions
//===========================================================================
void FemPeriodicPoissonSolver::build_matrix()
{
    int constexpr n_lower_diags = s_degree + 1;
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
    std::array<double, s_degree + 1> derivs_ptr;
    DSpan1D const derivs(derivs_ptr.data(), derivs_ptr.size());
    ddc::for_each(m_quad_coef.domain(), [&](ddc::DiscreteElement<QMeshX> const ix) {
        ddc::Coordinate<RDimX> const coord = coord_from_quad_point(ddc::coordinate(ix));
        ddc::DiscreteElement<BSplinesX> const jmin
                = ddc::discrete_space<BSplinesX>().eval_deriv(derivs, coord);
        for (int j = 0; j < s_degree + 1; ++j) {
            for (int k = 0; k < s_degree + 1; ++k) {
                int const j_idx = (j + jmin.uid()) % m_nbasis;
                int const k_idx = (k + jmin.uid()) % m_nbasis;
                double a_jk = m_fem_matrix->get_element(j_idx, k_idx);
                // Update element
                a_jk = a_jk + derivs(j) * derivs(k) * m_quad_coef(ix);

                m_fem_matrix->set_element(j_idx, k_idx, a_jk);
            }
        }
    });

    // Impose the boundary conditions
    BSDomainX const bspline_full_domain = ddc::discrete_space<BSplinesX>().full_domain();
    BSDomainX const bspline_dom
            = bspline_full_domain.take_first(ddc::DiscreteVector<BSplinesX>(m_nbasis));

    ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX>> int_vals(bspline_dom);
    ddc::discrete_space<BSplinesX>().integrals(int_vals);

    for (ddc::DiscreteElement<BSplinesX> const ix : bspline_dom) {
        int const i = ix.uid();
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
        ddc::ChunkSpan<double, BSDomainX> const phi_spline_coef,
        ddc::ChunkSpan<double, BSDomainX> const rho_spline_coef) const
{
    std::array<double, s_degree + 1> values_ptr;
    DSpan1D const values(values_ptr.data(), values_ptr.size());

    int const rhs_size = m_nbasis + 1;
    DSpan1D const phi_rhs(phi_spline_coef.data_handle(), rhs_size);
    for (int i(0); i < rhs_size; ++i) {
        phi_rhs(i) = 0.0;
    }

    // Fill phi_rhs(i) with \int rho(x) b_i(x) dx
    // Rk: phi_rhs no longer contains spline coefficients, but is the
    //     RHS of the matrix equation
    ddc::for_each(m_quad_coef.domain(), [&](ddc::DiscreteElement<QMeshX> const ix) {
        ddc::Coordinate<RDimX> const coord = coord_from_quad_point(ddc::coordinate(ix));
        ddc::DiscreteElement<BSplinesX> const jmin
                = ddc::discrete_space<BSplinesX>().eval_basis(values, coord);
        double const rho_val = m_spline_x_evaluator(coord, rho_spline_coef);
        for (int j = 0; j < s_degree + 1; ++j) {
            int const j_idx = (jmin.uid() + j) % m_nbasis;
            phi_rhs(j_idx) = phi_rhs(j_idx) + rho_val * values(j) * m_quad_coef(ix);
        }
    });

    // Solve the matrix equation to find the spline coefficients of phi
    m_fem_matrix->solve_inplace(phi_rhs);

    // Copy the first d coefficients into the last d coefficients
    // These coefficients refer to the same BSplines which cross the boundaries
    for (int i = 0; i < s_degree; i++) {
        phi_spline_coef(ddc::DiscreteElement<BSplinesX>(m_nbasis + i)) = phi_rhs(i);
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
        DSpanX const electrostatic_potential,
        DSpanX const electric_field,
        DViewSpXVx const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("PoissonSolver");
    assert(electrostatic_potential.domain() == ddc::get_domain<IDimX>(allfdistribu));
    IDomainX const dom_x = electrostatic_potential.domain();

    // Compute the RHS of the Poisson equation
    ddc::Chunk<double, IDomainX> rho(dom_x);
    m_compute_rho(rho, allfdistribu);

    //
    ddc::Chunk<double, BSDomainX> rho_spline_coef(m_spline_x_builder.spline_domain());
    m_spline_x_builder(rho_spline_coef, rho);
    ddc::Chunk<double, BSDomainX> phi_spline_coef(m_spline_x_builder.spline_domain());
    solve_matrix_system(phi_spline_coef, rho_spline_coef);

    //
    ddc::for_each(dom_x, [&](IndexX const ix) {
        electrostatic_potential(ix) = m_spline_x_evaluator(ddc::coordinate(ix), phi_spline_coef);
        electric_field(ix) = -m_spline_x_evaluator.deriv(ddc::coordinate(ix), phi_spline_coef);
    });
    Kokkos::Profiling::popRegion();
}