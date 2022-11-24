// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/matrix.hpp>
#include <sll/null_boundary_value.hpp>

#include <geometry.hpp>
#include <species_info.hpp>

#include "femnonperiodicpoissonsolver.hpp"

//===========================================================================
// Create a non-uniform spline evaluator from the provided evaluator
//===========================================================================
template <class BSplines>
static SplineEvaluator<NUBSplinesX> jit_build_nubsplinesx(
        SplineEvaluator<BSplines> const& spline_x_evaluator)
{
    static_assert(std::is_same_v<BSplines, UBSplinesX> || std::is_same_v<BSplines, NUBSplinesX>);
    if constexpr (std::is_same_v<BSplines, UBSplinesX>) {
        int ncells = discrete_space<UBSplinesX>().ncells();
        std::vector<CoordX> knots(ncells + 1);

        for (int i(0); i < ncells + 1; ++i) {
            knots[i] = CoordX(discrete_space<UBSplinesX>().get_knot(i));
        }
        init_discrete_space<NUBSplinesX>(knots);
        // Boundary values are never evaluated
        return SplineEvaluator<
                NUBSplinesX>(g_null_boundary<NUBSplinesX>, g_null_boundary<NUBSplinesX>);
    } else {
        return spline_x_evaluator;
    }
}

//===========================================================================
// Build an instance of FemNonPeriodicPoissonSolver
//===========================================================================
FemNonPeriodicPoissonSolver::FemNonPeriodicPoissonSolver(
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_spline_x_nu_evaluator(jit_build_nubsplinesx(spline_x_evaluator))
    , compute_rho(spline_vx_builder, spline_vx_evaluator)
    , m_quad_coef(DiscreteDomain<QMeshX>(
              DiscreteElement<QMeshX>(0),
              DiscreteVector<QMeshX>((m_degree + 1) * m_ncells)))
{
    static_assert(!SplineXBuilder::bsplines_type::is_periodic());
    typename BSplinesX::discrete_domain_type
            domain(typename BSplinesX::discrete_element_type {0},
                   typename BSplinesX::discrete_vector_type {m_ncells + 1});
    Chunk<Coordinate<QDimX>, typename BSplinesX::discrete_domain_type> knots(domain);

    for (auto i : domain) {
        knots(i) = quad_point_from_coord(discrete_space<NUBSplinesX>().get_knot(i.uid()));
    }

    // Calculate the integration coefficients
    GaussLegendre<QDimX> const gl(m_npts_gauss);
    std::vector<Coordinate<QDimX>> eval_pts_data(m_quad_coef.domain().size());
    ChunkSpan<Coordinate<QDimX>, DiscreteDomain<QMeshX>>
            eval_pts(eval_pts_data.data(), m_quad_coef.domain());
    gl.compute_points_and_weights_on_mesh(eval_pts, m_quad_coef.span_view(), knots.span_cview());

    init_discrete_space<QMeshX>(eval_pts_data);

    // Build the finite elements matrix
    build_matrix();
}

//===========================================================================
// Construct the finite element matrix with Bsplines as test functions
//===========================================================================
void FemNonPeriodicPoissonSolver::build_matrix()
{
    // Matrix contains all elements of the basis except the first
    // and last in order to impose Dirichlet conditions
    int const matrix_size = m_nbasis - 2;
    int constexpr n_lower_diags = m_degree;
    bool const positive_definite_symmetric = false;

    // Matrix with block is used instead of periodic to contain the
    // Dirichlet boundary conditions
    m_fem_matrix = Matrix::
            make_new_banded(matrix_size, n_lower_diags, n_lower_diags, positive_definite_symmetric);

    // Fill the banded part of the matrix
    double derivs_ptr[m_degree + 1];
    DSpan1D derivs(derivs_ptr, m_degree + 1);
    for_each(m_quad_coef.domain(), [&](DiscreteElement<QMeshX> const ix) {
        Coordinate<RDimX> coord = coord_from_quad_point(coordinate(ix));
        auto jmin = discrete_space<NUBSplinesX>().eval_deriv(derivs, coord);
        for (int j = 0; j < m_degree + 1; ++j) {
            for (int k = 0; k < m_degree + 1; ++k) {
                int const j_idx = (j + jmin.uid()) % m_nbasis - 1;
                int const k_idx = (k + jmin.uid()) % m_nbasis - 1;

                if (j_idx != -1 && j_idx != matrix_size && k_idx != -1 && k_idx != matrix_size) {
                    double a_jk = m_fem_matrix->get_element(j_idx, k_idx);
                    // Update element
                    a_jk += derivs(j) * derivs(k) * m_quad_coef(ix);

                    m_fem_matrix->set_element(j_idx, k_idx, a_jk);
                }
            }
        }
    });

    // Factorize the matrix ready to call solve
    m_fem_matrix->factorize();
}


//===========================================================================
//               Solve the Poisson equation
//---------------------------------------------------------------------------
void FemNonPeriodicPoissonSolver::solve_matrix_system(
        ChunkSpan<double, NUBSDomainX> phi_spline_coef,
        ChunkSpan<double, BSDomainX> rho_spline_coef) const
{
    double values_ptr[m_degree + 1];
    DSpan1D values(values_ptr, m_degree + 1);

    for (int i(0); i < m_nbasis; ++i) {
        phi_spline_coef(DiscreteElement<NUBSplinesX>(i)) = 0.0;
    }

    int const rhs_size = m_nbasis - 2;
    DSpan1D phi_rhs(phi_spline_coef.data() + 1, rhs_size);

    // Fill phi_rhs(i) with \int rho(x) b_i(x) dx
    // Rk: phi_rhs no longer contains spline coefficients, but is the
    //     RHS of the matrix equation
    for_each(m_quad_coef.domain(), [&](DiscreteElement<QMeshX> const ix) {
        Coordinate<RDimX> coord = coord_from_quad_point(coordinate(ix));
        auto jmin = discrete_space<BSplinesX>().eval_basis(values, coord);
        double const rho_val = m_spline_x_evaluator(coord, rho_spline_coef);
        for (int j = 0; j < m_degree + 1; ++j) {
            int const j_idx = (jmin.uid() + j) % m_nbasis - 1;
            if (j_idx != -1 && j_idx != rhs_size) {
                phi_rhs(j_idx) += rho_val * values(j) * m_quad_coef(ix);
            }
        }
    });

    // Solve the matrix equation to find the spline coefficients of phi
    m_fem_matrix->solve_inplace(phi_rhs);
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
void FemNonPeriodicPoissonSolver::operator()(
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
    Chunk<double, NUBSDomainX> phi_spline_coef(discrete_space<NUBSplinesX>().full_domain());
    solve_matrix_system(phi_spline_coef, rho_spline_coef);

    //
    for_each(dom_x, [&](IndexX const ix) {
        electrostatic_potential(ix) = m_spline_x_nu_evaluator(coordinate(ix), phi_spline_coef);
        electric_field(ix) = -m_spline_x_nu_evaluator.deriv(coordinate(ix), phi_spline_coef);
    });
}
