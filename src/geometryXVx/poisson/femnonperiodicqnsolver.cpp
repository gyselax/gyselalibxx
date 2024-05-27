// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <species_info.hpp>

#include "femnonperiodicqnsolver.hpp"

namespace {

//===========================================================================
// Create a non-uniform spline evaluator from the provided evaluator
//===========================================================================
NUBSplineXEvaluator_1d jit_build_nubsplinesx(SplineXEvaluator_1d const& spline_x_evaluator)
{
    static_assert(
            ddc::is_uniform_bsplines_v<
                    typename SplineXEvaluator_1d::
                            bsplines_type> || ddc::is_non_uniform_bsplines_v<typename SplineXEvaluator_1d::bsplines_type>);
    if constexpr (ddc::is_uniform_bsplines_v<typename SplineXEvaluator_1d::bsplines_type>) {
        int const ncells = ddc::discrete_space<BSplinesX>().ncells();
        std::vector<CoordX> knots(ncells + 1);

        for (int i(0); i < ncells + 1; ++i) {
            knots[i] = CoordX(ddc::discrete_space<BSplinesX>().get_knot(i));
        }
        ddc::init_discrete_space<NUBSplinesX>(knots);
    }
    // Boundary values are never evaluated
    return NUBSplineXEvaluator_1d(ddc::NullExtrapolationRule(), ddc::NullExtrapolationRule());
}

} // namespace

//===========================================================================
// Build an instance of FemNonPeriodicQNSolver
//===========================================================================
FemNonPeriodicQNSolver::FemNonPeriodicQNSolver(
        SplineXBuilder_1d const& spline_x_builder,
        SplineXEvaluator_1d const& spline_x_evaluator,
        IChargeDensityCalculator const& compute_rho)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_spline_x_nu_evaluator(jit_build_nubsplinesx(spline_x_evaluator))
    , m_compute_rho(compute_rho)
    , m_nbasis(ddc::discrete_space<BSplinesX>().nbasis())
    , m_ncells(ddc::discrete_space<BSplinesX>().ncells())
    , m_quad_coef_alloc(ddc::DiscreteDomain<QMeshX>(
              ddc::DiscreteElement<QMeshX>(0),
              ddc::DiscreteVector<QMeshX>(s_npts_gauss * m_ncells)))
{
    static_assert(!SplineXBuilder_1d::bsplines_type::is_periodic());
    BSDomainX const
            domain(ddc::DiscreteElement<BSplinesX>(0),
                   ddc::DiscreteVector<BSplinesX>(m_ncells + 1));
    ddc::Chunk<ddc::Coordinate<QDimX>, BSDomainX> knots(domain);

    for (ddc::DiscreteElement<BSplinesX> const i : domain) {
        knots(i) = quad_point_from_coord(ddc::discrete_space<NUBSplinesX>().get_knot(i.uid()));
    }

    // Calculate the integration coefficients
    GaussLegendre<QDimX> const gl(s_npts_gauss);
    std::vector<ddc::Coordinate<QDimX>> eval_pts_data(m_quad_coef_alloc.domain().size());
    ddc::ChunkSpan<ddc::Coordinate<QDimX>, ddc::DiscreteDomain<QMeshX>> const
            eval_pts(eval_pts_data.data(), m_quad_coef_alloc.domain());
    auto quad_coef_host = ddc::create_mirror_and_copy(m_quad_coef_alloc.span_view());
    gl.compute_points_and_weights_on_mesh(eval_pts, quad_coef_host.span_view(), knots.span_cview());
    ddc::parallel_deepcopy(m_quad_coef_alloc, quad_coef_host);

    ddc::init_discrete_space<QMeshX>(eval_pts_data);

    // Build the finite elements matrix
    build_matrix();
}

//===========================================================================
// Construct the finite element matrix with Bsplines as test functions
//===========================================================================
void FemNonPeriodicQNSolver::build_matrix()
{
    // Matrix contains all elements of the basis except the first
    // and last in order to impose Dirichlet conditions
    int const matrix_size = m_nbasis - 2;
    int constexpr n_lower_diags = s_degree;
    bool const positive_definite_symmetric = false;

    // Matrix with block is used instead of periodic to contain the
    // Dirichlet boundary conditions
    m_fem_matrix = Matrix::
            make_new_banded(matrix_size, n_lower_diags, n_lower_diags, positive_definite_symmetric);

    auto quad_coef_host = ddc::create_mirror_and_copy(m_quad_coef_alloc.span_cview());

    // Fill the banded part of the matrix
    std::array<double, s_degree + 1> derivs;
    ddc::for_each(m_quad_coef_alloc.domain(), [&](ddc::DiscreteElement<QMeshX> const ix) {
        ddc::Coordinate<RDimX> const coord = coord_from_quad_point(ddc::coordinate(ix));
        ddc::DiscreteElement<NUBSplinesX> const jmin
                = ddc::discrete_space<NUBSplinesX>().eval_deriv(derivs, coord);
        for (int j = 0; j < s_degree + 1; ++j) {
            for (int k = 0; k < s_degree + 1; ++k) {
                int const j_idx = (j + jmin.uid()) % m_nbasis - 1;
                int const k_idx = (k + jmin.uid()) % m_nbasis - 1;

                if (j_idx != -1 && j_idx != matrix_size && k_idx != -1 && k_idx != matrix_size) {
                    double a_jk = m_fem_matrix->get_element(j_idx, k_idx);
                    // Update element
                    a_jk += derivs[j] * derivs[k] * quad_coef_host(ix);

                    m_fem_matrix->set_element(j_idx, k_idx, a_jk);
                }
            }
        }
    });

    // Factorize the matrix ready to call solve
    m_fem_matrix->factorize();
}


//===========================================================================
//               Solve the Quasi-Neutrality equation
//---------------------------------------------------------------------------
void FemNonPeriodicQNSolver::solve_matrix_system(
        DNUBSSpanX const phi_spline_coef,
        DBSViewX const rho_spline_coef) const
{
    ddc::parallel_fill(phi_spline_coef, 0.0);

    int const rhs_size = m_nbasis - 2;
    int const nbasis_proxy = m_nbasis;
    SplineXEvaluator_1d spline_x_evaluator_proxy = m_spline_x_evaluator;
    ddc::ChunkSpan quad_coef = m_quad_coef_alloc.span_view();

    // Fill phi_rhs(i) with \int rho(x) b_i(x) dx
    // Rk: phi_rhs no longer contains spline coefficients, but is the
    //     RHS of the matrix equation
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            m_quad_coef_alloc.domain(),
            KOKKOS_LAMBDA(ddc::DiscreteElement<QMeshX> const ix) {
                ddc::Coordinate<RDimX> const coord = coord_from_quad_point(ddc::coordinate(ix));
                std::array<double, s_degree + 1> values;
                ddc::DiscreteElement<BSplinesX> const jmin
                        = ddc::discrete_space<BSplinesX>().eval_basis(values, coord);
                double const rho_val
                        = spline_x_evaluator_proxy(coord, rho_spline_coef.span_cview());
                for (int j = 0; j < s_degree + 1; ++j) {
                    int const j_idx = (jmin.uid() + j) % nbasis_proxy;
                    if (j_idx != 0 && j_idx != rhs_size + 1) {
                        Kokkos::atomic_fetch_add(
                                &phi_spline_coef(ddc::DiscreteElement<NUBSplinesX>(j_idx)),
                                rho_val * values[j] * quad_coef(ix));
                    }
                }
            });

    auto phi_spline_coef_host = ddc::create_mirror_and_copy(phi_spline_coef);
    DSpan1D const phi_rhs_host(phi_spline_coef_host.data_handle() + 1, rhs_size);

    // Solve the matrix equation to find the spline coefficients of phi
    m_fem_matrix->solve_inplace(phi_rhs_host);

    ddc::parallel_deepcopy(phi_spline_coef, phi_spline_coef_host);
}



//===========================================================================
// Resolution of Quasi-Neutrality equation through the finite difference method (FDM)
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
void FemNonPeriodicQNSolver::operator()(
        DSpanX electrostatic_potential,
        DSpanX electric_field,
        DViewSpXVx allfdistribu) const
{
    Kokkos::Profiling::pushRegion("QNSolver");
    assert(electrostatic_potential.domain() == ddc::get_domain<IDimX>(allfdistribu));
    IDomainX const dom_x = electrostatic_potential.domain();

    // Compute the RHS of the Quasi-Neutrality equation
    DFieldX rho(dom_x);
    m_compute_rho(rho, allfdistribu);

    //
    ddc::Chunk rho_spline_coef(m_spline_x_builder.spline_domain(), ddc::DeviceAllocator<double>());
    m_spline_x_builder(rho_spline_coef.span_view(), rho.span_cview());
    ddc::Chunk phi_spline_coef(
            ddc::discrete_space<NUBSplinesX>().full_domain(),
            ddc::DeviceAllocator<double>());
    solve_matrix_system(phi_spline_coef, rho_spline_coef.span_cview());

    //
    NUBSplineXEvaluator_1d spline_x_nu_evaluator_proxy = m_spline_x_nu_evaluator;
    ddc::ChunkSpan phi_spline_coef_view = phi_spline_coef.span_cview();
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_x,
            KOKKOS_LAMBDA(IndexX const ix) {
                electrostatic_potential(ix)
                        = spline_x_nu_evaluator_proxy(ddc::coordinate(ix), phi_spline_coef_view);
                electric_field(ix) = -spline_x_nu_evaluator_proxy
                                              .deriv(ddc::coordinate(ix), phi_spline_coef_view);
            });
    Kokkos::Profiling::popRegion();
}
