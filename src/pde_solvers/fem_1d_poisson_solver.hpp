// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/matrix.hpp>

#include "ipoisson_solver.hpp"


/**
 * A class to solve the following equation:
 * @f$ -\Delta \phi = \rho @f$
 * using a Finite Element Method.
 *
 * @tparam SplineEvaluator An evaluator which can be used to evaluate splines.
 */
template <class SplineBuilder, class SplineEvaluator>
class FEM1DPoissonSolver
    : public IPoissonSolver<
              typename SplineEvaluator::evaluation_domain_type,
              typename SplineEvaluator::batched_evaluation_domain_type,
              std::experimental::layout_right,
              typename SplineEvaluator::memory_space>
{
    static_assert(std::is_same_v<
                  typename SplineBuilder::batched_interpolation_domain_type,
                  typename SplineEvaluator::batched_evaluation_domain_type>);
    static_assert(std::is_same_v<
                  typename SplineBuilder::interpolation_mesh_type,
                  typename SplineEvaluator::evaluation_mesh_type>);

private:
    using base_type = IPoissonSolver<
            typename SplineEvaluator::evaluation_domain_type,
            typename SplineEvaluator::batched_evaluation_domain_type,
            std::experimental::layout_right,
            typename SplineEvaluator::memory_space>;

private:
    /// The interpolation mesh type
    using DimI = typename SplineBuilder::interpolation_mesh_type;

    using InputBSplines = typename SplineBuilder::bsplines_type;

    using RDimI = typename DimI::continuous_dimension_type;

    using CoordI = ddc::Coordinate<RDimI>;

    using IndexI = ddc::DiscreteElement<DimI>;

private:
    using chunk_span_type = typename base_type::chunk_span_type;

    using vector_span_type = typename base_type::vector_span_type;

    using fem_domain_type = typename base_type::laplacian_domain_type;

    using batch_domain_type = typename base_type::batch_domain_type;

    using batch_element_type = typename base_type::batch_element_type;

    using memory_space = typename base_type::memory_space;

    using exec_space = typename SplineEvaluator::exec_space;

public:
    struct QMeshI : ddc::NonUniformPointSampling<RDimI>
    {
    };

private:
    /// An index of the grid of quadrature points.
    using IndexQ = ddc::DiscreteElement<QMeshI>;

    /// A displacement on the grid of quadrature points.
    using VectQ = ddc::DiscreteVector<QMeshI>;

    /// The domain of the grid of quadrature points.
    using DomainQ = ddc::DiscreteDomain<QMeshI>;

    /**
     * @brief An array of values defined at the quadrature points.
     *
     * @tparam ElementType The type of the elements in the chunk.
     */
    using DQField = ddc::Chunk<double, DomainQ, ddc::KokkosAllocator<double, memory_space>>;

    /**
     * @brief An array of coordinates defined at the quadrature points.
     *
     * @tparam ElementType The type of the elements in the chunk.
     */
    using CoordQSpan
            = ddc::ChunkSpan<CoordI, DomainQ, std::experimental::layout_right, memory_space>;

public:
    /**
     * The internal type of the FEM basis. This is used if the provided BSplines are uniform as
     * non-uniform BSplines are required to correctly enforce Dirichlet boundary conditions.
     *
     * This type should be private but is public due to Kokkos restrictions.
     */
    struct HiddenFEMBasis : ddc::NonUniformBSplines<RDimI, InputBSplines::degree()>
    {
    };

private:
    using FEMBasis = std::conditional_t<
            ddc::is_uniform_bsplines_v<InputBSplines> && !InputBSplines::is_periodic(),
            HiddenFEMBasis,
            InputBSplines>;

    using FEMEvalExtrapolationRule = std::conditional_t<
            FEMBasis::is_periodic(),
            ddc::PeriodicExtrapolationRule<RDimI>,
            ddc::NullExtrapolationRule>;

    template <class Domain>
    struct FEMSplineEvaluatorBuilder;

    template <class... DimX>
    struct FEMSplineEvaluatorBuilder<ddc::DiscreteDomain<DimX...>>
    {
        using type = ddc::SplineEvaluator<
                Kokkos::DefaultExecutionSpace,
                Kokkos::DefaultExecutionSpace::memory_space,
                FEMBasis,
                DimI,
                FEMEvalExtrapolationRule,
                FEMEvalExtrapolationRule,
                DimX...>;
    };

    using FEMSplineEvaluator = typename FEMSplineEvaluatorBuilder<
            typename SplineEvaluator::batched_evaluation_domain_type>::type;

    using FEMBasisIndex = ddc::DiscreteElement<FEMBasis>;

    using FEMBSDomain = ddc::DiscreteDomain<FEMBasis>;

    using FEMBasisCoeffDomainType = typename FEMSplineEvaluator::batched_spline_domain_type;

    using FEMBSField = ddc::Chunk<double, FEMBSDomain, ddc::KokkosAllocator<double, memory_space>>;

    using FEMBasisCoeffType = ddc::
            Chunk<double, FEMBasisCoeffDomainType, ddc::KokkosAllocator<double, memory_space>>;

    using FEMBasisCoeffSpan = typename FEMBasisCoeffType::span_type;

private:
    using BSDomain = ddc::DiscreteDomain<InputBSplines>;
    using BSIndex = ddc::DiscreteElement<InputBSplines>;
    using BSVect = ddc::DiscreteVector<InputBSplines>;

    using BSCoeffDomainType = typename SplineEvaluator::batched_spline_domain_type;

    using full_index =
            typename SplineEvaluator::batched_evaluation_domain_type::discrete_element_type;

    using CoordField = ddc::Chunk<
            CoordI,
            typename FEMSplineEvaluator::batched_evaluation_domain_type,
            ddc::KokkosAllocator<CoordI, memory_space>>;

private:
    using RHSBsplines = InputBSplines;

    using RHSSplineType
            = ddc::Chunk<double, BSCoeffDomainType, ddc::KokkosAllocator<double, memory_space>>;
    using RHSSplineSpan = typename RHSSplineType::span_type;

    using RHSQuadTags
            = ddc::type_seq_merge_t<typename base_type::batch_tags, ddc::detail::TypeSeq<QMeshI>>;

    using RHSQuadratureDomain = ddc::detail::convert_type_seq_to_discrete_domain<RHSQuadTags>;

    using RHSQuadratureIndex = typename RHSQuadratureDomain::discrete_element_type;

private:
    // Spline degree in x direction
    static int constexpr s_degree = InputBSplines::degree();

    // Gauss points used for integration computation
    static int constexpr s_npts_gauss = InputBSplines::degree() + 1;

private:
    SplineBuilder const& m_spline_builder;

    SplineEvaluator m_spline_evaluator;

    FEMSplineEvaluator m_spline_fem_evaluator;

    int m_matrix_size;

    DQField m_quad_coef;

    std::unique_ptr<Matrix> m_fem_matrix;

public:
    /**
     * Construct the FemQNSolver operator.
     *
     * @param spline_builder A spline builder which calculates the coefficients of a spline representation.
     * @param spline_evaluator A spline evaluator which provides the value of a spline representation from its coefficients.
     */
    FEM1DPoissonSolver(SplineBuilder const& spline_builder, SplineEvaluator const& spline_evaluator)
        : m_spline_builder(spline_builder)
        , m_spline_evaluator(spline_evaluator)
        , m_spline_fem_evaluator(jit_build_nubsplinesx(spline_evaluator))
        , m_quad_coef(
                  DomainQ(IndexQ(0),
                          VectQ(s_npts_gauss * ddc::discrete_space<InputBSplines>().ncells())))
    {
        ddc::DiscreteDomain break_point_domain
                = ddc::discrete_space<InputBSplines>().break_point_domain();
        ddc::Chunk break_points(break_point_domain, ddc::HostAllocator<CoordI>());

        for (ddc::DiscreteElement const idx : break_point_domain) {
            break_points(idx) = ddc::coordinate(idx);
        }

        // Calculate the integration coefficients
        GaussLegendre<RDimI> const gl(s_npts_gauss);
        std::vector<CoordI> eval_pts_data(m_quad_coef.domain().size());
        host_t<CoordQSpan> const eval_pts(eval_pts_data.data(), m_quad_coef.domain());
        auto quad_coef_host = ddc::create_mirror_and_copy(m_quad_coef.span_view());
        gl.compute_points_and_weights_on_mesh(
                eval_pts,
                quad_coef_host.span_view(),
                break_points.span_cview());
        ddc::parallel_deepcopy(m_quad_coef, quad_coef_host);

        ddc::init_discrete_space<QMeshI>(eval_pts_data);
        // Build the finite elements matrix
        if constexpr (InputBSplines::is_periodic()) {
            build_periodic_matrix();
        } else {
            build_non_periodic_matrix();
        }
    }

    /**
     * @brief An operator which calculates the solution @f$\phi@f$ to Poisson's equation:
     * @f$ -\Delta \phi = \rho @f$
     *
     * @param[out] phi The solution to Poisson's equation.
     * @param[in] rho The right-hand side of Poisson's equation.
     *
     * @return A reference to the solution to Poisson's equation.
     */
    chunk_span_type operator()(chunk_span_type phi, chunk_span_type rho) const override
    {
        batch_domain_type batch_dom(phi.domain());
        FEMBasisCoeffDomainType
                phi_coefs_domain(batch_dom, ddc::discrete_space<FEMBasis>().full_domain());
        FEMBasisCoeffType phi_coefs_alloc(phi_coefs_domain);
        FEMBasisCoeffSpan phi_coefs(phi_coefs_alloc);
        solve_matrix_system(phi_coefs, rho);

        CoordField eval_pts_alloc(phi.domain());
        ddc::ChunkSpan eval_pts = eval_pts_alloc.span_view();

        ddc::parallel_for_each(
                exec_space(),
                eval_pts.domain(),
                KOKKOS_LAMBDA(full_index const idx) {
                    eval_pts(idx) = ddc::coordinate(ddc::select<DimI>(idx));
                });

        FEMSplineEvaluator spline_nu_evaluator_proxy = m_spline_fem_evaluator;

        m_spline_fem_evaluator(phi, eval_pts.span_cview(), phi_coefs.span_cview());

        return phi;
    }

    /**
     * @brief An operator which calculates the solution @f$\phi@f$ to Poisson's equation and
     * its derivative:
     * @f$ - \Delta \phi = \rho @f$
     * @f$ E = - \nabla \phi @f$
     *
     * @param[out] phi The solution to Poisson's equation.
     * @param[out] E The negative derivative of the solution to Poisson's equation.
     * @param[in] rho The right-hand side of Poisson's equation.
     *
     * @return A reference to the solution to Poisson's equation.
     */
    chunk_span_type operator()(chunk_span_type phi, vector_span_type E, chunk_span_type rho)
            const override
    {
        batch_domain_type batch_dom(phi.domain());
        FEMBasisCoeffDomainType
                phi_coefs_domain(batch_dom, ddc::discrete_space<FEMBasis>().full_domain());
        FEMBasisCoeffType phi_coefs_alloc(phi_coefs_domain);
        FEMBasisCoeffSpan phi_coefs(phi_coefs_alloc);
        solve_matrix_system(phi_coefs, rho);

        CoordField eval_pts_alloc(phi.domain());
        ddc::ChunkSpan eval_pts = eval_pts_alloc.span_view();

        ddc::parallel_for_each(
                exec_space(),
                eval_pts.domain(),
                KOKKOS_LAMBDA(full_index const idx) {
                    eval_pts(idx) = ddc::coordinate(ddc::select<DimI>(idx));
                });

        FEMSplineEvaluator spline_nu_evaluator_proxy = m_spline_fem_evaluator;

        m_spline_fem_evaluator(phi, eval_pts.span_cview(), phi_coefs.span_cview());
        m_spline_fem_evaluator.deriv(E, eval_pts.span_cview(), phi_coefs.span_cview());

        ddc::parallel_for_each(
                exec_space(),
                phi.domain(),
                KOKKOS_LAMBDA(full_index const idx) { E(idx) = -E(idx); });
        return phi;
    }

private:
    void build_periodic_matrix()
    {
        int constexpr n_lower_diags = s_degree + 1;
        int const nbasis = ddc::discrete_space<FEMBasis>().nbasis();
        m_matrix_size = nbasis + 1;
        bool const positive_definite_symmetric = false;

        // Matrix with block is used instead of periodic to contain the
        // Dirichlet boundary conditions
        m_fem_matrix = Matrix::make_new_block_with_banded_region(
                m_matrix_size,
                n_lower_diags,
                n_lower_diags,
                positive_definite_symmetric,
                n_lower_diags);

        auto quad_coef_host = ddc::create_mirror_and_copy(m_quad_coef.span_cview());

        // Fill the banded part of the matrix
        std::array<double, s_degree + 1> derivs_alloc;
        DSpan1D derivs = as_span(derivs_alloc);
        ddc::for_each(m_quad_coef.domain(), [&](IndexQ const ix) {
            CoordI const coord = ddc::coordinate(ix);
            FEMBasisIndex const jmin = ddc::discrete_space<FEMBasis>().eval_deriv(derivs, coord);
            for (int j = 0; j < s_degree + 1; ++j) {
                for (int k = 0; k < s_degree + 1; ++k) {
                    int const j_idx = (j + jmin.uid()) % nbasis;
                    int const k_idx = (k + jmin.uid()) % nbasis;
                    double a_jk = m_fem_matrix->get_element(j_idx, k_idx);
                    // Update element
                    a_jk += derivs[j] * derivs[k] * quad_coef_host(ix);

                    m_fem_matrix->set_element(j_idx, k_idx, a_jk);
                }
            }
        });

        // Impose the boundary conditions
        FEMBSDomain const bspline_full_domain = ddc::discrete_space<FEMBasis>().full_domain();
        FEMBSDomain const bspline_dom
                = bspline_full_domain.take_first(ddc::DiscreteVector<FEMBasis>(nbasis));

        host_t<FEMBSField> int_vals(bspline_dom);
        ddc::discrete_space<InputBSplines>().integrals(int_vals.span_view());

        for (FEMBasisIndex const ix : bspline_dom) {
            int const i = ix.uid();
            m_fem_matrix->set_element(nbasis, i, int_vals(ix));
            m_fem_matrix->set_element(i, nbasis, int_vals(ix));
        }

        // Factorize the matrix ready to call solve
        m_fem_matrix->factorize();
    }

    void build_non_periodic_matrix()
    {
        // Matrix contains all elements of the basis except the first
        // and last in order to impose Dirichlet conditions
        int const nbasis = ddc::discrete_space<FEMBasis>().nbasis();
        m_matrix_size = nbasis - 2;
        int constexpr n_lower_diags = s_degree;
        bool const positive_definite_symmetric = false;

        // Matrix with block is used instead of periodic to contain the
        // Dirichlet boundary conditions
        m_fem_matrix = Matrix::make_new_banded(
                m_matrix_size,
                n_lower_diags,
                n_lower_diags,
                positive_definite_symmetric);

        auto quad_coef_host = ddc::create_mirror_and_copy(m_quad_coef.span_cview());

        // Fill the banded part of the matrix
        std::array<double, s_degree + 1> derivs_alloc;
        DSpan1D derivs = as_span(derivs_alloc);
        ddc::for_each(m_quad_coef.domain(), [&](IndexQ const ix) {
            CoordI const coord = ddc::coordinate(ix);
            FEMBasisIndex const jmin = ddc::discrete_space<FEMBasis>().eval_deriv(derivs, coord);
            for (int j = 0; j < s_degree + 1; ++j) {
                for (int k = 0; k < s_degree + 1; ++k) {
                    int const j_idx = (j + jmin.uid()) % nbasis - 1;
                    int const k_idx = (k + jmin.uid()) % nbasis - 1;

                    if (j_idx != -1 && j_idx != m_matrix_size && k_idx != -1
                        && k_idx != m_matrix_size) {
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


public:
    /**
     * [SHOULD BE PRIVATE (Kokkos limitation)]
     *
     * @param[out] phi_spline_coef The spline coefficients which will describe the result on the FEM basis.
     * @param[in] rho The function on the right hand side of the equation at the interpolation points.
     */
    void solve_matrix_system(FEMBasisCoeffSpan phi_spline_coef, chunk_span_type rho) const
    {
        // Calculate the spline representation of the RHS.
        BSCoeffDomainType rho_spline_coef_domain(
                batch_domain_type(rho.domain()),
                m_spline_builder.spline_domain());
        RHSSplineType rho_spline_coef_alloc(rho_spline_coef_domain);
        RHSSplineSpan rho_spline_coef(rho_spline_coef_alloc);
        m_spline_builder(rho_spline_coef, rho.span_cview());

        ddc::DiscreteDomain<FEMBasis> fem_domain = ddc::discrete_space<FEMBasis>().full_domain();
        int const nbasis_proxy = ddc::discrete_space<FEMBasis>().nbasis();
        SplineEvaluator spline_evaluator_proxy = m_spline_evaluator;
        ddc::ChunkSpan quad_coef_proxy = m_quad_coef.span_view();

        FEMBasisIndex last_basis_element(nbasis_proxy - 1);

        ddc::parallel_fill(phi_spline_coef, 0.0);

        // Create the rhs as an alias for phi_spline_coef as the matrix equation is solved in place.
        ddc::ChunkSpan rhs(phi_spline_coef);

        batch_domain_type batch_domain(rho.domain());
        RHSQuadratureDomain rhs_build_domain(batch_domain, m_quad_coef.domain());

        // Fill phi_rhs(i) with \int rho(x) b_i(x) dx
        // Rk: phi_rhs no longer contains spline coefficients, but is the
        //     RHS of the matrix equation
        ddc::parallel_for_each(
                exec_space(),
                rhs_build_domain,
                KOKKOS_LAMBDA(RHSQuadratureIndex const idx) {
                    batch_element_type ib(idx);
                    IndexQ iq(idx);
                    CoordI const coord = ddc::coordinate(iq);
                    std::array<double, s_degree + 1> values_alloc;
                    DSpan1D values = as_span(values_alloc);
                    FEMBasisIndex const jmin
                            = ddc::discrete_space<FEMBasis>().eval_basis(values, coord);
                    double const rho_val
                            = spline_evaluator_proxy(coord, rho_spline_coef[ib].span_cview());
                    for (int j = 0; j < s_degree + 1; ++j) {
                        FEMBasisIndex j_idx = jmin + j;
                        while (j_idx > last_basis_element) {
                            j_idx -= nbasis_proxy;
                        }
                        Kokkos::atomic_fetch_add(
                                &rhs(ib, j_idx),
                                rho_val * values[j] * quad_coef_proxy(iq));
                    }
                });

        auto phi_spline_coef_host_alloc = ddc::create_mirror_and_copy(phi_spline_coef);
        ddc::ChunkSpan phi_spline_coef_host = phi_spline_coef_host_alloc.span_view();

        int constexpr n_implicit_min_bcs(!InputBSplines::is_periodic());

        ddc::for_each(batch_domain, [&](batch_element_type ib) {
            FEMBSDomain solve_domain(
                    fem_domain.front() + n_implicit_min_bcs,
                    ddc::DiscreteVector<FEMBasis>(m_matrix_size));
            DSpan1D const phi_rhs_host = phi_spline_coef_host[ib][solve_domain].allocation_mdspan();

            // Solve the matrix equation to find the spline coefficients of phi
            m_fem_matrix->solve_inplace(phi_rhs_host);
        });

        ddc::parallel_deepcopy(phi_spline_coef, phi_spline_coef_host);

        if constexpr (!InputBSplines::is_periodic()) {
            // Apply Dirichlet BCs
            ddc::parallel_fill(phi_spline_coef[fem_domain.front()], 0.0);
            ddc::parallel_fill(phi_spline_coef[fem_domain.back()], 0.0);
        } else {
            FEMBasisIndex first_repeat_bspline(ddc::discrete_space<FEMBasis>().nbasis());
            // Copy the first d coefficients into the last d coefficients
            // These coefficients refer to the same InputBSplines which cross the boundaries
            ddc::parallel_for_each(
                    exec_space(),
                    batch_domain,
                    KOKKOS_LAMBDA(batch_element_type ib) {
                        for (int i = 0; i < s_degree; i++) {
                            phi_spline_coef(ib, first_repeat_bspline + i)
                                    = phi_spline_coef(ib, FEMBasisIndex(i));
                        }
                    });
        }
    }

private:
    //===========================================================================
    // Create a non-uniform spline evaluator from the provided evaluator
    //===========================================================================
    static FEMSplineEvaluator jit_build_nubsplinesx(SplineEvaluator const& spline_evaluator)
    {
        if constexpr (InputBSplines::is_periodic()) {
            return spline_evaluator;
        } else {
            static_assert(
                    (ddc::is_uniform_bsplines_v<typename SplineEvaluator::bsplines_type>)
                    || ddc::is_non_uniform_bsplines_v<typename SplineEvaluator::bsplines_type>);
            if constexpr (ddc::is_uniform_bsplines_v<typename SplineEvaluator::bsplines_type>) {
                ddc::DiscreteDomain break_point_domain
                        = ddc::discrete_space<InputBSplines>().break_point_domain();
                std::vector<CoordI> break_points(break_point_domain.size());

                for (ddc::DiscreteElement const idx : break_point_domain) {
                    break_points[(idx - break_point_domain.front()).value()] = ddc::coordinate(idx);
                }
                ddc::init_discrete_space<FEMBasis>(break_points);
            }
            // Boundary values are never evaluated
            return FEMSplineEvaluator(ddc::NullExtrapolationRule(), ddc::NullExtrapolationRule());
        }
    }
};
