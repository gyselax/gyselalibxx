// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/matrix.hpp>

#include "ddc_aliases.hpp"
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
                  typename SplineBuilder::interpolation_discrete_dimension_type,
                  typename SplineEvaluator::evaluation_discrete_dimension_type>);

private:
    using base_type = IPoissonSolver<
            typename SplineEvaluator::evaluation_domain_type,
            typename SplineEvaluator::batched_evaluation_domain_type,
            std::experimental::layout_right,
            typename SplineEvaluator::memory_space>;

private:
    /// The interpolation mesh type
    using GridPDEDim = typename SplineBuilder::interpolation_discrete_dimension_type;

    using InputBSplines = typename SplineBuilder::bsplines_type;

    using PDEDim = typename GridPDEDim::continuous_dimension_type;

    using CoordPDEDim = Coord<PDEDim>;

    using IdxPDEDim = Idx<GridPDEDim>;

private:
    using field_type = typename base_type::field_type;

    using vector_field_type = typename base_type::vector_field_type;

    using fem_idx_range_type = typename base_type::laplacian_idx_range_type;

    using batch_idx_range_type = typename base_type::batch_idx_range_type;

    using batch_index_type = typename base_type::batch_index_type;

    using memory_space = typename base_type::memory_space;

    using exec_space = typename SplineEvaluator::exec_space;

public:
    /// The grid of quadrature points along the PDEDim direction
    struct GridPDEDimQ : NonUniformGridBase<PDEDim>
    {
    };

private:
    /// An index of the grid of quadrature points.
    using IdxQ = Idx<GridPDEDimQ>;

    /// A displacement on the grid of quadrature points.
    using IdxStepQ = IdxStep<GridPDEDimQ>;

    /// The index range of the grid of quadrature points.
    using IdxRangeQ = IdxRange<GridPDEDimQ>;

    /**
     * @brief An array of values defined at the quadrature points.
     *
     * @tparam ElementType The type of the elements in the chunk.
     */
    using DQField = DFieldMem<IdxRangeQ, ddc::KokkosAllocator<double, memory_space>>;

    /**
     * @brief An array of coordinates defined at the quadrature points.
     *
     * @tparam ElementType The type of the elements in the chunk.
     */
    using CoordQField
            = Field<CoordPDEDim, IdxRangeQ, std::experimental::layout_right, memory_space>;

public:
    /**
     * The internal type of the FEM basis. This is used if the provided BSplines are uniform as
     * non-uniform BSplines are required to correctly enforce Dirichlet boundary conditions.
     *
     * This type should be private but is public due to Kokkos restrictions.
     */
    struct HiddenFEMBSplines : ddc::NonUniformBSplines<PDEDim, InputBSplines::degree()>
    {
    };

private:
    using FEMBSplines = std::conditional_t<
            ddc::is_uniform_bsplines_v<InputBSplines> && !InputBSplines::is_periodic(),
            HiddenFEMBSplines,
            InputBSplines>;

    using FEMEvalExtrapolationRule = std::conditional_t<
            FEMBSplines::is_periodic(),
            ddc::PeriodicExtrapolationRule<PDEDim>,
            ddc::NullExtrapolationRule>;

    template <class IdxRange>
    struct FEMSplineEvaluatorBuilder;

    template <class... DimX>
    struct FEMSplineEvaluatorBuilder<IdxRange<DimX...>>
    {
        using type = ddc::SplineEvaluator<
                Kokkos::DefaultExecutionSpace,
                Kokkos::DefaultExecutionSpace::memory_space,
                FEMBSplines,
                GridPDEDim,
                FEMEvalExtrapolationRule,
                FEMEvalExtrapolationRule,
                DimX...>;
    };

    using FEMSplineEvaluator = typename FEMSplineEvaluatorBuilder<
            typename SplineEvaluator::batched_evaluation_domain_type>::type;

    using FEMBSplinesIdx = Idx<FEMBSplines>;

    using FEMBSplinesIdxRange = IdxRange<FEMBSplines>;

    using BatchedFEMBSplinesIdxRange = typename FEMSplineEvaluator::batched_spline_domain_type;

    using FEMBSplinesCoeffMem
            = DFieldMem<FEMBSplinesIdxRange, ddc::KokkosAllocator<double, memory_space>>;

    using BatchedFEMBSplinesCoeffMem
            = DFieldMem<BatchedFEMBSplinesIdxRange, ddc::KokkosAllocator<double, memory_space>>;

    using BatchedFEMBSplinesCoeff = typename BatchedFEMBSplinesCoeffMem::span_type;

private:
    using BSplinesIdxRange = IdxRange<InputBSplines>;
    using BSplinesIdx = Idx<InputBSplines>;

    using BatchedBSplinesIdxRange = typename SplineEvaluator::batched_spline_domain_type;

    using full_index =
            typename SplineEvaluator::batched_evaluation_domain_type::discrete_element_type;

    using CoordField = FieldMem<
            CoordPDEDim,
            typename FEMSplineEvaluator::batched_evaluation_domain_type,
            ddc::KokkosAllocator<CoordPDEDim, memory_space>>;

private:
    using RHSBSplines = InputBSplines;

    using RHSSplineCoeffMem
            = DFieldMem<BatchedBSplinesIdxRange, ddc::KokkosAllocator<double, memory_space>>;
    using RHSSplineCoeff = typename RHSSplineCoeffMem::span_type;

    using RHSQuadTags = ddc::
            type_seq_merge_t<typename base_type::batch_tags, ddc::detail::TypeSeq<GridPDEDimQ>>;

    using RHSQuadratureIdxRange = ddc::detail::convert_type_seq_to_discrete_domain_t<RHSQuadTags>;

    using RHSQuadratureIdx = typename RHSQuadratureIdxRange::discrete_element_type;

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
        , m_quad_coef(IdxRangeQ(
                  IdxQ(0),
                  IdxStepQ(s_npts_gauss * ddc::discrete_space<InputBSplines>().ncells())))
    {
        using break_point_grid = ddc::knot_discrete_dimension_t<InputBSplines>;
        IdxRange<break_point_grid> break_point_idx_range
                = ddc::discrete_space<InputBSplines>().break_point_domain();
        host_t<FieldMem<CoordPDEDim, IdxRange<break_point_grid>>> break_points(
                break_point_idx_range);

        for (Idx<break_point_grid> const idx : break_point_idx_range) {
            break_points(idx) = ddc::coordinate(idx);
        }

        // Calculate the integration coefficients
        GaussLegendre<PDEDim> const gl(s_npts_gauss);
        std::vector<CoordPDEDim> eval_pts_data(get_idx_range(m_quad_coef).size());
        host_t<CoordQField> const eval_pts(eval_pts_data.data(), get_idx_range(m_quad_coef));
        auto quad_coef_host = ddc::create_mirror_and_copy(get_field(m_quad_coef));
        gl.compute_points_and_weights_on_mesh(
                eval_pts,
                get_field(quad_coef_host),
                get_const_field(break_points));
        ddc::parallel_deepcopy(m_quad_coef, quad_coef_host);

        ddc::init_discrete_space<GridPDEDimQ>(eval_pts_data);
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
    field_type operator()(field_type phi, field_type rho) const override
    {
        batch_idx_range_type batch_dom(get_idx_range(phi));
        BatchedFEMBSplinesIdxRange
                phi_coefs_idx_range(batch_dom, ddc::discrete_space<FEMBSplines>().full_domain());
        BatchedFEMBSplinesCoeffMem phi_coefs_alloc(phi_coefs_idx_range);
        BatchedFEMBSplinesCoeff phi_coefs(phi_coefs_alloc);
        solve_matrix_system(phi_coefs, rho);

        CoordField eval_pts_alloc(get_idx_range(phi));
        ddc::ChunkSpan eval_pts = get_field(eval_pts_alloc);

        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(eval_pts),
                KOKKOS_LAMBDA(full_index const idx) {
                    eval_pts(idx) = ddc::coordinate(ddc::select<GridPDEDim>(idx));
                });

        FEMSplineEvaluator spline_nu_evaluator_proxy = m_spline_fem_evaluator;

        m_spline_fem_evaluator(phi, get_const_field(eval_pts), get_const_field(phi_coefs));

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
    field_type operator()(field_type phi, vector_field_type E, field_type rho) const override
    {
        batch_idx_range_type batch_dom(get_idx_range(phi));
        BatchedFEMBSplinesIdxRange
                phi_coefs_idx_range(batch_dom, ddc::discrete_space<FEMBSplines>().full_domain());
        BatchedFEMBSplinesCoeffMem phi_coefs_alloc(phi_coefs_idx_range);
        BatchedFEMBSplinesCoeff phi_coefs(phi_coefs_alloc);
        solve_matrix_system(phi_coefs, rho);

        CoordField eval_pts_alloc(get_idx_range(phi));
        ddc::ChunkSpan eval_pts = get_field(eval_pts_alloc);

        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(eval_pts),
                KOKKOS_LAMBDA(full_index const idx) {
                    eval_pts(idx) = ddc::coordinate(ddc::select<GridPDEDim>(idx));
                });

        FEMSplineEvaluator spline_nu_evaluator_proxy = m_spline_fem_evaluator;

        m_spline_fem_evaluator(phi, get_const_field(eval_pts), get_const_field(phi_coefs));
        m_spline_fem_evaluator.deriv(E, get_const_field(eval_pts), get_const_field(phi_coefs));

        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(phi),
                KOKKOS_LAMBDA(full_index const idx) { E(idx) = -E(idx); });
        return phi;
    }

private:
    void build_periodic_matrix()
    {
        int constexpr n_lower_diags = s_degree + 1;
        int const nbasis = ddc::discrete_space<FEMBSplines>().nbasis();
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

        auto quad_coef_host = ddc::create_mirror_and_copy(get_const_field(m_quad_coef));

        // Fill the banded part of the matrix
        std::array<double, s_degree + 1> derivs_alloc;
        DSpan1D derivs = as_span(derivs_alloc);
        ddc::for_each(get_idx_range(m_quad_coef), [&](IdxQ const ix) {
            CoordPDEDim const coord = ddc::coordinate(ix);
            FEMBSplinesIdx const jmin
                    = ddc::discrete_space<FEMBSplines>().eval_deriv(derivs, coord);
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
        FEMBSplinesIdxRange const bspline_full_idx_range
                = ddc::discrete_space<FEMBSplines>().full_domain();
        FEMBSplinesIdxRange const bspline_dom
                = bspline_full_idx_range.take_first(IdxStep<FEMBSplines>(nbasis));

        host_t<FEMBSplinesCoeffMem> int_vals(bspline_dom);
        ddc::discrete_space<InputBSplines>().integrals(get_field(int_vals));

        for (FEMBSplinesIdx const ix : bspline_dom) {
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
        int const nbasis = ddc::discrete_space<FEMBSplines>().nbasis();
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

        auto quad_coef_host = ddc::create_mirror_and_copy(get_const_field(m_quad_coef));

        // Fill the banded part of the matrix
        std::array<double, s_degree + 1> derivs_alloc;
        DSpan1D derivs = as_span(derivs_alloc);
        ddc::for_each(get_idx_range(m_quad_coef), [&](IdxQ const ix) {
            CoordPDEDim const coord = ddc::coordinate(ix);
            FEMBSplinesIdx const jmin
                    = ddc::discrete_space<FEMBSplines>().eval_deriv(derivs, coord);
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
    void solve_matrix_system(BatchedFEMBSplinesCoeff phi_spline_coef, field_type rho) const
    {
        // Calculate the spline representation of the RHS.
        BatchedBSplinesIdxRange rho_spline_coef_idx_range(
                batch_idx_range_type(get_idx_range(rho)),
                get_spline_idx_range(m_spline_builder));
        RHSSplineCoeffMem rho_spline_coef_alloc(rho_spline_coef_idx_range);
        RHSSplineCoeff rho_spline_coef(rho_spline_coef_alloc);
        m_spline_builder(rho_spline_coef, get_const_field(rho));

        IdxRange<FEMBSplines> fem_idx_range = ddc::discrete_space<FEMBSplines>().full_domain();
        int const nbasis_proxy = ddc::discrete_space<FEMBSplines>().nbasis();
        SplineEvaluator spline_evaluator_proxy = m_spline_evaluator;
        ddc::ChunkSpan quad_coef_proxy = get_field(m_quad_coef);

        FEMBSplinesIdx last_basis_element(nbasis_proxy - 1);

        ddc::parallel_fill(phi_spline_coef, 0.0);

        // Create the rhs as an alias for phi_spline_coef as the matrix equation is solved in place.
        ddc::ChunkSpan rhs(phi_spline_coef);

        batch_idx_range_type batch_idx_range(get_idx_range(rho));
        RHSQuadratureIdxRange rhs_build_idx_range(batch_idx_range, get_idx_range(m_quad_coef));

        // Fill phi_rhs(i) with \int rho(x) b_i(x) dx
        // Rk: phi_rhs no longer contains spline coefficients, but is the
        //     RHS of the matrix equation
        ddc::parallel_for_each(
                exec_space(),
                rhs_build_idx_range,
                KOKKOS_LAMBDA(RHSQuadratureIdx const idx) {
                    batch_index_type ib(idx);
                    IdxQ iq(idx);
                    CoordPDEDim const coord = ddc::coordinate(iq);
                    std::array<double, s_degree + 1> values_alloc;
                    DSpan1D values = as_span(values_alloc);
                    FEMBSplinesIdx const jmin
                            = ddc::discrete_space<FEMBSplines>().eval_basis(values, coord);
                    double const rho_val
                            = spline_evaluator_proxy(coord, rho_spline_coef[ib].span_cview());
                    for (int j = 0; j < s_degree + 1; ++j) {
                        FEMBSplinesIdx j_idx = jmin + j;
                        while (j_idx > last_basis_element) {
                            j_idx -= nbasis_proxy;
                        }
                        Kokkos::atomic_fetch_add(
                                &rhs(ib, j_idx),
                                rho_val * values[j] * quad_coef_proxy(iq));
                    }
                });

        auto phi_spline_coef_host_alloc = ddc::create_mirror_and_copy(phi_spline_coef);
        ddc::ChunkSpan phi_spline_coef_host = get_field(phi_spline_coef_host_alloc);

        int constexpr n_implicit_min_bcs(!InputBSplines::is_periodic());

        ddc::for_each(batch_idx_range, [&](batch_index_type ib) {
            FEMBSplinesIdxRange solve_idx_range(
                    fem_idx_range.front() + n_implicit_min_bcs,
                    IdxStep<FEMBSplines>(m_matrix_size));
            DSpan1D const phi_rhs_host
                    = phi_spline_coef_host[ib][solve_idx_range].allocation_mdspan();

            // Solve the matrix equation to find the spline coefficients of phi
            m_fem_matrix->solve_inplace(phi_rhs_host);
        });

        ddc::parallel_deepcopy(phi_spline_coef, phi_spline_coef_host);

        if constexpr (!InputBSplines::is_periodic()) {
            // Apply Dirichlet BCs
            ddc::parallel_fill(phi_spline_coef[fem_idx_range.front()], 0.0);
            ddc::parallel_fill(phi_spline_coef[fem_idx_range.back()], 0.0);
        } else {
            FEMBSplinesIdx first_repeat_bspline(ddc::discrete_space<FEMBSplines>().nbasis());
            // Copy the first d coefficients into the last d coefficients
            // These coefficients refer to the same InputBSplines which cross the boundaries
            ddc::parallel_for_each(
                    exec_space(),
                    batch_idx_range,
                    KOKKOS_LAMBDA(batch_index_type ib) {
                        for (int i = 0; i < s_degree; i++) {
                            phi_spline_coef(ib, first_repeat_bspline + i)
                                    = phi_spline_coef(ib, FEMBSplinesIdx(i));
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
                using break_point_grid = ddc::knot_discrete_dimension_t<InputBSplines>;
                IdxRange<break_point_grid> break_point_idx_range
                        = ddc::discrete_space<InputBSplines>().break_point_domain();
                std::vector<CoordPDEDim> break_points(break_point_idx_range.size());

                for (Idx<break_point_grid> const idx : break_point_idx_range) {
                    break_points[(idx - break_point_idx_range.front()).value()]
                            = ddc::coordinate(idx);
                }
                ddc::init_discrete_space<FEMBSplines>(break_points);
            }
            // Boundary values are never evaluated
            return FEMSplineEvaluator(ddc::NullExtrapolationRule(), ddc::NullExtrapolationRule());
        }
    }
};
