// SPDX-License-Identifier: MIT
#pragma once

#include "polarpoissonlikeassembler.hpp"

/**
* @brief Define a polar PDE solver for a Poisson-like equation.
 *
 * Solve the following Partial Differential Equation
 *
 * (1) @f$  L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho @f$, in  @f$ \Omega@f$,
 *
 * @f$  \phi = 0 @f$, on  @f$ \partial \Omega@f$,
 * 
 * As finite element basis functions we will use polar b-splines which are divided into two types:
 * 1) Basis splines that can be written as a tensor product of 1d basis splines
 *    ("non-singular B-splines")
 * 2) Basis splines that cover the centre point and are defined as a linear combination
 *    of basis splines of type 1 ("singular B-splines")
 * 
 * (see in Emily Bourne's thesis "Non-Uniform Numerical Schemes for the Modelling of Turbulence
 * in the 5D GYSELA Code". December 2022.)
 *
 * @tparam GridR The radial grid type.
 * @tparam GridR The poloidal grid type.
 * @tparam PolarBSplinesRTheta The type of the 2D polar B-splines (on the coordinate
 *          system @f$(r,\theta)@f$ including B-splines which traverse the O point).
 * @tparam SplineRThetaEvaluatorNullBound The type of the 2D (cross-product) spline evaluator.
 * @tparam IdxRangeFull The full index range of @f$ \phi @f$ including any batch dimensions.
 */
template <
        class GridR,
        class GridTheta,
        class PolarBSplinesRTheta,
        class SplineRThetaEvaluatorNullBound,
        class IdxRangeFull = IdxRange<GridR, GridTheta>>
class PolarSplineFEMPoissonLikeSolver
{
    // TODO: Add a batch loop to operator()
    static_assert(
            std::is_same_v<IdxRangeFull, IdxRange<GridR, GridTheta>>,
            "PolarSplineFEMPoissonLikeSolver is not yet batched");

public:
    /// The radial dimension
    using R = typename GridR::continuous_dimension_type;
    /// The poloidal dimension
    using Theta = typename GridTheta::continuous_dimension_type;

    static_assert(R::IS_CONTRAVARIANT);
    static_assert(Theta::IS_CONTRAVARIANT);

private:
    /// The radial dimension
    using R_cov = typename R::Dual;
    /// The poloidal dimension
    using Theta_cov = typename Theta::Dual;


public:
    /**
     * @brief Tag the first dimension for the quadrature mesh.
     */
    struct QDimRMesh : NonUniformGridBase<R>
    {
    };
    /**
     * @brief Tag the second dimension for the quadrature mesh.
     */
    struct QDimThetaMesh : NonUniformGridBase<Theta>
    {
    };

    /// The tag for the batch dimension for the equation. This is public due to Cuda.
    struct InternalBatchDim
    {
    };

private:
    using CoordRTheta = Coord<R, Theta>;
    /// The 1D B-splines in the radial direction
    using BSplinesR = typename PolarBSplinesRTheta::BSplinesR_tag;
    /// The 1D B-splines in the poloidal direction
    using BSplinesTheta = typename PolarBSplinesRTheta::BSplinesTheta_tag;

    using KnotsR = ddc::knot_discrete_dimension_t<BSplinesR>;
    using KnotsTheta = ddc::knot_discrete_dimension_t<BSplinesTheta>;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;

    /// The type of an index range over the polar B-splines
    using IdxRangeBSPolar = IdxRange<PolarBSplinesRTheta>;
    using IdxBSPolar = Idx<PolarBSplinesRTheta>;
    using IdxStepBSPolar = IdxStep<PolarBSplinesRTheta>;

    using IdxRangeBSR = IdxRange<BSplinesR>;
    using IdxRangeBSTheta = IdxRange<BSplinesTheta>;
    using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;

    using IdxBSR = Idx<BSplinesR>;
    using IdxBSTheta = Idx<BSplinesTheta>;
    using IdxBSRTheta = Idx<BSplinesR, BSplinesTheta>;

    using IdxStepBSR = IdxStep<BSplinesR>;
    using IdxStepBSTheta = IdxStep<BSplinesTheta>;
    using IdxStepBSRTheta = IdxStep<BSplinesR, BSplinesTheta>;

    using IdxRangeBatchedBSRTheta
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<IdxRangeFull>,
                    ddc::detail::TypeSeq<GridR, GridTheta>,
                    ddc::detail::TypeSeq<BSplinesR, BSplinesTheta>>>;

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFull, IdxRange<GridR>, IdxRange<GridTheta>>;

    /**
     * @brief Tag the quadrature index range in the first dimension.
     */
    using IdxRangeQuadratureR = IdxRange<QDimRMesh>;
    /**
     * @brief Tag the quadrature index range in the second dimension.
     */
    using IdxRangeQuadratureTheta = IdxRange<QDimThetaMesh>;
    /**
     * @brief Tag the quadrature index range.
     */
    using IdxRangeQuadratureRTheta = IdxRange<QDimRMesh, QDimThetaMesh>;
    /**
     * @brief Tag the elements (index) of the quadrature index range in the first dimension.
     */
    using IdxQuadratureR = Idx<QDimRMesh>;
    /**
     * @brief Tag the elements (index) of the quadrature index range in the second dimension.
     */
    using IdxQuadratureTheta = Idx<QDimThetaMesh>;
    /**
     * @brief Tag the elements (index) of the quadrature index range.
     */
    using IdxQuadratureRTheta = Idx<QDimRMesh, QDimThetaMesh>;
    /**
     * @brief Tag a vector on the first dimension of the quadrature mesh.
     */
    using IdxStepQuadratureR = IdxStep<QDimRMesh>;
    /**
     * @brief Tag a vector on the second dimension of the quadrature mesh.
     */
    using IdxStepQuadratureTheta = IdxStep<QDimThetaMesh>;

    using ConstSpline2D = DConstField<IdxRangeBatchedBSRTheta>;
    using PolarSplineMemRTheta = DFieldMem<IdxRange<PolarBSplinesRTheta>>;
    using PolarSplineRTheta = DField<IdxRange<PolarBSplinesRTheta>>;

    using CoordFieldMemRTheta = FieldMem<CoordRTheta, IdxRangeRTheta>;
    using CoordFieldRTheta = Field<CoordRTheta, IdxRangeRTheta>;
    using DFieldRTheta = DField<IdxRangeRTheta>;

public:
    // TODO: Are these types needed?
    /**
    * @brief Object storing a value and a value of the derivative
    * of a 1D function.
    */
    struct EvalDeriv1DType
    {
        /// The value of the function @f$f(x)@f$.
        double value;
        /// The derivative of the function @f$\partial_x f(x)@f$.
        double derivative;
    };

    /**
    * @brief Object storing a value and a value of the derivatives
    * in each direction of a 2D function.
    */
    struct EvalDeriv2DType
    {
        /// The value of the function @f$f(r, \theta)@f$.
        double value;
        /// The gradient of the function @f$\nabla f(r, \theta)@f$.
        DVector<R_cov, Theta_cov> derivative;
    };

private:
    static constexpr int s_n_gauss_legendre_r = BSplinesR::degree() + 1;
    static constexpr int s_n_gauss_legendre_theta = BSplinesTheta::degree() + 1;
    // The number of cells (in the radial direction) in which both types of basis splines can be found
    static constexpr int m_n_overlap_cells = PolarBSplinesRTheta::continuity + 1;

    // Number of cells over which a radial B-splines has its support
    // This is the case for b-splines which are not affected by the higher knot multiplicity at the boundary.
    static constexpr int m_n_non_zero_bases_r = BSplinesR::degree() + 1;

    // Number of cells over which a poloidal B-splines has its support
    static constexpr int m_n_non_zero_bases_theta = BSplinesTheta::degree() + 1;

    const int m_nbasis_theta;

    // Domains
    IdxRangeBSPolar m_idxrange_fem_non_singular;
    IdxRangeBSR m_idxrange_bsplines_r;
    IdxRangeBSTheta m_idxrange_bsplines_theta;

    IdxRangeQuadratureR m_idxrange_quadrature_r;
    IdxRangeQuadratureTheta m_idxrange_quadrature_theta;
    IdxRangeQuadratureRTheta m_idxrange_quadrature_singular;
    IdxRangeQuadratureRTheta m_idxrange_quadrature;

    FieldMem<double, IdxRangeQuadratureRTheta> m_int_volume_alloc;

    PolarSplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            PolarBSplinesRTheta,
            ddc::NullExtrapolationRule>
            m_polar_spline_evaluator;
    std::unique_ptr<MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>
            m_gko_matrix;
    mutable PolarSplineMemRTheta m_phi_spline_coef_alloc;
    mutable DFieldMem<IdxRange<InternalBatchDim, PolarBSplinesRTheta>> m_x_init_alloc;

    const int m_batch_idx {0}; // TODO: Remove when batching is supported
public:
    /**
     * @brief Instantiate a polar Poisson-like solver using FEM with B-splines.
     *
     * The equation we are studying:
     *
     * (1) @f$  L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho @f$, in  @f$ \Omega@f$,
     *
     *  @f$  \phi = 0 @f$, on  @f$ \partial \Omega@f$.
     *
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @param[in] spline_evaluator
     *      An evaluator for evaluating 2D splines on @f$(r,\theta)@f$.
     * @param[in] max_iter
     *      The maximum number of iterations possible for the batched CSR solver.
     * @param[in] res_tol
     *      The residual tolerance for the batched CSR solver. Be careful! the relative residual
     *      provided here, will be used as "implicit residual" in ginkgo solver.
     * @param[in] batch_solver_logger
     *      Indicates whether log information such as the residual and the number of iterations
     *      should be monitored.
     * @param[in] preconditioner_max_block_size
     *      The maximum size of the Jacobi preconditioner used by the batched CSR solver.
     *
     * @tparam Mapping A class describing a mapping from curvilinear coordinates to Cartesian coordinates.
     */
    template <class Mapping>
    PolarSplineFEMPoissonLikeSolver(
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
            std::optional<int> max_iter = std::nullopt,
            std::optional<double> res_tol = std::nullopt,
            std::optional<bool> batch_solver_logger = std::nullopt,
            std::optional<int> preconditioner_max_block_size = std::nullopt)
        : m_nbasis_theta(ddc::discrete_space<BSplinesTheta>().nbasis())
        , m_idxrange_fem_non_singular(
                  ddc::discrete_space<PolarBSplinesRTheta>().tensor_bspline_idx_range().remove_last(
                          IdxStepBSPolar {m_nbasis_theta}))
        , m_idxrange_bsplines_r(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                  IdxStepBSR {m_n_overlap_cells}))
        , m_idxrange_bsplines_theta(ddc::discrete_space<BSplinesTheta>().full_domain().take_first(
                  IdxStepBSTheta {m_nbasis_theta}))
        , m_idxrange_quadrature_r(
                  Idx<QDimRMesh>(0),
                  IdxStep<QDimRMesh>(
                          s_n_gauss_legendre_r * ddc::discrete_space<BSplinesR>().ncells()))
        , m_idxrange_quadrature_theta(
                  Idx<QDimThetaMesh>(0),
                  IdxStep<QDimThetaMesh>(
                          s_n_gauss_legendre_theta * ddc::discrete_space<BSplinesTheta>().ncells()))
        , m_idxrange_quadrature_singular(
                  m_idxrange_quadrature_r.take_first(
                          IdxStep<QDimRMesh> {m_n_overlap_cells * s_n_gauss_legendre_r}),
                  m_idxrange_quadrature_theta)
        , m_idxrange_quadrature(m_idxrange_quadrature_r, m_idxrange_quadrature_theta)
        , m_polar_spline_evaluator(ddc::NullExtrapolationRule())
        , m_phi_spline_coef_alloc(ddc::discrete_space<PolarBSplinesRTheta>().full_domain())
        , m_x_init_alloc(
                  "x_init",
                  IdxRange<InternalBatchDim, PolarBSplinesRTheta>(
                          Idx<InternalBatchDim, PolarBSplinesRTheta>(0, 0),
                          IdxStep<InternalBatchDim, PolarBSplinesRTheta>(
                                  1,
                                  ddc::discrete_space<PolarBSplinesRTheta>().nbasis()
                                          - ddc::discrete_space<BSplinesTheta>().nbasis())))
    {
        static_assert(has_jacobian_v<Mapping>);
        //initialise x_init
        ddc::parallel_fill(get_field(m_x_init_alloc), 0);
        // Get break points
        IdxRange<KnotsR> idxrange_r_edges = ddc::discrete_space<BSplinesR>().break_point_domain();
        IdxRange<KnotsTheta> idxrange_theta_edges
                = ddc::discrete_space<BSplinesTheta>().break_point_domain();
        host_t<FieldMem<Coord<R>, IdxRange<KnotsR>>> breaks_r(idxrange_r_edges);
        host_t<FieldMem<Coord<Theta>, IdxRange<KnotsTheta>>> breaks_theta(idxrange_theta_edges);

        ddcHelper::dump_coordinates(Kokkos::DefaultHostExecutionSpace(), get_field(breaks_r));
        ddcHelper::dump_coordinates(Kokkos::DefaultHostExecutionSpace(), get_field(breaks_theta));

        // Define quadrature points and weights
        GaussLegendre<QDimRMesh, s_n_gauss_legendre_r> gl_coeffs_r(get_const_field(breaks_r));
        GaussLegendre<QDimThetaMesh, s_n_gauss_legendre_theta> gl_coeffs_theta(
                get_const_field(breaks_theta));
        m_int_volume_alloc = compute_coeffs_on_mapping(
                Kokkos::DefaultExecutionSpace(),
                mapping,
                gauss_legendre_quadrature_coefficients<
                        Kokkos::DefaultExecutionSpace>(gl_coeffs_r, gl_coeffs_theta));

        using PoissonAssembler = PolarSplineFEMPoissonLikeAssembler<
                GridR,
                GridTheta,
                PolarBSplinesRTheta,
                SplineRThetaEvaluatorNullBound,
                QDimRMesh,
                QDimThetaMesh>;
        PoissonAssembler assembler(get_field(m_int_volume_alloc));
        assembler(
                m_gko_matrix,
                coeff_alpha,
                coeff_beta,
                mapping,
                spline_evaluator,
                max_iter,
                res_tol,
                batch_solver_logger,
                preconditioner_max_block_size);
    }

    /**
     * @brief Solve the Poisson-like equation.
     *
     * This operator returns the coefficients associated with the B-Splines
     * of the solution @f$\phi@f$.
     *
     * @param[in] rhs
     *      The rhs @f$ \rho@f$ of the Poisson-like equation.
     *      The type is templated but we can use the PoissonLikeRHSFunction
     *      class. It must be an object with an operator() which evaluates a
     *      CoordRTheta and can be called from GPU.
     * @param[inout] spline
     *      The spline representation of the solution @f$\phi@f$.
     */
    template <class RHSFunction>
    void operator()(RHSFunction const& rhs, PolarSplineRTheta spline) const
    {
        Kokkos::Profiling::pushRegion("PolarPoissonRHS");

        static_assert(
                std::is_invocable_r_v<double, RHSFunction, CoordRTheta>,
                "RHSFunction must have an operator() which takes a coordinate and returns a "
                "double");
        assert(get_idx_range(spline) == ddc::discrete_space<PolarBSplinesRTheta>().full_domain());
        IdxRange<InternalBatchDim> batch_idx_range(get_idx_range(m_x_init_alloc));

        assert(batch_idx_range.size() == 1);

        Idx<InternalBatchDim> batch_idx = batch_idx_range.front();

        // Create b for rhs
        DFieldMem<IdxRange<InternalBatchDim, PolarBSplinesRTheta>> b_alloc(
                get_idx_range(m_x_init_alloc));
        DField<IdxRange<InternalBatchDim, PolarBSplinesRTheta>> b = get_field(b_alloc);

        // Get initial guess
        DField<IdxRange<InternalBatchDim, PolarBSplinesRTheta>> x_init = get_field(m_x_init_alloc);

        DConstField<IdxRangeQuadratureRTheta> int_volume = get_const_field(m_int_volume_alloc);

        IdxRangeBSPolar idx_range_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();

        IdxRangeQuadratureRTheta idx_range_quad_singular = m_idxrange_quadrature_singular;

        // Fill b
        // Multi-level parallelism is needed as idx_range_singular.size() ~= 3 but b is on GPU
        Kokkos::parallel_for(
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        idx_range_singular.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    IdxBSPolar idx
                            = idx_range_singular.front() + IdxStepBSPolar(team.league_rank());
                    double teamSum = 0;
                    Kokkos::parallel_reduce(
                            Kokkos::TeamThreadMDRange(
                                    team,
                                    idx_range_quad_singular.template extent<QDimRMesh>(),
                                    idx_range_quad_singular.template extent<QDimThetaMesh>()),
                            [&](int r_thread_index, int theta_thread_index, double& sum) {
                                IdxQuadratureRTheta idx_quad = idx_range_quad_singular.front()
                                                               + IdxStep<QDimRMesh, QDimThetaMesh>(
                                                                       r_thread_index,
                                                                       theta_thread_index);
                                const CoordRTheta coord(ddc::coordinate(idx_quad));
                                sum += rhs(coord) * get_polar_bspline_vals(coord, idx)
                                       * int_volume(idx_quad);
                            },
                            teamSum);

                    b(batch_idx, idx) = teamSum;
                });

        IdxRangeQuadratureRTheta full_quad_idx_range = m_idxrange_quadrature;
        IdxRangeQuadratureTheta full_quad_idx_range_theta(full_quad_idx_range);

        const std::source_location location = std::source_location::current();ddc::parallel_for_each(location.function_name(),
                m_idxrange_fem_non_singular,
                KOKKOS_LAMBDA(IdxBSPolar const idx) {
                    const IdxBSRTheta idx_2d(PolarBSplinesRTheta::get_2d_index(idx));
                    const IdxBSR idx_r(idx_2d);
                    const IdxBSTheta idx_theta(idx_2d);

                    auto& bspl_r = ddc::discrete_space<BSplinesR>();
                    auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

                    // Find the cells on which the bspline is non-zero
                    const Idx<KnotsR> start_non_zero_r(
                            Kokkos::
                                    max(bspl_r.break_point_domain().front(),
                                        bspl_r.get_first_support_knot(idx_r)));
                    const Idx<KnotsR> end_non_zero_r(
                            Kokkos::
                                    min(bspl_r.break_point_domain().back(),
                                        bspl_r.get_last_support_knot(idx_r)));

                    const Idx<KnotsTheta> start_non_zero_theta(
                            bspl_theta.get_first_support_knot(idx_theta));
                    const Idx<KnotsTheta> end_non_zero_theta(
                            bspl_theta.get_last_support_knot(idx_theta));
                    const IdxRangeQuadratureRTheta quad_range
                            = detail_poisson::get_quadrature_between_knots<
                                    QDimRMesh,
                                    QDimThetaMesh,
                                    BSplinesR,
                                    BSplinesTheta>(
                                    start_non_zero_r,
                                    end_non_zero_r,
                                    start_non_zero_theta,
                                    end_non_zero_theta,
                                    full_quad_idx_range.front());

                    // Calculate the weak integral
                    b(batch_idx, idx) = 0.0;
                    for (IdxQuadratureTheta idx_quad_theta :
                         ddc::select<QDimThetaMesh>(quad_range)) {
                        // Manage periodicity
                        if (!full_quad_idx_range_theta.contains(idx_quad_theta)) {
                            idx_quad_theta -= full_quad_idx_range_theta.extents();
                        }
                        for (IdxQuadratureR idx_quad_r : ddc::select<QDimRMesh>(quad_range)) {
                            IdxQuadratureRTheta idx_quad(idx_quad_r, idx_quad_theta);
                            CoordRTheta coord(ddc::coordinate(idx_quad));
                            b(batch_idx, idx) += rhs(coord) * get_polar_bspline_vals(coord, idx)
                                                 * int_volume(idx_quad);
                        }
                    }
                });

        Kokkos::Profiling::popRegion();

        // Solve the matrix equation
        Kokkos::Profiling::pushRegion("PolarPoissonSolve");
        m_gko_matrix->solve(x_init.allocation_kokkos_view(), b.allocation_kokkos_view());
        //-----------------
        IdxStepBSPolar radial_boundary_splines(m_nbasis_theta);
        IdxRangeBSPolar polar_bspl_idx_range
                = ddc::discrete_space<PolarBSplinesRTheta>().full_domain().remove_last(
                        radial_boundary_splines);
        IdxRangeBSPolar bc_polar_bspl_idx_range
                = ddc::discrete_space<PolarBSplinesRTheta>().full_domain().take_last(
                        radial_boundary_splines);

        // Fill the spline
        ddc::parallel_for_each(location.function_name(),
                polar_bspl_idx_range,
                KOKKOS_LAMBDA(IdxBSPolar const idx) { spline(idx) = x_init(batch_idx, idx); });
        ddc::parallel_fill(spline[bc_polar_bspl_idx_range], 0.0);
        Kokkos::Profiling::popRegion();
    }

    /**
     * @brief Solve the Poisson-like equation.
     *
     * This operator uses the other operator () and returns the values on
     * the grid of the solution @f$\phi@f$.
     *
     * @param[in] rhs
     *      The rhs @f$ \rho@f$ of the Poisson-like equation.
     *      The type is templated but we can use the PoissonLikeRHSFunction
     *      class. It must be an object with an operator() which evaluates a
     *      CoordRTheta and can be called from GPU.
     * @param[inout] phi
     *      The values of the solution @f$\phi@f$ on the given coords_eval.
     */
    template <class RHSFunction>
    void operator()(RHSFunction const& rhs, DFieldRTheta phi) const
    {
        static_assert(
                std::is_invocable_r_v<double, RHSFunction, CoordRTheta>,
                "RHSFunction must have an operator() which takes a coordinate and returns a "
                "double");

        (*this)(rhs, get_field(m_phi_spline_coef_alloc));
        CoordFieldMemRTheta coords_eval_alloc(get_idx_range(phi));
        m_polar_spline_evaluator(phi, get_const_field(m_phi_spline_coef_alloc));
    }

    /**
     * @brief Get the value of the specified polar bspline at the specified point.
     *
     * @param[in] coord
     *      The coordinate where the value of the polar bspline should be calculated.
     * @param[in] idx
     *      The polar bspline of interest.
     * @return The value of the polar bspline at the coordinate.
     */
    static KOKKOS_INLINE_FUNCTION double get_polar_bspline_vals(CoordRTheta coord, IdxBSPolar idx)
    {
        double val;
        detail_poisson::
                get_polar_bspline_vals_and_derivs<PolarBSplinesRTheta, false>(val, coord, idx);
        return val;
    }
};
