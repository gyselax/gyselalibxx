// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "gauss_legendre_integration.hpp"
#include "math_tools.hpp"
#include "matrix_batch_csr.hpp"
#include "metric_tensor_evaluator.hpp"
#include "polar_spline_evaluator.hpp"
#include "quadrature_coeffs_nd.hpp"
#include "view.hpp"
#include "volume_quadrature_nd.hpp"


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

    class InternalBatchDim;

private:
    using CoordRTheta = Coord<R, Theta>;
    /// The 1D B-splines in the radial direction
    using BSplinesR = typename PolarBSplinesRTheta::BSplinesR_tag;
    /// The 1D B-splines in the poloidal direction
    using BSplinesTheta = typename PolarBSplinesRTheta::BSplinesTheta_tag;

    using KnotsR = ddc::knot_discrete_dimension_t<BSplinesR>;
    using KnotsTheta = ddc::knot_discrete_dimension_t<BSplinesTheta>;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRTheta = Idx<GridR, GridTheta>;

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

    const int m_nbasis_r;
    const int m_nbasis_theta;

    // Matrix size is equal to the number of Polar bspline
    const int m_matrix_size;

    // Domains
    IdxRangeBSPolar m_idxrange_fem_non_singular;
    IdxRangeBSR m_idxrange_bsplines_r;
    IdxRangeBSTheta m_idxrange_bsplines_theta;

    IdxRangeQuadratureR m_idxrange_quadrature_r;
    IdxRangeQuadratureTheta m_idxrange_quadrature_theta;
    IdxRangeQuadratureRTheta m_idxrange_quadrature;
    IdxRangeQuadratureRTheta m_idxrange_quadrature_singular;

    // Gauss-Legendre points and weights
    FieldMem<double, IdxRangeQuadratureR> m_weights_r;
    FieldMem<double, IdxRangeQuadratureTheta> m_weights_theta;

    FieldMem<double, IdxRangeQuadratureRTheta> m_int_volume_alloc;

    PolarSplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            PolarBSplinesRTheta,
            ddc::NullExtrapolationRule>
            m_polar_spline_evaluator;
    std::unique_ptr<MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>
            m_gko_matrix;
    mutable PolarSplineMemRTheta m_phi_spline_coef;
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
        : m_nbasis_r(ddc::discrete_space<BSplinesR>().nbasis() - m_n_overlap_cells - 1)
        , m_nbasis_theta(ddc::discrete_space<BSplinesTheta>().nbasis())
        , m_matrix_size(ddc::discrete_space<PolarBSplinesRTheta>().nbasis() - m_nbasis_theta)
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
        , m_weights_r(m_idxrange_quadrature_r)
        , m_weights_theta(m_idxrange_quadrature_theta)
        , m_polar_spline_evaluator(ddc::NullExtrapolationRule())
        , m_phi_spline_coef(ddc::discrete_space<PolarBSplinesRTheta>().full_domain())
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

        // Number of elements in the matrix that correspond to the splines
        // that cover the singular point
        constexpr int n_elements_singular
                = PolarBSplinesRTheta::n_singular_basis() * PolarBSplinesRTheta::n_singular_basis();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // polar splines at the singular point and the other splines
        const int n_elements_overlap = 2
                                       * (PolarBSplinesRTheta::n_singular_basis()
                                          * BSplinesR::degree() * m_nbasis_theta);
        const int n_stencil_theta
                = m_nbasis_theta * min(int(1 + 2 * BSplinesTheta::degree()), m_nbasis_theta);
        const int n_stencil_r = m_nbasis_r * (1 + 2 * BSplinesR::degree())
                                - (1 + BSplinesR::degree()) * BSplinesR::degree();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // non-central splines. These have a tensor product structure
        const int n_elements_stencil = n_stencil_r * n_stencil_theta;

        const int batch_size = 1;

        const int n_matrix_elements = n_elements_singular + n_elements_overlap + n_elements_stencil;

        //CSR data storage
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>
                values_csr_host("values_csr", batch_size, n_matrix_elements);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace>
                col_idx_csr_host("idx_csr", n_matrix_elements);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
                nnz_per_row_csr("nnz_per_row_csr", m_matrix_size + 1);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace>
                nnz_per_row_csr_host("nnz_per_row_csr", m_matrix_size + 1);

        m_gko_matrix = std::make_unique<
                MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>(
                1,
                m_matrix_size,
                n_matrix_elements,
                max_iter,
                res_tol,
                batch_solver_logger,
                preconditioner_max_block_size);
        auto [values, col_idx, nnz_per_row] = m_gko_matrix->get_batch_csr();
        init_nnz_per_line(nnz_per_row);
        Kokkos::deep_copy(nnz_per_row_csr_host, nnz_per_row);

        compute_singular_elements(
                coeff_alpha,
                coeff_beta,
                mapping,
                spline_evaluator,
                values_csr_host,
                col_idx_csr_host,
                nnz_per_row_csr_host);
        compute_overlapping_singular_elements(
                coeff_alpha,
                coeff_beta,
                mapping,
                spline_evaluator,
                values_csr_host,
                col_idx_csr_host,
                nnz_per_row_csr_host);
        compute_stencil_elements(
                coeff_alpha,
                coeff_beta,
                mapping,
                spline_evaluator,
                values_csr_host,
                col_idx_csr_host,
                nnz_per_row_csr_host);

        assert(nnz_per_row_csr_host(m_matrix_size) == n_matrix_elements);
        Kokkos::deep_copy(values, values_csr_host);
        Kokkos::deep_copy(col_idx, col_idx_csr_host);
        Kokkos::deep_copy(nnz_per_row, nnz_per_row_csr_host);
        m_gko_matrix->setup_solver();
    }

    /**
     * @brief Computes the matrix element corresponding to the singular area.
     *        ie: the region enclosing the O-point.
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
     * @param[out] values_csr_host
     *             A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[out] col_idx_csr_host
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix).
     * @param[inout] nnz_per_row_csr_host
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    template <class Mapping>
    void compute_singular_elements(
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
            Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> const values_csr_host,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> const col_idx_csr_host,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> const nnz_per_row_csr_host)
    {
        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();
        IdxRangeQuadratureRTheta idx_range_quad_singular = m_idxrange_quadrature_singular;

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = get_field(m_int_volume_alloc);

        Kokkos::Profiling::pushRegion("PolarPoissonFillFemMatrix");
        // Calculate the matrix elements corresponding to the B-splines which cover the singular point
        ddc::for_each(idxrange_singular, [&](IdxBSPolar const idx_test) {
            ddc::for_each(idxrange_singular, [&](IdxBSPolar const idx_trial) {
                // Calculate the weak integral
                double const element = ddc::parallel_transform_reduce(
                        Kokkos::DefaultExecutionSpace(),
                        idx_range_quad_singular,
                        0.0,
                        ddc::reducer::sum<double>(),
                        KOKKOS_LAMBDA(Idx<QDimRMesh, QDimThetaMesh> const& idx_quad) {
                            return weak_integral_element(
                                    idx_test,
                                    idx_trial,
                                    idx_quad,
                                    coeff_alpha,
                                    coeff_beta,
                                    spline_evaluator,
                                    mapping,
                                    int_volume_proxy);
                        });
                const int row_idx = idx_test - idxrange_singular.front();
                const int col_idx = idx_trial - idxrange_singular.front();
                const int csr_idx_singular_area = nnz_per_row_csr_host(row_idx + 1);
                //Fill the C matrix corresponding to the splines on the singular point
                col_idx_csr_host(csr_idx_singular_area) = col_idx;
                values_csr_host(m_batch_idx, csr_idx_singular_area) = element;
                nnz_per_row_csr_host(row_idx + 1)++;
            });
        });
    }

    /**
     * @brief Computes the matrix element corresponding to singular elements overlapping
     *        with regular grid.
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
     * @param[out] values_csr_host
     *             A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[out] col_idx_csr_host
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @param[inout] nnz_per_row_csr_host
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    template <class Mapping>
    void compute_overlapping_singular_elements(
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
            Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> const values_csr_host,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> const col_idx_csr_host,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> const nnz_per_row_csr_host)
    {
        // Create index ranges associated with the 2D splines
        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();
        IdxRangeBSR central_radial_bspline_idx_range(
                m_idxrange_bsplines_r.take_first(IdxStep<BSplinesR> {BSplinesR::degree()}));

        IdxRangeBSRTheta idxrange_non_singular_near_centre(
                central_radial_bspline_idx_range,
                m_idxrange_bsplines_theta);

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = get_field(m_int_volume_alloc);
        IdxRangeQuadratureRTheta
                full_quad_idx_range(m_idxrange_quadrature_r, m_idxrange_quadrature_theta);

        // Calculate the matrix elements where bspline products overlap the B-splines which cover the singular point
        ddc::for_each(idxrange_singular, [&](IdxBSPolar const idx_test) {
            ddc::for_each(idxrange_non_singular_near_centre, [&](IdxBSRTheta const idx_trial) {
                const IdxBSPolar idx_trial_polar(to_polar(idx_trial));
                const IdxBSR idx_trial_r(idx_trial);
                const IdxBSTheta idx_trial_theta(idx_trial);

                auto& bspl_r = ddc::discrete_space<BSplinesR>();
                auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

                // Find the index range covering the cells where both the test and trial functions are non-zero
                const Idx<KnotsR> start_non_zero_r(
                        ::max(bspl_r.break_point_domain().front(),
                              bspl_r.get_first_support_knot(idx_trial_r)));
                const Idx<KnotsR> end_non_zero_r(
                        ::min(bspl_r.get_last_support_knot(IdxBSR(PolarBSplinesRTheta::continuity)),
                              bspl_r.get_last_support_knot(idx_trial_r)));

                const Idx<KnotsTheta> start_non_zero_theta(
                        bspl_theta.get_first_support_knot(idx_trial_theta));
                const Idx<KnotsTheta> end_non_zero_theta(
                        bspl_theta.get_last_support_knot(idx_trial_theta));

                const IdxRangeQuadratureRTheta quad_range = get_quadrature_between_knots(
                        start_non_zero_r,
                        end_non_zero_r,
                        start_non_zero_theta,
                        end_non_zero_theta,
                        m_idxrange_quadrature.front());


                assert(quad_range.size() > 0);
                // Calculate the weak integral
                double element = ddc::parallel_transform_reduce(
                        Kokkos::DefaultExecutionSpace(),
                        quad_range,
                        0.0,
                        ddc::reducer::sum<double>(),
                        KOKKOS_LAMBDA(IdxQuadratureRTheta idx_quad) {
                            // Manage periodicity
                            if (!full_quad_idx_range.contains(idx_quad)) {
                                idx_quad -= full_quad_idx_range.template extent<QDimThetaMesh>();
                            }

                            return weak_integral_element<Mapping>(
                                    idx_test,
                                    idx_trial_polar,
                                    idx_quad,
                                    coeff_alpha,
                                    coeff_beta,
                                    spline_evaluator,
                                    mapping,
                                    int_volume_proxy);
                        });

                int const row_idx = idx_test - idxrange_singular.front();
                int const col_idx = idx_trial_polar - idxrange_singular.front();
                //a_ij
                col_idx_csr_host(nnz_per_row_csr_host(row_idx + 1)) = col_idx;
                values_csr_host(m_batch_idx, nnz_per_row_csr_host(row_idx + 1)) = element;
                nnz_per_row_csr_host(row_idx + 1)++;
                //a_ji
                col_idx_csr_host(nnz_per_row_csr_host(col_idx + 1)) = row_idx;
                values_csr_host(m_batch_idx, nnz_per_row_csr_host(col_idx + 1)) = element;
                nnz_per_row_csr_host(col_idx + 1)++;
            });
        });
    }

    /**
     * @brief Computes the matrix element corresponding to the regular stencil
     *        ie: out to singular or overlapping areas.
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
     * @param[out] values_csr_host
     *             A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[out] col_idx_csr_host
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @param[inout] nnz_per_row_csr_host
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    template <class Mapping>
    void compute_stencil_elements(
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
            Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> const values_csr_host,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> const col_idx_csr_host,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> const nnz_per_row_csr_host)
    {
        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();

        // Get index range for basis elements (last element removed due to homogeneous Dirichlet)
        IdxRangeBSR full_idx_range_r
                = ddc::discrete_space<BSplinesR>().full_domain().remove_last(IdxStepBSR(1));

        // Calculate the matrix elements following a stencil
        ddc::for_each(m_idxrange_fem_non_singular, [&](IdxBSPolar const idx_test_polar) {
            const IdxBSRTheta idx_test(PolarBSplinesRTheta::get_2d_index(idx_test_polar));
            const IdxBSR idx_test_r(idx_test);
            const IdxBSTheta idx_test_theta(idx_test);

            // Calculate the index of the elements that are already filled
            IdxRangeBSTheta remaining_theta(
                    idx_test_theta,
                    IdxStep<BSplinesTheta> {BSplinesTheta::degree() + 1});
            ddc::for_each(remaining_theta, [&](IdxBSTheta const idx_trial_theta) {
                IdxBSRTheta idx_trial(idx_test_r, idx_trial_theta);
                IdxBSPolar idx_trial_polar(to_polar(theta_mod(idx_trial)));
                double element = get_matrix_stencil_element(
                        idx_test,
                        idx_trial,
                        coeff_alpha,
                        coeff_beta,
                        spline_evaluator,
                        mapping);
                int const int_polar_idx_test = idx_test_polar - idxrange_singular.front();
                if (idx_test_polar == idx_trial_polar) {
                    const int idx = nnz_per_row_csr_host(int_polar_idx_test + 1);
                    col_idx_csr_host(idx) = int_polar_idx_test;
                    values_csr_host(m_batch_idx, idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_test + 1)++;
                } else {
                    int const int_polar_idx_trial = idx_trial_polar - idxrange_singular.front();

                    const int aij_idx = nnz_per_row_csr_host(int_polar_idx_test + 1);
                    col_idx_csr_host(aij_idx) = int_polar_idx_trial;
                    values_csr_host(m_batch_idx, aij_idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_test + 1)++;

                    const int aji_idx = nnz_per_row_csr_host(int_polar_idx_trial + 1);
                    col_idx_csr_host(aji_idx) = int_polar_idx_test;
                    values_csr_host(m_batch_idx, aji_idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_trial + 1)++;
                }
            });
            IdxStepBSR n_remaining_r(
                    ::min(IdxStepBSR(BSplinesR::degree()), full_idx_range_r.back() - idx_test_r));
            IdxRangeBSR remaining_r(idx_test_r + 1, n_remaining_r);
            IdxRangeBSTheta relevant_theta(
                    idx_test_theta + ddc::discrete_space<BSplinesTheta>().nbasis()
                            - BSplinesTheta::degree(),
                    IdxStepBSTheta {2 * BSplinesTheta::degree() + 1});

            IdxRangeBSRTheta trial_idx_range(remaining_r, relevant_theta);

            ddc::for_each(trial_idx_range, [&](IdxBSRTheta const idx_trial) {
                IdxBSPolar idx_trial_polar(to_polar(theta_mod(idx_trial)));
                double element = get_matrix_stencil_element(
                        idx_test,
                        idx_trial,
                        coeff_alpha,
                        coeff_beta,
                        spline_evaluator,
                        mapping);
                int const int_polar_idx_test = idx_test_polar - idxrange_singular.front();
                if (idx_test_polar == idx_trial_polar) {
                    const int idx = nnz_per_row_csr_host(int_polar_idx_test + 1);
                    col_idx_csr_host(idx) = int_polar_idx_test;
                    values_csr_host(m_batch_idx, idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_test + 1)++;
                } else {
                    int const int_polar_idx_trial = idx_trial_polar - idxrange_singular.front();
                    const int aij_idx = nnz_per_row_csr_host(int_polar_idx_test + 1);
                    col_idx_csr_host(aij_idx) = int_polar_idx_trial;
                    values_csr_host(m_batch_idx, aij_idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_test + 1)++;

                    const int aji_idx = nnz_per_row_csr_host(int_polar_idx_trial + 1);
                    col_idx_csr_host(aji_idx) = int_polar_idx_test;
                    values_csr_host(m_batch_idx, aji_idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_trial + 1)++;
                }
            });
        });

        Kokkos::Profiling::popRegion();
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
     *      CoordRTheta and can be called from CPU.
     * @param[inout] spline
     *      The spline representation of the solution @f$\phi@f$, also used as initial data for the iterative solver.
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

        DConstField<IdxRangeQuadratureRTheta> int_volume = get_field(m_int_volume_alloc);

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

        ddc::parallel_for_each(
                m_idxrange_fem_non_singular,
                KOKKOS_LAMBDA(IdxBSPolar const idx) {
                    const IdxBSRTheta idx_2d(PolarBSplinesRTheta::get_2d_index(idx));
                    const IdxBSR idx_r(idx_2d);
                    const IdxBSTheta idx_theta(idx_2d);

                    auto& bspl_r = ddc::discrete_space<BSplinesR>();
                    auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

                    // Find the cells on which the bspline is non-zero
                    const Idx<KnotsR> start_non_zero_r(
                            ::max(bspl_r.break_point_domain().front(),
                                  bspl_r.get_first_support_knot(idx_r)));
                    const Idx<KnotsR> end_non_zero_r(
                            ::min(bspl_r.break_point_domain().back(),
                                  bspl_r.get_last_support_knot(idx_r)));

                    const Idx<KnotsTheta> start_non_zero_theta(
                            bspl_theta.get_first_support_knot(idx_theta));
                    const Idx<KnotsTheta> end_non_zero_theta(
                            bspl_theta.get_last_support_knot(idx_theta));

                    const IdxRangeQuadratureRTheta quad_range = get_quadrature_between_knots(
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
        ddc::parallel_for_each(
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
     *      CoordRTheta and can be called from CPU.
     * @param[inout] phi
     *      The values of the solution @f$\phi@f$ on the given coords_eval, also used as initial data for the iterative solver.
     */
    template <class RHSFunction>
    void operator()(RHSFunction const& rhs, DFieldRTheta phi) const
    {
        static_assert(
                std::is_invocable_r_v<double, RHSFunction, CoordRTheta>,
                "RHSFunction must have an operator() which takes a coordinate and returns a "
                "double");

        (*this)(rhs, get_field(m_phi_spline_coef));
        CoordFieldMemRTheta coords_eval_alloc(get_idx_range(phi));
        m_polar_spline_evaluator(phi, get_const_field(m_phi_spline_coef));
    }

    /**
     * @brief Computes a quadrature summand corresponding to the 
     *        inner product.
     *
     * @param[in] idx_test
     *      The index of the test basis spline.
     * @param[in] idx_trial
     *      The index of the trial basis spline.
     * @param[in] idx_quad
     *      The index for the point in the quadrature scheme.
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] spline_evaluator
     *      An evaluator for evaluating 2D splines on @f$(r,\theta)@f$.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @param[in] int_volume
     *      The integral volume associated with each point used in the quadrature scheme.
     * @return
      inner product of the test and trial spline is computed using a 
     * quadrature. This function returns one summand of the quadrature for 
     * the quadrature point given by the indices.
     */
    template <class Mapping>
    static KOKKOS_FUNCTION double weak_integral_element(
            IdxBSPolar idx_test,
            IdxBSPolar idx_trial,
            IdxQuadratureRTheta idx_quad,
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
            Mapping const& mapping,
            DField<IdxRangeQuadratureRTheta> int_volume)
    {
        // Calculate coefficients at quadrature point
        CoordRTheta coord(ddc::coordinate(idx_quad));
        const double alpha = spline_evaluator(coord, coeff_alpha);
        const double beta = spline_evaluator(coord, coeff_beta);

        // Define the value and gradient of the test and trial basis functions
        double basis_val_test_space;
        double basis_val_trial_space;
        DVector<R_cov, Theta_cov> basis_derivs_test_space
                = get_polar_bspline_vals_and_derivs(basis_val_test_space, coord, idx_test);
        DVector<R_cov, Theta_cov> basis_derivs_trial_space
                = get_polar_bspline_vals_and_derivs(basis_val_trial_space, coord, idx_trial);

        MetricTensorEvaluator<Mapping, CoordRTheta> get_metric_tensor(mapping);

        Tensor inv_metric_tensor = get_metric_tensor.inverse(coord);

        // Assemble the weak integral element
        return int_volume(idx_quad)
               * (alpha
                          * tensor_mul(
                                  index<'i'>(basis_derivs_test_space),
                                  index<'i', 'j'>(inv_metric_tensor),
                                  index<'j'>(basis_derivs_trial_space))
                  + beta * basis_val_test_space * basis_val_trial_space);
    }

    /**
     * @brief Computes the matrix element corresponding to two tensor product splines
     *        with index idx_test and idx_trial
     *
     * @param[in] idx_test
     *      The index for polar B-spline in the test space.
     * @param[in] idx_trial
     *      The index for polar B-spline in the trial space.
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] evaluator
     *      An evaluator for evaluating 2D splines on @f$ (r, \theta) @f$.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @return 
     *      The value of the matrix element.
     */
    template <class Mapping>
    double get_matrix_stencil_element(
            IdxBSRTheta idx_test,
            IdxBSRTheta idx_trial,
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        IdxRangeQuadratureRTheta
                full_quad_idx_range(m_idxrange_quadrature_r, m_idxrange_quadrature_theta);
        const IdxBSR idx_test_r(idx_test);
        const IdxBSR idx_trial_r(idx_trial);
        const IdxBSTheta idx_test_theta(theta_mod(IdxBSTheta(idx_test)));
        const IdxBSTheta idx_trial_theta(theta_mod(IdxBSTheta(idx_trial)));

        auto& bspl_r = ddc::discrete_space<BSplinesR>();
        auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

        const Idx<KnotsR> start_non_zero_r(
                ::max(bspl_r.break_point_domain().front(),
                      ::max(bspl_r.get_first_support_knot(idx_test_r),
                            bspl_r.get_first_support_knot(idx_trial_r))));
        const Idx<KnotsR> end_non_zero_r(
                ::min(bspl_r.break_point_domain().back(),
                      ::min(bspl_r.get_last_support_knot(idx_test_r),
                            bspl_r.get_last_support_knot(idx_trial_r))));

        IdxStep<KnotsTheta> span_theta(BSplinesTheta::degree() + 1);

        Idx<KnotsTheta> first_support_knot_theta_test
                = bspl_theta.get_first_support_knot(idx_test_theta);
        Idx<KnotsTheta> first_support_knot_theta_trial
                = bspl_theta.get_first_support_knot(idx_trial_theta);
        Idx<KnotsTheta> last_support_knot_theta_test = first_support_knot_theta_test + span_theta;
        Idx<KnotsTheta> last_support_knot_theta_trial = first_support_knot_theta_trial + span_theta;

        if (first_support_knot_theta_test > last_support_knot_theta_trial) {
            first_support_knot_theta_trial += ddc::discrete_space<BSplinesTheta>().nbasis();
            last_support_knot_theta_trial += ddc::discrete_space<BSplinesTheta>().nbasis();
        } else if (last_support_knot_theta_test < first_support_knot_theta_trial) {
            first_support_knot_theta_test += ddc::discrete_space<BSplinesTheta>().nbasis();
            last_support_knot_theta_test += ddc::discrete_space<BSplinesTheta>().nbasis();
        }
        const Idx<KnotsTheta> start_non_zero_theta(
                ::max(first_support_knot_theta_test, first_support_knot_theta_trial));
        const Idx<KnotsTheta> end_non_zero_theta(
                ::min(last_support_knot_theta_test, last_support_knot_theta_trial));

        const IdxRangeQuadratureRTheta quad_range = get_quadrature_between_knots(
                start_non_zero_r,
                end_non_zero_r,
                start_non_zero_theta,
                end_non_zero_theta,
                m_idxrange_quadrature.front());

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = get_field(m_int_volume_alloc);

        const IdxBSPolar idx_test_polar(to_polar(idx_test));
        const IdxBSPolar idx_trial_polar(to_polar(idx_trial));

        return ddc::parallel_transform_reduce(
                quad_range,
                0.0,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(IdxQuadratureRTheta idx_quad) {
                    // Manage periodicity
                    if (!full_quad_idx_range.contains(idx_quad)) {
                        idx_quad -= full_quad_idx_range.template extent<QDimThetaMesh>();
                    }
                    assert(full_quad_idx_range.contains(idx_quad));
                    return weak_integral_element(
                            idx_test_polar,
                            idx_trial_polar,
                            idx_quad,
                            coeff_alpha,
                            coeff_beta,
                            evaluator,
                            mapping,
                            int_volume_proxy);
                });
    }

    /**
     * @brief Calculates the modulo idx_theta in relation to cells number along  @f$ \theta @f$ direction .
     *
     * @param[in] idx_theta @f$ \theta @f$ index.
     *
     * @return The corresponding indice modulo @f$ \theta @f$ direction cells number
     */
    static KOKKOS_FUNCTION IdxStepBSTheta theta_mod(IdxStepBSTheta idx_theta)
    {
        int n_theta = ddc::discrete_space<BSplinesTheta>().nbasis();
        while (idx_theta < 0)
            idx_theta += n_theta;
        while (idx_theta >= n_theta)
            idx_theta -= n_theta;
        return idx_theta;
    }

    /**
     * @brief Calculates the index which is inside the poloidal domain using the periodicity properties.
     *
     * @param[in] idx A multi-dimensional index including the polar bspline index.
     *
     * @return The corresponding index inside the domain.
     */
    template <class IdxType>
    static KOKKOS_INLINE_FUNCTION IdxType theta_mod(IdxType idx)
    {
        static_assert(ddc::is_discrete_element_v<IdxType>);
        static_assert(ddc::in_tags_v<BSplinesTheta, ddc::to_type_seq_t<IdxType>>);
        IdxRangeBSTheta idx_range_theta
                = ddc::discrete_space<BSplinesTheta>().full_domain().take_first(
                        IdxStepBSTheta(ddc::discrete_space<BSplinesTheta>().nbasis()));
        while (ddc::select<BSplinesTheta>(idx) < idx_range_theta.front())
            idx += idx_range_theta.extents();
        while (ddc::select<BSplinesTheta>(idx) > idx_range_theta.back())
            idx -= idx_range_theta.extents();
        assert(idx_range_theta.contains(ddc::select<BSplinesTheta>(idx)));
        return idx;
    }

    /**
     * @brief Get the value and derivative of the specified polar bspline at the specified quadrature point.
     *
     * This method calculates the value and the derivatives of polar bsplines. It is templated by
     * calculate_derivs to avoid code duplication between get_polar_bspline_vals_and_derivs and
     * get_polar_bspline_vals. The calling method should not need to use the template parameter.
     *
     * @param[out] val
     *      The value of the specified polar bspline at the specified point.
     * @param[in] coord
     *      The coordinate where the value of the polar bspline should be calculated.
     * @param[in] idx
     *      The polar bspline of interest.
     * @return The derivative of the polar bspline (only returned if calculate_derivs is true).
     */
    template <bool calculate_derivs = true>
    static KOKKOS_FUNCTION auto get_polar_bspline_vals_and_derivs(
            double& val,
            CoordRTheta coord,
            IdxBSPolar idx)
    {
        std::array<double, PolarBSplinesRTheta::n_singular_basis()> singular_data;
        std::array<double, m_n_non_zero_bases_r * m_n_non_zero_bases_theta> data;
        // Values of the polar basis splines around the singular point
        // at a given coordinate
        DSpan1D singular_vals(singular_data.data(), PolarBSplinesRTheta::n_singular_basis());
        // Values of the polar basis splines, that do not cover the singular point,
        // at a given coordinate
        DSpan2D vals(data.data(), m_n_non_zero_bases_r, m_n_non_zero_bases_theta);

        auto& polar_bspl = ddc::discrete_space<PolarBSplinesRTheta>();

        if (idx < IdxBSPolar(PolarBSplinesRTheta::n_singular_basis())) {
            IdxStepBSPolar offset
                    = idx
                      - PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>()
                                .front();
            polar_bspl.eval_basis(singular_vals, vals, coord);
            val = singular_vals[offset.value()];
            if constexpr (calculate_derivs) {
                polar_bspl.eval_deriv_r(singular_vals, vals, coord);
                double r_deriv = singular_vals[offset.value()];
                polar_bspl.eval_deriv_theta(singular_vals, vals, coord);
                double theta_deriv = singular_vals[offset.value()];
                DVector<R_cov, Theta_cov> derivs(r_deriv, theta_deriv);
                return derivs;
            } else {
                return;
            }
        } else {
            IdxBSRTheta idx_front = polar_bspl.eval_basis(singular_vals, vals, coord);
            IdxStepBSRTheta offset = PolarBSplinesRTheta::get_2d_index(idx) - idx_front;
            IdxStepBSR ir(offset);
            IdxStepBSTheta itheta(theta_mod(IdxStepBSTheta(offset)));

            val = vals(ir, itheta);
            if constexpr (calculate_derivs) {
                polar_bspl.eval_deriv_r(singular_vals, vals, coord);
                double r_deriv = vals(ir, itheta);
                polar_bspl.eval_deriv_theta(singular_vals, vals, coord);
                double theta_deriv = vals(ir, itheta);
                DVector<R_cov, Theta_cov> derivs(r_deriv, theta_deriv);
                return derivs;
            } else {
                return;
            }
        }
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
        get_polar_bspline_vals_and_derivs<false>(val, coord, idx);
        return val;
    }

    /**
     * @brief Fills the nnz data structure by computing the number of non-zero per line.
     * This number is linked to the weak formulation and depends on @f$(r,\theta)@f$ splines.
     * After this function the array will contain:
     * nnz_per_row[0] = 0.
     * nnz_per_row[1] = 0.
     * nnz_per_row[2] = number of non-zero elements in line 0.
     * nnz_per_row[3] = number of non-zero elements in lines 0-1.
     * ..._per_row
     * nnz_per_row[matrix_size] = number of non-zero elements in lines 0-(matrix_size-1).
     *
     * @param[out]  nnz_per_row A 1D Kokkos view of length matrix_size+1 which stores the sum of the non-zeros in the matrix on all lines up
     *                          to the one in.
     */
    void init_nnz_per_line(Kokkos::View<int*, Kokkos::LayoutRight> nnz_per_row) const
    {
        Kokkos::Profiling::pushRegion("PolarPoissonInitNnz");
        size_t const mat_size = nnz_per_row.extent(0) - 1;
        IdxStepBSPolar radial_boundary_splines(m_nbasis_theta);
        IdxRangeBSPolar polar_bspl_idx_range
                = ddc::discrete_space<PolarBSplinesRTheta>().full_domain().remove_last(
                        radial_boundary_splines);

        Kokkos::deep_copy(nnz_per_row, 0.);

        assert(mat_size == polar_bspl_idx_range.size());
        // nnz per line
        Field<int, IdxRangeBSPolar>
                nnz(Kokkos::subview(nnz_per_row, std::pair<int, int>(2, mat_size + 1)),
                    polar_bspl_idx_range.remove_last(IdxStepBSPolar(1)));

        size_t constexpr n_singular_basis = PolarBSplinesRTheta::n_singular_basis();
        size_t constexpr degree = BSplinesR::degree();
        size_t constexpr stencil_overlap = 2 * degree + 1;
        size_t const nbasis_theta_proxy = m_nbasis_theta;

        // The number of radial splines which overlap with a given radial spline to its
        // left or to its right
        IdxStepBSR n_one_side_overlap(BSplinesR::degree());

        int nnz_for_sigular_rows = n_singular_basis + degree * nbasis_theta_proxy;

        // rows representing the bsplines which cover the singular domain
        // These overlap with other singular domain splines and degree radial splines
        Kokkos::parallel_for(
                "overlap singular radial",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(1, n_singular_basis + 1),
                KOKKOS_LAMBDA(const int k) { nnz_per_row(k + 1) = k * nnz_for_sigular_rows; });

        int nnz_sum = nnz_for_sigular_rows * n_singular_basis;

        IdxRangeBSR idxrange_bsplines_r_overlapping_singular
                = m_idxrange_bsplines_r.take_first(n_one_side_overlap);
        IdxRangeBSTheta idxrange_bsplines_theta = m_idxrange_bsplines_theta;
        IdxRangeBSRTheta idxrange_singular_overlap(
                idxrange_bsplines_r_overlapping_singular,
                idxrange_bsplines_theta);
        // rows representing the bsplines which can be represented as a tensor product of
        // radial and poloidal bsplines but which still overlap with the singular domain
        Kokkos::parallel_for(
                "inner overlap",
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        idxrange_bsplines_r_overlapping_singular.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    IdxStepBSR ir(team.league_rank());
                    // number of radial bases which overlap with this element
                    int nr = degree + ir + 1;

                    // nnz to this line = nnz_singular + n_singular * ir*ntheta + ntheta * n_stencil_theta \sum_{k=degree+1}^{nr-1} k
                    int nnz_sum_to_r = nnz_sum + n_singular_basis * ir * nbasis_theta_proxy
                                       + (nr * (nr - 1) - degree * (degree + 1))
                                                 * nbasis_theta_proxy * stencil_overlap / 2;

                    // Loop over poloidal dimensions
                    Kokkos::parallel_for(
                            Kokkos::TeamThreadRange(team, nbasis_theta_proxy),
                            [&](int const& thread_index) {
                                IdxStepBSTheta itheta(thread_index);
                                IdxBSRTheta k_2d = idxrange_singular_overlap.front() + ir + itheta;
                                IdxBSPolar k(to_polar(k_2d));
                                nnz(k) = nnz_sum_to_r
                                         + (thread_index + 1)
                                                   * (n_singular_basis + nr * stencil_overlap);
                            });
                });

        {
            int nr = idxrange_bsplines_r_overlapping_singular.size();
            nnz_sum += (n_singular_basis * nr + nr * (nr + 2 * degree + 1) * stencil_overlap / 2)
                       * nbasis_theta_proxy;
        }

        IdxRangeBSR idxrange_bsplines_r_stencil
                = m_idxrange_bsplines_r.remove(n_one_side_overlap, n_one_side_overlap + 1);
        IdxRangeBSRTheta idxrange_stencil(idxrange_bsplines_r_stencil, idxrange_bsplines_theta);
        IdxBSPolar idx_stencil_front = to_polar(idxrange_stencil.front());
        // Stencil for tensor product bsplines which only overlap with other tensor product bsplines
        Kokkos::parallel_for(
                "Inner Stencil",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, idxrange_stencil.size()),
                KOKKOS_LAMBDA(const int k) {
                    IdxBSPolar kp(idx_stencil_front + IdxStepBSPolar(k));
                    nnz(kp) = nnz_sum + (k + 1) * stencil_overlap * stencil_overlap;
                });

        nnz_sum += idxrange_stencil.size() * stencil_overlap * stencil_overlap;

        // Get the outer radial bsplines excluding the bspline used for homogeneous Hermite boundary conditions
        IdxRangeBSR outer_bsplines_r = m_idxrange_bsplines_r.take_last(n_one_side_overlap + 1)
                                               .remove_last(IdxStepBSR(1));
        IdxRangeBSRTheta outer_bsplines_2d(outer_bsplines_r, idxrange_bsplines_theta);
        IdxRangeBSPolar outer_bsplines(
                to_polar(outer_bsplines_2d.front()),
                IdxStepBSPolar(outer_bsplines_2d.size() - 1));
        assert(outer_bsplines.back() == get_idx_range(nnz).back());
        // Approaching the external boundary the overlapping possibilities between two radial splines decrease
        Kokkos::parallel_for(
                "outer overlap",
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        outer_bsplines_r.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    IdxStepBSR ir(team.league_rank());
                    int nr = 2 * degree - ir;
                    IdxBSR k_r(outer_bsplines_2d.front() + ir);
                    int nnz_sum_to_r = nnz_sum
                                       + (degree * (2 * degree + 1) - nr * (nr + 1) / 2)
                                                 * stencil_overlap * nbasis_theta_proxy;
                    // Loop over poloidal dimensions
                    Kokkos::parallel_for(
                            Kokkos::TeamThreadRange(team, nbasis_theta_proxy),
                            [&](int const& thread_index) {
                                IdxStepBSTheta itheta(thread_index);
                                IdxBSRTheta k_2d = outer_bsplines_2d.front() + ir + itheta;
                                IdxBSPolar k(to_polar(k_2d));
                                if (outer_bsplines.contains(k)) {
                                    nnz(k) = nnz_sum_to_r
                                             + (thread_index + 1) * nr * stencil_overlap;
                                }
                            });
                });

        Kokkos::Profiling::popRegion();
    }

    /**
     * Convert a 2D (r,theta) bspline index into a polar bspline index.
     *
     * @param[in] idx The 2D (r,theta) bspline index.
     * @return The polar bspline index.
     */
    static KOKKOS_INLINE_FUNCTION IdxBSPolar to_polar(IdxBSRTheta idx)
    {
        return PolarBSplinesRTheta::template get_polar_index<PolarBSplinesRTheta>(idx);
    }

    /**
     * @brief Compute the quadrature range between a provided set of knots.
     *
     * Compute the range of quadrature points which are found between a set of knots
     * in both the radial and poloidal directions. In order to return a contiguous range
     * the result may include indices which are outside the domain. A modulo operator
     * should be applied before using the indices.
     *
     * @param[in] start_knot_r
     *      The index of the knot describing the lower bound of the domain of interest
     *      in the radial direction.
     * @param[in] end_knot_r
     *      The index of the knot describing the upper bound of the domain of interest
     *      in the radial direction.
     * @param[in] start_knot_theta
     *      The index of the knot describing the lower bound of the domain of interest
     *      in the poloidal direction.
     * @param[in] end_knot_theta
     *      The index of the knot describing the upper bound of the domain of interest
     *      in the poloidal direction.
     * @param[in] idx_quad_front
     *      The first index of the index range of the quadrature points.
     * @return 
     *      The range of quadrature points in the specified domain.
     */
    static KOKKOS_FUNCTION IdxRangeQuadratureRTheta get_quadrature_between_knots(
            Idx<KnotsR> start_knot_r,
            Idx<KnotsR> end_knot_r,
            Idx<KnotsTheta> start_knot_theta,
            Idx<KnotsTheta> end_knot_theta,
            IdxQuadratureRTheta idx_quad_front)
    {
        const IdxRange<KnotsR> k_range_r(start_knot_r, end_knot_r - start_knot_r);
        const IdxRange<KnotsTheta>
                k_range_theta(start_knot_theta, end_knot_theta - start_knot_theta);

        IdxStep<KnotsR> k_r_offset
                = k_range_r.front() - ddc::discrete_space<BSplinesR>().break_point_domain().front();
        IdxQuadratureR q_r_offset
                = IdxQuadratureR(idx_quad_front) + k_r_offset.value() * s_n_gauss_legendre_r;
        IdxStepQuadratureR q_r_len(k_range_r.extents().value() * s_n_gauss_legendre_r);
        IdxRangeQuadratureR q_range_r(q_r_offset, q_r_len);

        IdxStep<KnotsTheta> k_theta_offset
                = k_range_theta.front()
                  - ddc::discrete_space<BSplinesTheta>().break_point_domain().front();
        if (k_theta_offset < 0)
            k_theta_offset += ddc::discrete_space<BSplinesTheta>().nbasis();
        IdxQuadratureTheta q_theta_offset = IdxQuadratureTheta(idx_quad_front)
                                            + k_theta_offset.value() * s_n_gauss_legendre_theta;
        IdxStepQuadratureTheta q_theta_len(
                k_range_theta.extents().value() * s_n_gauss_legendre_theta);
        IdxRangeQuadratureTheta q_range_theta(q_theta_offset, q_theta_len);
        assert(q_range_r.extents() > 0);
        assert(q_range_theta.extents() > 0);

        return IdxRangeQuadratureRTheta(q_range_r, q_range_theta);
    }
};
