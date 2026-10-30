// SPDX-License-Identifier: MIT
#pragma once

#include <tuple>

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "fem_poisson_assembler_utils.hpp"
#include "gauss_legendre_integration.hpp"
#include "math_tools.hpp"
#include "matrix_batch_csr.hpp"
#include "metric_tensor_evaluator.hpp"
#include "polar_spline_evaluator.hpp"
#include "quadrature_coeffs_nd.hpp"
#include "tensor_product_bspline_poisson_assembler.hpp"
#include "vector_field.hpp"
#include "view.hpp"
#include "volume_quadrature_nd.hpp"


namespace detail_poisson {

/**
 * @brief Get the value and derivatives of a polar bspline which covers the singular point
 *      at the specified quadrature point.
 *
 * The polar bsplines which cover the singular point are linear combinations of tensor product
 * bsplines so they are evaluated using PolarBSplines::eval_basis_and_n_derivs. The polar
 * bsplines which do not cover the singular point should be passed as a 2D index (see the
 * overload below).
 *
 * @tparam calculate_derivs If true, the value and the derivatives are calculated. If false,
 *     only the value is calculated.
 *
 * @param[in] idx_quad
 *      The quadrature point where the polar bspline should be evaluated.
 * @param[in] idx
 *      The index of a polar bspline which covers the singular point.
 * @return The value of the polar bspline if calculate_derivs is false. Otherwise a tuple
 *      containing the value and the derivatives of the polar bspline.
 */
template <bool calculate_derivs = true, class PolarBSplinesRTheta, class IdxQuad>
KOKKOS_FUNCTION std::conditional_t<
        calculate_derivs,
        std::tuple<
                double,
                DVector<typename PolarBSplinesRTheta::R::Dual,
                        typename PolarBSplinesRTheta::Theta::Dual>>,
        double>
get_polar_bspline_vals_and_derivs(IdxQuad idx_quad, Idx<PolarBSplinesRTheta> idx)
{
    using R = PolarBSplinesRTheta::R;
    using Theta = PolarBSplinesRTheta::Theta;
    using BSplinesR = PolarBSplinesRTheta::BSplinesR_tag;
    using BSplinesTheta = PolarBSplinesRTheta::BSplinesTheta_tag;
    constexpr std::size_t n_singular = PolarBSplinesRTheta::n_singular_basis();

    KOKKOS_ASSERT(idx < Idx<PolarBSplinesRTheta>(n_singular));

    constexpr int n_non_zero_bases_r = BSplinesR::degree() + 1;
    constexpr int n_non_zero_bases_theta = BSplinesTheta::degree() + 1;
    // The maximum derivative order in each direction
    constexpr std::size_t n_derivs = calculate_derivs ? 1 : 0;
    constexpr std::size_t n_deriv_combinations = (n_derivs + 1) * (n_derivs + 1);
    std::array<double, n_singular * n_deriv_combinations> singular_data;
    std::array<double, n_non_zero_bases_r * n_non_zero_bases_theta * n_deriv_combinations> data;
    // Values and derivatives of the polar basis splines around the singular point
    // at a given coordinate
    DSpan3D singular_vals(singular_data.data(), n_singular, n_derivs + 1, n_derivs + 1);
    // Values and derivatives of the polar basis splines, that do not cover the singular point,
    // at a given coordinate
    DSpan4D
            vals(data.data(),
                 n_non_zero_bases_r,
                 n_non_zero_bases_theta,
                 n_derivs + 1,
                 n_derivs + 1);

    Coord<R, Theta> coord = ddc::coordinate(idx_quad);
    ddc::discrete_space<PolarBSplinesRTheta>()
            .eval_basis_and_n_derivs(singular_vals, vals, coord, n_derivs, n_derivs);

    IdxStep<PolarBSplinesRTheta> offset
            = idx - PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>().front();
    const double val = singular_vals(offset.value(), 0, 0);
    if constexpr (calculate_derivs) {
        return {val,
                DVector<typename R::Dual, typename Theta::Dual>(
                        singular_vals(offset.value(), 1, 0),
                        singular_vals(offset.value(), 0, 1))};
    } else {
        return val;
    }
}

} // namespace detail_poisson


/**
 * @brief An operator to assemble a Poisson-like stiffness matrix using polar B-splines.
 * 
 * Assemble the finite element stiffness matrix with entries 
 * @f$ \int \alpha \nabla B_i \cdot \nabla B_j + \beta B_i B_j dx @f$
 * for @f$0 < i,j < N @f$ where @f$ N @f$ is the number of polar splines. We only consider 
 * splines that vanish on the boundary.
 * 
 * For more details see the `PolarSplineFEMPoissonLikeSolver` and Emily Bourne's thesis 
 * ("Non-Uniform Numerical Schemes for the Modelling of Turbulence
 * in the 5D GYSELA Code". December 2022.)
 * @tparam GridR The radial grid type.
 * @tparam GridR The poloidal grid type.
 * @tparam PolarBSplinesRTheta The type of the 2D polar B-splines (on the coordinate
 *          system @f$(r,\theta)@f$ including B-splines which traverse the O point).
 * @tparam QDimRMesh The radial quadrature grid type.
 * @tparam QDimThetaMesh The poloidal quadrature grid type.
 * @tparam MetricTensorType The type of the metric tensor associated with the mapping from the
 *          logical domain to the physical domain where the equation is defined.
 * @tparam IdxRangeFull The full index range of @f$ \phi @f$ including any batch dimensions.
 *          One matrix is assembled for each index in the batch dimensions.
 */

template <
        typename GridR,
        typename GridTheta,
        typename PolarBSplinesRTheta,
        typename QDimRMesh,
        typename QDimThetaMesh,
        concepts::MetricTensor MetricTensorType,
        typename IdxRangeFull = IdxRange<GridR, GridTheta>>
class PolarSplineFEMPoissonLikeAssembler
{
public:
    /// The radial dimension
    using R = typename GridR::continuous_dimension_type;
    /// The poloidal dimension
    using Theta = typename GridTheta::continuous_dimension_type;

    static_assert(R::IS_CONTRAVARIANT);
    static_assert(Theta::IS_CONTRAVARIANT);
    static_assert(
            std::is_same_v<typename MetricTensorType::CoordArg, Coord<R, Theta>>,
            "The metric tensor must be evaluated on (r, theta) coordinates.");

private:
    /// The radial dimension
    using R_cov = typename R::Dual;
    /// The poloidal dimension
    using Theta_cov = typename Theta::Dual;

private:
    /// The 1D B-splines in the radial direction
    using BSplinesR = typename PolarBSplinesRTheta::BSplinesR_tag;
    /// The 1D B-splines in the poloidal direction
    using BSplinesTheta = typename PolarBSplinesRTheta::BSplinesTheta_tag;

    using KnotsR = ddc::knot_discrete_dimension_t<BSplinesR>;
    using KnotsTheta = ddc::knot_discrete_dimension_t<BSplinesTheta>;

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

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFull, GridR, GridTheta>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;

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

    /// The assembler for the inner products of tensor basis functions
    using TensorAssembler = TensorProductBSplinePoissonAssembler<
            GridR,
            GridTheta,
            BSplinesR,
            BSplinesTheta,
            QDimRMesh,
            QDimThetaMesh,
            MetricTensorType,
            IdxRangeFull>;

public:
    /// The values and derivatives of the 1D B-splines at the quadrature points.
    using TensorBSplineStash = typename TensorAssembler::TensorBSplineStash;

private:
    static constexpr int s_n_gauss_legendre_r = BSplinesR::degree() + 1;
    // The number of cells (in the radial direction) in which both types of basis splines can be found
    static constexpr int s_n_overlap_cells = PolarBSplinesRTheta::continuity + 1;

    const int m_nbasis_theta;

    // Matrix size is equal to the number of Polar bspline
    const int m_matrix_size;

    // Domains
    IdxRangeBSR m_idxrange_bsplines_r;
    IdxRangeBSTheta m_idxrange_bsplines_theta;

    // The index of the first row/column of the tensor-tensor sub-matrix in the matrix
    const int m_tensor_idx_offset;

    IdxRangeBatch m_idx_range_batch;

    MetricTensorType m_metric_tensor;

    // The assembler for the tensor-tensor sub-matrix. This also stores the quadrature and
    // B-spline data shared by all sub-matrices.
    TensorAssembler m_tensor_assembler;

public:
    /**
     * @brief Instantiate the assembler operator.
     *
     * @param int_volume The initialised field of Jacobian values of the mapping.
     * @param metric_tensor The metric tensor associated with the mapping from the logical
     *          domain to the physical domain where the equation is defined.
     * @param idx_range_full The full index range of @f$ \phi @f$ including any batch dimensions.
     * @param include_boundary_splines True if the outer radial bsplines should be included in
     *          the matrix assembly, false otherwise.
     */
    explicit PolarSplineFEMPoissonLikeAssembler(
            Field<double, IdxRangeQuadratureRTheta> int_volume,
            MetricTensorType const& metric_tensor,
            IdxRangeFull idx_range_full,
            bool include_boundary_splines = true)
        : m_nbasis_theta(ddc::discrete_space<BSplinesTheta>().nbasis())
        , m_matrix_size(get_matrix_size(include_boundary_splines))
        , m_idxrange_bsplines_r(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                  IdxStepBSR {s_n_overlap_cells}))
        , m_idxrange_bsplines_theta(ddc::discrete_space<BSplinesTheta>().full_domain().take_first(
                  IdxStepBSTheta {m_nbasis_theta}))
        , m_tensor_idx_offset(PolarBSplinesRTheta::n_singular_basis())
        , m_idx_range_batch(idx_range_full)
        , m_metric_tensor(metric_tensor)
        // The radial bsplines exclude the bspline used for homogeneous Hermite boundary conditions
        , m_tensor_assembler(
                  int_volume,
                  metric_tensor,
                  idx_range_full,
                  IdxRangeBSRTheta(
                          m_idxrange_bsplines_r.remove_last(IdxStepBSR(include_boundary_splines)),
                          m_idxrange_bsplines_theta))
    {
    }

    /**
     * @brief Get the TensorBSplineStash of the values and derivatives of the 1D B-splines at the quadrature points.
     *
     * @return The stashed values of the 1D B-splines.
     */
    TensorBSplineStash get_tensor_bspline_stash() const
    {
        return m_tensor_assembler.get_tensor_bspline_stash();
    }

    /**
     * @brief Compute the sparsity pattern for the stiffness matrix.
     * 
     * @param[out] gko_matrix The pointer to the assembled matrix.
     * @param[in] max_iter
     *      The maximum number of iterations possible for the batched CSR solver.
     * @param[in] res_tol
     *      The residual tolerance for the batched CSR solver. The solver stops when
     *      ||r||/||b|| <= res_tol where r is the residual computed by the CG recurrence.
     * @param[in] batch_solver_logger
     *      Indicates whether log information such as the residual and the number of iterations
     *      should be monitored.
     * @param[in] preconditioner_max_block_size
     *      The maximum size of the Jacobi preconditioner used by the batched CSR solver.
     * @param[in] residual_norm
     *      The norm used to evaluate the residual in the stopping criterion.
     *
     * @tparam Mapping A class describing a mapping from curvilinear coordinates to Cartesian coordinates.
     */
    void setup_sparse_matrix(
            std::unique_ptr<
                    MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>&
                    gko_matrix,
            std::optional<int> max_iter = std::nullopt,
            std::optional<double> res_tol = std::nullopt,
            std::optional<bool> batch_solver_logger = std::nullopt,
            std::optional<int> preconditioner_max_block_size = std::nullopt,
            std::optional<GkoMatrixResidualNorm> residual_norm = std::nullopt)
    {
        const int n_matrix_elements
                = n_singular_elements() + m_tensor_assembler.n_tensor_tensor_elements();

        //CSR data storage
        gko_matrix = std::make_unique<
                MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>(
                m_idx_range_batch.size(),
                m_matrix_size,
                n_matrix_elements,
                max_iter,
                res_tol,
                batch_solver_logger,
                preconditioner_max_block_size,
                residual_norm);
        auto [values, col_idx, nnz_per_row] = gko_matrix->get_batch_csr();
        init_nnz_per_line(nnz_per_row);
        compute_singular_singular_col_idx(col_idx, nnz_per_row);
        compute_singular_tensor_col_idx(col_idx, nnz_per_row);
        m_tensor_assembler
                .compute_col_idx(col_idx, get_tensor_nnz_per_row(nnz_per_row), m_tensor_idx_offset);
    }

    /**
     * @brief Get the number of non-zero elements in the matrix which involve a singular
     * basis function.
     *
     * These are the elements corresponding to the inner products of singular basis functions
     * and singular basis functions, and to the inner products of singular basis functions and
     * tensor basis functions (and vice versa).
     *
     * @return The number of non-zero elements involving a singular basis function.
     */
    int n_singular_elements() const
    {
        // Number of elements in the matrix that correspond to the splines
        // that cover the singular point
        constexpr int n_elements_singular
                = PolarBSplinesRTheta::n_singular_basis() * PolarBSplinesRTheta::n_singular_basis();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // polar bsplines which traverse the singular point and other polar bsplines
        const int n_elements_overlap = 2
                                       * (PolarBSplinesRTheta::n_singular_basis()
                                          * BSplinesR::degree() * m_nbasis_theta);
        return n_elements_singular + n_elements_overlap;
    }

    /**
     * @brief Get the elements of the nnz data structure of the matrix which correspond to the
     * rows of the tensor-tensor sub-matrix.
     *
     * @param[in] nnz_per_row A 1D Kokkos view of length matrix_size+1 containing the nnz data
     *                  structure of the matrix.
     * @return A 1D Kokkos view of length n+1, where n is the number of tensor basis functions.
     */
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> get_tensor_nnz_per_row(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> nnz_per_row)
            const
    {
        int const n_tensor_basis = m_tensor_assembler.get_tensor_basis_idx_range().size();
        return Kokkos::subview(
                nnz_per_row,
                std::pair<int, int>(m_tensor_idx_offset, m_tensor_idx_offset + n_tensor_basis + 1));
    }

    /**
     * @brief Assemble the stiffness matrix.
     *
     * @param[out] gko_matrix The pointer to the assembled matrix.
     * @param[in] coeff_alpha
     *      A callable object with signature `double operator()(IdxBatch, CoordRTheta)`
     *      returning the value of @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable object with signature `double operator()(IdxBatch, CoordRTheta)`
     *      returning the value of @f$ \beta @f$ at the given batch and coordinate.
     *
     * @tparam CoeffAlpha A callable type for evaluating @f$ \alpha @f$ at a batch and coordinate.
     * @tparam CoeffBeta A callable type for evaluating @f$ \beta @f$ at a batch and coordinate.
     */
    template <class CoeffAlpha, class CoeffBeta>
    void operator()(
            std::unique_ptr<
                    MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>> const&
                    gko_matrix,
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta)
    {
        //CSR data storage
        auto [values, col_idx, nnz_per_row] = gko_matrix->get_batch_csr();

        compute_singular_singular_elements(
                coeff_alpha,
                coeff_beta,
                m_metric_tensor,
                values,
                col_idx,
                nnz_per_row);
        compute_singular_tensor_elements(
                coeff_alpha,
                coeff_beta,
                m_metric_tensor,
                values,
                col_idx,
                nnz_per_row);
        m_tensor_assembler.compute_elements(
                coeff_alpha,
                coeff_beta,
                values,
                col_idx,
                get_tensor_nnz_per_row(nnz_per_row),
                m_tensor_idx_offset);

        gko_matrix->setup_solver();
    }

    /**
     * @brief Computes the column indices of the matrix elements corresponding to the
     * inner products of singular basis functions and singular basis functions.
     *
     * @param[out] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix).
     * @param[in] nnz_per_row_csr
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    void compute_singular_singular_col_idx(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    col_idx_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    nnz_per_row_csr)
    {
        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();
        const int n_singular = idxrange_singular.size();

        const std::source_location location = std::source_location::current();
        // Calculate the matrix elements corresponding to the B-splines which cover the singular point
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::MDRangePolicy<
                        Kokkos::Rank<2>,
                        Kokkos::DefaultExecutionSpace>({0, 0}, {n_singular, n_singular}),
                KOKKOS_LAMBDA(const int row_idx, const int col_idx) {
                    const int csr_idx_singular_area = nnz_per_row_csr(row_idx) + col_idx;
                    //Fill the dense matrix corresponding to the b-splines on the singular point
                    col_idx_csr(csr_idx_singular_area) = col_idx;
                });
    }

    /**
     * @brief Computes the column indices of the matrix elements corresponding to the
     * inner products of singular basis functions and tensor basis functions.
     *
     * @param[out] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @param[in] nnz_per_row_csr
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    void compute_singular_tensor_col_idx(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    col_idx_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    nnz_per_row_csr)
    {
        // Create index ranges associated with the 2D splines
        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();
        IdxRangeBSR central_radial_bspline_idx_range(
                m_idxrange_bsplines_r.take_first(IdxStep<BSplinesR> {BSplinesR::degree()}));

        IdxRangeBSRTheta idxrange_non_singular_near_centre(
                central_radial_bspline_idx_range,
                m_idxrange_bsplines_theta);

        const int n_singular = idxrange_singular.size();
        // Number of tensor product polar bsplines which overlap with polar bsplines which
        // traverse the singular point
        const int n_overlapping_singular = idxrange_non_singular_near_centre.size();

        const std::source_location location = std::source_location::current();

        // Calculate the matrix elements where bspline products overlap the B-splines which cover the singular point
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        n_singular * n_overlapping_singular,
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const IdxStepBSPolar idx_step_test(team.league_rank() / n_overlapping_singular);
                    const IdxStepBSPolar idx_step_trial(
                            team.league_rank() % n_overlapping_singular);

                    const int row_idx = idx_step_test.value();
                    const int col_idx = idx_step_trial.value() + n_singular;

                    //a_ij
                    col_idx_csr(nnz_per_row_csr(row_idx) + col_idx) = col_idx;
                    //a_ji
                    col_idx_csr(nnz_per_row_csr(col_idx) + row_idx) = row_idx;
                });
    }

    /**
     * @brief Computes the matrix elements corresponding to the inner products of singular
     * basis functions and singular basis functions.
     *
     * @param[in] coeff_alpha
     *      A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning
     *      @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning
     *      @f$ \beta @f$ at the given batch and coordinate.
     * @param[in] metric_tensor
     *      The metric tensor associated with the mapping from the logical domain to the
     *      physical domain where the equation is defined.
     * @param[out] values_csr
     *             A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[in] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix).
     * @param[in] nnz_per_row_csr
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    template <class CoeffAlpha, class CoeffBeta>
    void compute_singular_singular_elements(
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            MetricTensorType const& metric_tensor,
            Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    values_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    col_idx_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    nnz_per_row_csr)
    {
        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();
        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_tensor_assembler.get_int_volume();

        IdxRangeQuadratureRTheta full_quad_idx_range = get_idx_range(int_volume_proxy);
        IdxRangeQuadratureRTheta idx_range_quad_singular(
                IdxRangeQuadratureR(full_quad_idx_range)
                        .take_first(IdxStep<QDimRMesh> {s_n_overlap_cells * s_n_gauss_legendre_r}),
                IdxRangeQuadratureTheta(full_quad_idx_range));

        IdxRangeBatch idx_range_batch = m_idx_range_batch;
        const int n_batch = idx_range_batch.size();
        const int n_singular = idxrange_singular.size();

        const std::source_location location = std::source_location::current();
        // Calculate the matrix elements corresponding to the B-splines which cover the singular point
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        n_batch * n_singular * n_singular,
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const int batch_linear_idx = team.league_rank() / (n_singular * n_singular);
                    const IdxBatch idx_batch = ddcHelper::
                            get_idx_from_linear_index(idx_range_batch, batch_linear_idx);
                    const int row_col_idx = team.league_rank() % (n_singular * n_singular);
                    const int row_idx = row_col_idx / n_singular;
                    const int col_idx = row_col_idx % n_singular;
                    IdxBSPolar idx_test = idxrange_singular.front() + IdxStepBSPolar(row_idx);
                    IdxBSPolar idx_trial = idxrange_singular.front() + IdxStepBSPolar(col_idx);
                    double element = 0;
                    // Calculate the weak integral
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
                                auto [test_mass_weight, test_stiffness_weight] = get_test_weights(
                                        idx_test,
                                        idx_quad,
                                        idx_batch,
                                        coeff_alpha,
                                        coeff_beta,
                                        metric_tensor,
                                        int_volume_proxy);

                                auto [basis_val_trial_space, basis_derivs_trial_space]
                                        = detail_poisson::get_polar_bspline_vals_and_derivs(
                                                idx_quad,
                                                idx_trial);

                                // Assemble the weak integral element
                                sum += scalar_product(
                                               test_stiffness_weight,
                                               basis_derivs_trial_space)
                                       + test_mass_weight * basis_val_trial_space;
                            },
                            element);
                    const int csr_idx_singular_area = nnz_per_row_csr(row_idx) + col_idx;
                    //Fill the dense matrix corresponding to the b-splines on the singular point
                    values_csr(batch_linear_idx, csr_idx_singular_area) = element;
                });
    }

    /**
     * @brief Computes the matrix elements corresponding to the inner products of singular
     * basis functions and tensor basis functions.
     *
     * @param[in] coeff_alpha
     *      A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning
     *      @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning
     *      @f$ \beta @f$ at the given batch and coordinate.
     * @param[in] metric_tensor
     *      The metric tensor associated with the mapping from the logical domain to the
     *      physical domain where the equation is defined.
     * @param[out] values_csr
     *             A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[in] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @param[in] nnz_per_row_csr
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    template <class CoeffAlpha, class CoeffBeta>
    void compute_singular_tensor_elements(
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            MetricTensorType const& metric_tensor,
            Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    values_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    col_idx_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    nnz_per_row_csr)
    {
        // Create index ranges associated with the 2D splines
        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();
        IdxRangeBSR central_radial_bspline_idx_range(
                m_idxrange_bsplines_r.take_first(IdxStep<BSplinesR> {BSplinesR::degree()}));

        IdxRangeBSRTheta idxrange_non_singular_near_centre(
                central_radial_bspline_idx_range,
                m_idxrange_bsplines_theta);

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_tensor_assembler.get_int_volume();
        TensorBSplineStash bspline_stash = get_tensor_bspline_stash();
        IdxRangeQuadratureRTheta full_quad_idx_range = get_idx_range(int_volume_proxy);

        IdxQuadratureRTheta idxrange_quadrature_front = full_quad_idx_range.front();

        IdxRangeBatch idx_range_batch = m_idx_range_batch;
        const int n_batch = idx_range_batch.size();
        // Number of polar bsplines which traverse the singular point
        const int n_singular = idxrange_singular.size();
        // Number of tensor product polar bsplines which overlap with polar bsplines which
        // traverse the singular point
        const int n_overlapping_singular = idxrange_non_singular_near_centre.size();

        const std::source_location location = std::source_location::current();

        // Calculate the matrix elements where bspline products overlap the B-splines which cover the singular point
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        n_batch * n_singular * n_overlapping_singular,
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const int batch_linear_idx(
                            team.league_rank() / (n_singular * n_overlapping_singular));
                    const IdxBatch idx_batch = ddcHelper::
                            get_idx_from_linear_index(idx_range_batch, batch_linear_idx);
                    const int idx_test_trial(
                            team.league_rank() % (n_singular * n_overlapping_singular));
                    const IdxStepBSPolar idx_step_test(idx_test_trial / n_overlapping_singular);
                    const IdxStepBSPolar idx_step_trial(idx_test_trial % n_overlapping_singular);
                    IdxBSPolar idx_test = idxrange_singular.front() + idx_step_test;
                    const IdxBSPolar idx_trial_polar(idxrange_singular.back() + 1 + idx_step_trial);
                    IdxBSRTheta idx_trial = PolarBSplinesRTheta::get_2d_index(idx_trial_polar);
                    IdxBSR idx_trial_r(idx_trial);
                    IdxBSTheta idx_trial_theta(idx_trial);

                    auto& bspl_r = ddc::discrete_space<BSplinesR>();
                    auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

                    // Find the index range covering the cells where both the test and trial functions are non-zero
                    const Idx<KnotsR> start_non_zero_r(
                            Kokkos::
                                    max(bspl_r.break_point_domain().front(),
                                        bspl_r.get_first_support_knot(idx_trial_r)));
                    const Idx<KnotsR> end_non_zero_r(
                            Kokkos::
                                    min(bspl_r.get_last_support_knot(
                                                IdxBSR(PolarBSplinesRTheta::continuity)),
                                        bspl_r.get_last_support_knot(idx_trial_r)));

                    const Idx<KnotsTheta> start_non_zero_theta(
                            bspl_theta.get_first_support_knot(idx_trial_theta));
                    const Idx<KnotsTheta> end_non_zero_theta(
                            bspl_theta.get_last_support_knot(idx_trial_theta));
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
                                    idxrange_quadrature_front);
                    assert(quad_range.size() > 0);
                    double element = 0;
                    // Calculate the weak integral
                    Kokkos::parallel_reduce(
                            Kokkos::TeamThreadMDRange(
                                    team,
                                    quad_range.template extent<QDimRMesh>(),
                                    quad_range.template extent<QDimThetaMesh>()),
                            [&](int r_thread_index, int theta_thread_index, double& sum) {
                                IdxQuadratureRTheta idx_quad = quad_range.front()
                                                               + IdxStep<QDimRMesh, QDimThetaMesh>(
                                                                       r_thread_index,
                                                                       theta_thread_index);
                                // Manage periodicity
                                if (!full_quad_idx_range.contains(idx_quad)) {
                                    idx_quad
                                            -= full_quad_idx_range.template extent<QDimThetaMesh>();
                                }

                                auto [test_mass_weight, test_stiffness_weight] = get_test_weights(
                                        idx_test,
                                        idx_quad,
                                        idx_batch,
                                        coeff_alpha,
                                        coeff_beta,
                                        metric_tensor,
                                        int_volume_proxy);

                                auto [basis_val_trial_space, basis_derivs_trial_space]
                                        = bspline_stash.vals_and_derivs(idx_quad, idx_trial);

                                // Assemble the weak integral element
                                sum += scalar_product(
                                               test_stiffness_weight,
                                               basis_derivs_trial_space)
                                       + test_mass_weight * basis_val_trial_space;
                            },
                            element);

                    const int row_idx = idx_step_test.value();
                    const int col_idx = idx_step_trial.value() + n_singular;

                    //a_ij
                    values_csr(batch_linear_idx, nnz_per_row_csr(row_idx) + col_idx) = element;
                    //a_ji
                    values_csr(batch_linear_idx, nnz_per_row_csr(col_idx) + row_idx) = element;
                });
    }

    /**
     * @brief Computes the test_stiffness_weight and test_mass_weight for a given quadrature point
     * for a polar test basis spline which covers the singular point.
     *
     * See detail_poisson::get_test_weights for more details.
     *
     * @param[in] idx_test
     *      The index of the test polar basis spline.
     * @param[in] idx_quad
     *      The index for the point in the quadrature scheme.
     * @param[in] idx_batch
     *      The index of the batch for which the element is being calculated.
     * @param[in] coeff_alpha
     *      A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning
     *      @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning
     *      @f$ \beta @f$ at the given batch and coordinate.
     * @param[in] metric_tensor
     *      The metric tensor associated with the mapping from the logical domain to the
     *      physical domain where the equation is defined.
     * @param[in] int_volume
     *      The integral volume associated with each point used in the quadrature scheme.
     * @return A tuple containing test_mass_weight and test_stiffness_weight.
     */
    template <class CoeffAlpha, class CoeffBeta>
    static KOKKOS_FUNCTION std::tuple<double, DVector<R, Theta>> get_test_weights(
            IdxBSPolar idx_test,
            IdxQuadratureRTheta idx_quad,
            IdxBatch idx_batch,
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            MetricTensorType const& metric_tensor,
            DField<IdxRangeQuadratureRTheta> int_volume)
    {
        auto [basis_val_test_space, basis_derivs_test_space]
                = detail_poisson::get_polar_bspline_vals_and_derivs(idx_quad, idx_test);
        return detail_poisson::get_test_weights(
                basis_val_test_space,
                basis_derivs_test_space,
                idx_quad,
                idx_batch,
                coeff_alpha,
                coeff_beta,
                metric_tensor,
                int_volume);
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
        Kokkos::deep_copy(nnz_per_row, 0.);
        init_singular_nnz_per_line(nnz_per_row);
        m_tensor_assembler.init_nnz_per_line(get_tensor_nnz_per_row(nnz_per_row));
    }

    /**
     * @brief Fills the nnz data structure for the rows of the singular basis functions and
     * adds the contribution of the singular columns to the rows of the tensor basis functions.
     *
     * The rows of the tensor basis functions near the singular point contain one non-zero
     * element for each singular basis function. As these elements are found at the start of
     * the row, all subsequent rows are shifted by the same amount. This function adds these
     * contributions to the values already saved in nnz_per_row so it can be called before or
     * after TensorProductBSplinePoissonAssembler::init_nnz_per_line.
     *
     * @param[inout] nnz_per_row A 1D Kokkos view of length matrix_size+1 which stores the sum of
     *                  the non-zeros in the matrix on all lines up to the one in.
     */
    void init_singular_nnz_per_line(Kokkos::View<int*, Kokkos::LayoutRight> nnz_per_row) const
    {
        int const mat_size = static_cast<int>(nnz_per_row.extent(0)) - 1;

        int constexpr n_singular_basis = PolarBSplinesRTheta::n_singular_basis();
        int constexpr degree = BSplinesR::degree();
        int const nbasis_theta_proxy = m_nbasis_theta;

        // The number of tensor basis functions which overlap with the singular basis functions
        int const n_overlapping_singular = degree * nbasis_theta_proxy;

        int const nnz_for_singular_rows = n_singular_basis + n_overlapping_singular;

        // rows representing the bsplines which cover the singular domain
        // These overlap with other singular domain splines and degree radial splines
        Kokkos::parallel_for(
                "overlap singular radial",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(1, n_singular_basis + 1),
                KOKKOS_LAMBDA(const int k) { nnz_per_row(k) = k * nnz_for_singular_rows; });

        int const nnz_singular_rows = nnz_for_singular_rows * n_singular_basis;

        // rows representing the bsplines which can be represented as a tensor product of
        // radial and poloidal bsplines. The rows which overlap with the singular domain contain
        // n_singular_basis additional elements.
        Kokkos::parallel_for(
                "singular columns",
                Kokkos::RangePolicy<
                        Kokkos::DefaultExecutionSpace>(n_singular_basis + 1, mat_size + 1),
                KOKKOS_LAMBDA(const int k) {
                    nnz_per_row(k)
                            += nnz_singular_rows
                               + n_singular_basis
                                         * Kokkos::
                                                 min(k - n_singular_basis, n_overlapping_singular);
                });
    }

    /**
     * @brief Convert a 2D (r,theta) bspline index into a polar bspline index.
     *
     * @param[in] idx The 2D (r,theta) bspline index.
     * @return The polar bspline index.
     */
    static KOKKOS_INLINE_FUNCTION IdxBSPolar to_polar(IdxBSRTheta idx)
    {
        return PolarBSplinesRTheta::template get_polar_index<PolarBSplinesRTheta>(idx);
    }

private:
    static std::size_t get_matrix_size(bool include_boundary_splines)
    {
        if constexpr (include_boundary_splines)
            return ddc::discrete_space<PolarBSplinesRTheta>().nbasis()
                   - ddc::discrete_space<BSplinesTheta>().nbasis();
        else
            return ddc::discrete_space<PolarBSplinesRTheta>().nbasis();
    }
};
