// SPDX-License-Identifier: MIT
#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <optional>
#include <source_location>
#include <tuple>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "fem_poisson_assembler_utils.hpp"
#include "math_tools.hpp"
#include "matrix_batch_csr.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"
#include "view.hpp"


/**
 * @brief An operator to assemble a Poisson-like stiffness matrix using tensor product B-splines.
 *
 * The stiffness matrix has entries
 * @f$ \int \alpha \nabla B_i \cdot \nabla B_j + \beta B_i B_j dx @f$
 * where @f$ B_i @f$ and @f$ B_j @f$ are tensor product B-splines.
 *
 * The B-splines and the quadrature points both have a tensor product structure. The values
 * and derivatives of the 1D B-splines at the 1D quadrature points are therefore stored (see
 * detail_poisson::TensorBSplineStash) so that the tensor product B-splines can be evaluated
 * at the quadrature points without being recalculated.
 *
 * @tparam GridR The radial grid type.
 * @tparam GridTheta The poloidal grid type.
 * @tparam BSplinesR The 1D B-splines in the radial direction.
 * @tparam BSplinesTheta The 1D B-splines in the poloidal direction.
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
        typename BSplinesR,
        typename BSplinesTheta,
        typename QDimRMesh,
        typename QDimThetaMesh,
        concepts::MetricTensor MetricTensorType,
        typename IdxRangeFull = IdxRange<GridR, GridTheta>>
class TensorProductBSplinePoissonAssembler
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
    static_assert(
            !R::PERIODIC && Theta::PERIODIC,
            "TensorProductBSplinePoissonAssembler is currently only implemented for a ring "
            "geometry");

    /// The index range of the tensor product B-splines
    using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;

    /// The index range of the batch dimensions
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFull, GridR, GridTheta>;

    /// The index range of the quadrature points
    using IdxRangeQuadratureRTheta = IdxRange<QDimRMesh, QDimThetaMesh>;

    /// The values and derivatives of the 1D B-splines at the quadrature points.
    using TensorBSplineStash = detail_poisson::
            TensorBSplineStash<QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta>;

private:
    using KnotsR = ddc::knot_discrete_dimension_t<BSplinesR>;
    using KnotsTheta = ddc::knot_discrete_dimension_t<BSplinesTheta>;

    using IdxRangeBSR = IdxRange<BSplinesR>;
    using IdxRangeBSTheta = IdxRange<BSplinesTheta>;

    using IdxBSR = Idx<BSplinesR>;
    using IdxBSTheta = Idx<BSplinesTheta>;
    using IdxBSRTheta = Idx<BSplinesR, BSplinesTheta>;

    using IdxStepBSR = IdxStep<BSplinesR>;
    using IdxStepBSTheta = IdxStep<BSplinesTheta>;

    using IdxBatch = typename IdxRangeBatch::discrete_element_type;

    using IdxRangeQuadratureR = IdxRange<QDimRMesh>;
    using IdxRangeQuadratureTheta = IdxRange<QDimThetaMesh>;
    using IdxQuadratureRTheta = Idx<QDimRMesh, QDimThetaMesh>;

    /**
     * @brief The index range of the stash of the non-zero radial B-splines at the radial
     * quadrature points. The B-spline index is the offset from the first non-zero B-spline
     * at the quadrature point.
     */
    using IdxRangeBasisStashR = IdxRange<QDimRMesh, BSplinesR, ddc::Deriv<R>>;
    /**
     * @brief The index range of the stash of the non-zero poloidal B-splines at the poloidal
     * quadrature points. The B-spline index is the offset from the first non-zero B-spline
     * at the quadrature point.
     */
    using IdxRangeBasisStashTheta = IdxRange<QDimThetaMesh, BSplinesTheta, ddc::Deriv<Theta>>;

    using QuadField = DField<IdxRangeQuadratureRTheta>;

    static constexpr int s_n_gauss_legendre_r = BSplinesR::degree() + 1;
    static constexpr int s_n_gauss_legendre_theta = BSplinesTheta::degree() + 1;

    // Number of cells over which a radial B-splines has its support
    // This is the case for b-splines which are not affected by the higher knot multiplicity at the boundary.
    static constexpr int s_n_non_zero_bases_r = BSplinesR::degree() + 1;

    // Number of cells over which a poloidal B-splines has its support
    static constexpr int s_n_non_zero_bases_theta = BSplinesTheta::degree() + 1;

    // Maximum number of quadrature points in the support of a tensor product B-spline
    static constexpr int s_n_quad_support_r = s_n_gauss_legendre_r * s_n_non_zero_bases_r;
    static constexpr int s_n_quad_support_theta
            = s_n_gauss_legendre_theta * s_n_non_zero_bases_theta;

    // A view on team scratch memory
    using ScratchView = Kokkos::View<
            double*,
            Kokkos::LayoutRight,
            Kokkos::DefaultExecutionSpace::scratch_memory_space,
            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    // The tensor basis functions in the FEM basis. The order of the rows and columns of the
    // tensor-tensor matrix is the layout_right order of this index range.
    IdxRangeBSRTheta m_idxrange_tensor_basis;

    IdxRangeQuadratureR m_idxrange_quadrature_r;
    IdxRangeQuadratureTheta m_idxrange_quadrature_theta;
    IdxRangeQuadratureRTheta m_idxrange_quadrature;

    IdxRangeBatch m_idx_range_batch;

    QuadField m_int_volume;

    MetricTensorType m_metric_tensor;

    // Values and derivatives of the non-zero 1D B-splines at the 1D quadrature points
    DFieldMem<IdxRangeBasisStashR> m_basis_vals_r;
    DFieldMem<IdxRangeBasisStashTheta> m_basis_vals_theta;
    FieldMem<IdxBSR, IdxRangeQuadratureR> m_first_basis_r;
    FieldMem<IdxBSTheta, IdxRangeQuadratureTheta> m_first_basis_theta;

    // The nnz data structure of the matrix (see init_tensor_tensor_nnz_per_line)
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m_tensor_nnz_per_row;

public:
    /**
     * @brief Instantiate the assembler operator.
     *
     * @param int_volume The initialised field of Jacobian values of the mapping.
     * @param metric_tensor The metric tensor associated with the mapping from the logical
     *          domain to the physical domain where the equation is defined.
     * @param idx_range_full The full index range of @f$ \phi @f$ including any batch dimensions.
     * @param idx_range_tensor_basis The index range of the tensor product B-splines in the
     *          FEM basis. B-splines which are excluded from this index range (e.g. to impose
     *          boundary conditions) are not included in the matrix.
     */
    explicit TensorProductBSplinePoissonAssembler(
            QuadField int_volume,
            MetricTensorType const& metric_tensor,
            IdxRangeFull idx_range_full,
            IdxRangeBSRTheta idx_range_tensor_basis)
        : m_idxrange_tensor_basis(idx_range_tensor_basis)
        , m_idxrange_quadrature_r(
                  Idx<QDimRMesh>(0),
                  IdxStep<QDimRMesh>(
                          s_n_gauss_legendre_r * ddc::discrete_space<BSplinesR>().ncells()))
        , m_idxrange_quadrature_theta(
                  Idx<QDimThetaMesh>(0),
                  IdxStep<QDimThetaMesh>(
                          s_n_gauss_legendre_theta * ddc::discrete_space<BSplinesTheta>().ncells()))
        , m_idxrange_quadrature(m_idxrange_quadrature_r, m_idxrange_quadrature_theta)
        , m_idx_range_batch(idx_range_full)
        , m_int_volume(int_volume)
        , m_metric_tensor(metric_tensor)
        , m_basis_vals_r(IdxRangeBasisStashR(
                  m_idxrange_quadrature_r,
                  ddc::discrete_space<BSplinesR>().full_domain().take_first(
                          IdxStepBSR(s_n_non_zero_bases_r)),
                  IdxRange<ddc::Deriv<R>>(Idx<ddc::Deriv<R>>(0), IdxStep<ddc::Deriv<R>>(2))))
        , m_basis_vals_theta(IdxRangeBasisStashTheta(
                  m_idxrange_quadrature_theta,
                  ddc::discrete_space<BSplinesTheta>().full_domain().take_first(
                          IdxStepBSTheta(s_n_non_zero_bases_theta)),
                  IdxRange<ddc::Deriv<
                          Theta>>(Idx<ddc::Deriv<Theta>>(0), IdxStep<ddc::Deriv<Theta>>(2))))
        , m_first_basis_r(m_idxrange_quadrature_r)
        , m_first_basis_theta(m_idxrange_quadrature_theta)
        , m_tensor_nnz_per_row(
                  "m_tensor_nnz_per_row "
                  "(TensorProductBSplinePoissonAssembler::TensorProductBSplinePoissonAssembler)",
                  m_idxrange_tensor_basis.size() + 1)
    {
        stash_basis_eval(get_field(m_basis_vals_r), get_field(m_first_basis_r));
        stash_basis_eval(get_field(m_basis_vals_theta), get_field(m_first_basis_theta));
        init_tensor_tensor_nnz_per_line(m_tensor_nnz_per_row);
    }

    /**
     * @brief Get the TensorBSplineStash of the values and derivatives of the 1D B-splines at the quadrature points.
     *
     * @return The stashed values of the 1D B-splines.
     */
    TensorBSplineStash get_tensor_bspline_stash() const
    {
        return TensorBSplineStash {
                get_const_field(m_basis_vals_r),
                get_const_field(m_basis_vals_theta),
                get_const_field(m_first_basis_r),
                get_const_field(m_first_basis_theta)};
    }

    /**
     * @brief Get the index range of the tensor product B-splines in the FEM basis.
     *
     * @return The index range of the tensor product B-splines.
     */
    IdxRangeBSRTheta get_tensor_basis_idx_range() const
    {
        return m_idxrange_tensor_basis;
    }

    /**
     * @brief Get the integral volume associated with each quadrature point.
     *
     * @return The integral volume.
     */
    QuadField get_int_volume() const
    {
        return m_int_volume;
    }

    /**
     * @brief Compute the sparsity pattern for the stiffness matrix.
     *
     * The rows and columns of the matrix are the tensor product B-splines in the FEM basis
     * ordered following the layout_right order of their index range.
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
     */
    void setup_sparse_matrix(
            std::unique_ptr<
                    MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>&
                    gko_matrix,
            std::optional<int> max_iter = std::nullopt,
            std::optional<double> res_tol = std::nullopt,
            std::optional<bool> batch_solver_logger = std::nullopt,
            std::optional<int> preconditioner_max_block_size = std::nullopt,
            std::optional<GkoMatrixResidualNorm> residual_norm = std::nullopt) const
    {
        //CSR data storage
        gko_matrix = std::make_unique<
                MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>(
                m_idx_range_batch.size(),
                m_idxrange_tensor_basis.size(),
                n_tensor_tensor_elements(),
                max_iter,
                res_tol,
                batch_solver_logger,
                preconditioner_max_block_size,
                residual_norm);
        auto [values, col_idx, nnz_per_row] = gko_matrix->get_batch_csr();
        Kokkos::deep_copy(nnz_per_row, 0);
        init_nnz_per_line(nnz_per_row);
        compute_col_idx(col_idx, nnz_per_row, 0);
    }

    /**
     * @brief Assemble the stiffness matrix.
     *
     * @param[out] gko_matrix The pointer to the matrix initialised by setup_sparse_matrix.
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
            CoeffBeta const& coeff_beta) const
    {
        //CSR data storage
        auto [values, col_idx, nnz_per_row] = gko_matrix->get_batch_csr();

        compute_elements(coeff_alpha, coeff_beta, values, col_idx, nnz_per_row, 0);

        gko_matrix->setup_solver();
    }

    /**
     * @brief Fill the stash with the values and derivatives of the non-zero 1D B-splines at
     * the 1D quadrature points.
     *
     * @param[out] basis_vals The values and derivatives of the non-zero 1D B-splines. The
     *                  B-spline index is the offset from the first non-zero B-spline.
     * @param[out] first_basis The first non-zero B-spline at each quadrature point.
     *
     * @tparam QDimMesh The 1D quadrature grid.
     * @tparam BSplines The 1D B-splines.
     * @tparam DerivDim The derivative dimension.
     */
    template <class QDimMesh, class BSplines, class DerivDim>
    static void stash_basis_eval(
            DField<IdxRange<QDimMesh, BSplines, DerivDim>> basis_vals,
            Field<Idx<BSplines>, IdxRange<QDimMesh>> first_basis)
    {
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                get_idx_range(first_basis),
                KOKKOS_LAMBDA(Idx<QDimMesh> idx_quad) {
                    constexpr std::size_t n_non_zero = BSplines::degree() + 1;
                    std::array<double, n_non_zero * 2> data;
                    DSpan2D derivs(data.data(), n_non_zero, 2);
                    first_basis(idx_quad) = ddc::discrete_space<BSplines>().eval_basis_and_n_derivs(
                            derivs,
                            ddc::coordinate(idx_quad),
                            1);
                    DField<IdxRange<BSplines, DerivDim>> basis_vals_at_quad = basis_vals[idx_quad];
                    for (std::size_t i(0); i < n_non_zero; ++i) {
                        for (std::size_t k(0); k < 2; ++k) {
                            basis_vals_at_quad(IdxStep<BSplines>(i), IdxStep<DerivDim>(k))
                                    = derivs(i, k);
                        }
                    }
                });
    }

    /**
     * @brief Get the number of non-zero elements in the matrix corresponding to the inner
     * products of tensor basis functions and tensor basis functions.
     *
     * @return The number of non-zero elements in the tensor-tensor matrix.
     */
    int n_tensor_tensor_elements() const
    {
        const int nbasis_r = m_idxrange_tensor_basis.template extent<BSplinesR>().value();
        const int nbasis_theta = m_idxrange_tensor_basis.template extent<BSplinesTheta>().value();
        const int n_stencil_theta
                = nbasis_theta * std::min(int(1 + 2 * BSplinesTheta::degree()), nbasis_theta);
        const int n_stencil_r = nbasis_r * (1 + 2 * BSplinesR::degree())
                                - (1 + BSplinesR::degree()) * BSplinesR::degree();
        // These elements have a tensor product structure
        return n_stencil_r * n_stencil_theta;
    }

    /**
     * @brief Add the number of non-zero elements per line of the tensor-tensor matrix to
     * a nnz data structure.
     *
     * After this function the array will contain:
     * nnz_per_row[0] = unchanged (0 for a standalone matrix).
     * nnz_per_row[1] += number of non-zero elements in line 0 of the tensor-tensor matrix.
     * ...
     * nnz_per_row[n] += number of non-zero elements in lines 0-(n-1) of the tensor-tensor matrix.
     *
     * @param[inout] nnz_per_row A 1D Kokkos view of length n+1, where n is the number of tensor
     *                  basis functions, containing the elements of the nnz data structure which
     *                  correspond to the rows of the tensor-tensor matrix.
     */
    void init_nnz_per_line(Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
                                   nnz_per_row) const
    {
        assert(nnz_per_row.size() == m_tensor_nnz_per_row.size());
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> tensor_nnz_per_row
                = m_tensor_nnz_per_row;
        const std::source_location location = std::source_location::current();
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, tensor_nnz_per_row.size()),
                KOKKOS_LAMBDA(const int idx) { nnz_per_row(idx) += tensor_nnz_per_row(idx); });
    }

    /**
     * @brief Computes the column indices of the matrix elements corresponding to the
     * inner products of tensor basis functions and tensor basis functions.
     *
     * The tensor-tensor elements are found at the end of each row of the tensor-tensor
     * matrix.
     *
     * @param[out] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component of
     *             the full matrix (only for one matrix).
     * @param[in] nnz_per_row_csr
     *             A 1D Kokkos view of length n+1, where n is the number of tensor basis functions,
     *             containing the elements of the nnz data structure of the full matrix which
     *             correspond to the rows of the tensor-tensor matrix.
     * @param[in] first_col_offset
     *             The index of the first column of the tensor-tensor matrix in the full matrix.
     */
    void compute_col_idx(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    col_idx_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    nnz_per_row_csr,
            int first_col_offset) const
    {
        IdxRangeBSRTheta idx_range_tensor_basis = m_idxrange_tensor_basis;

        IdxRangeBSR idx_range_fem_r(idx_range_tensor_basis);
        IdxRangeBSTheta full_idx_range_theta(idx_range_tensor_basis);

        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                tensor_nnz_per_row
                = m_tensor_nnz_per_row;

        const std::source_location location = std::source_location::current();

        // Calculate the matrix elements following a stencil
        // Teams loop over rows
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        idx_range_tensor_basis.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    // The row of the tensor-tensor matrix
                    const int row_idx = team.league_rank();

                    // Calculate the index of the test b-spline and its radial and poloidal components
                    const IdxBSRTheta idx_test(
                            ddcHelper::get_idx_from_linear_index(idx_range_tensor_basis, row_idx));
                    const IdxBSR idx_test_r(idx_test);
                    const IdxBSTheta idx_test_theta(idx_test);

                    // Calculate the offsets to locate all trial b-splines which are non-zero
                    // somewhere in the support of the test b-spline.
                    IdxStepBSR idx_step_trial_r_offset_min(
                            Kokkos::
                                    max(IdxStepBSR(-BSplinesR::degree()),
                                        idx_range_fem_r.front() - idx_test_r));
                    IdxStepBSR idx_step_trial_r_offset_max(
                            Kokkos::
                                    min(IdxStepBSR(BSplinesR::degree() + 1),
                                        idx_range_fem_r.back() + 1 - idx_test_r));
                    IdxStepBSTheta idx_step_trial_theta_offset_min(-BSplinesTheta::degree());
                    IdxStepBSTheta idx_step_trial_theta_offset_max(BSplinesTheta::degree() + 1);

                    // Calculate the position of the first tensor-tensor element in the row.
                    // These elements are found at the end of the row.
                    int csr_idx = nnz_per_row_csr(row_idx + 1)
                                  - (tensor_nnz_per_row(row_idx + 1) - tensor_nnz_per_row(row_idx));

                    // Loop over the radial offset to a trial b-spline
                    for (IdxStepBSR idx_step_trial_r(idx_step_trial_r_offset_min);
                         idx_step_trial_r < idx_step_trial_r_offset_max;
                         ++idx_step_trial_r) {
                        const IdxBSR idx_trial_r = idx_test_r + idx_step_trial_r;

                        // As theta is periodic the theta component must be split in 2 to
                        // guarantee that col_idx and values (from the CSR format) are
                        // filled in order.
                        // Calculate the start and end of the theta offset
                        IdxBSTheta first_idx_trial_theta = detail_poisson::
                                mod_add(idx_test_theta,
                                        idx_step_trial_theta_offset_min,
                                        full_idx_range_theta);
                        const IdxBSTheta last_idx_trial_theta = detail_poisson::
                                mod_add(idx_test_theta,
                                        idx_step_trial_theta_offset_max,
                                        full_idx_range_theta);
                        IdxBSTheta first_periodic_idx_trial_theta = full_idx_range_theta.back() + 1;
                        // If the start is after the end then the periodicity wraps the b-splines
                        // first_periodic_idx_trial_theta is modified to fill the wrapped area (at the
                        // right hand side of the block) last.
                        if (first_idx_trial_theta > last_idx_trial_theta) {
                            first_periodic_idx_trial_theta = first_idx_trial_theta;
                            first_idx_trial_theta = full_idx_range_theta.front();
                        }
                        // Loop over the poloidal offset to a trial b-spline
                        for (IdxBSTheta idx_trial_theta(first_idx_trial_theta);
                             idx_trial_theta < last_idx_trial_theta;
                             ++idx_trial_theta) {
                            const IdxBSRTheta idx_trial(idx_trial_r, idx_trial_theta);
                            int const col_idx = first_col_offset
                                                + ddcHelper::get_linear_index_from_idx(
                                                        idx_range_tensor_basis,
                                                        idx_trial);
                            col_idx_csr(csr_idx) = col_idx;
                            csr_idx++;
                        }
                        // Loop over the poloidal offset to a trial b-spline for any elements that
                        // are to the right of zeros following the elements already set
                        for (IdxBSTheta idx_trial_theta(first_periodic_idx_trial_theta);
                             idx_trial_theta < full_idx_range_theta.back() + 1;
                             ++idx_trial_theta) {
                            const IdxBSRTheta idx_trial(idx_trial_r, idx_trial_theta);
                            int const col_idx = first_col_offset
                                                + ddcHelper::get_linear_index_from_idx(
                                                        idx_range_tensor_basis,
                                                        idx_trial);
                            col_idx_csr(csr_idx) = col_idx;
                            csr_idx++;
                        }
                    }
                    assert(csr_idx == nnz_per_row_csr(row_idx + 1));
                });
    }

    /**
     * @brief Computes the matrix elements corresponding to the inner products of tensor
     * basis functions and tensor basis functions.
     *
     * @param[in] coeff_alpha
     *      A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning
     *      @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning
     *      @f$ \beta @f$ at the given batch and coordinate.
     * @param[out] values_csr
     *             A 2D Kokkos view which stores the values of non-zero elements of the full matrix
     *             for the whole batch.
     * @param[in] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component of
     *             the full matrix (only for one matrix).
     * @param[in] nnz_per_row_csr
     *             A 1D Kokkos view of length n+1, where n is the number of tensor basis functions,
     *             containing the elements of the nnz data structure of the full matrix which
     *             correspond to the rows of the tensor-tensor matrix.
     * @param[in] first_col_offset
     *             The index of the first column of the tensor-tensor matrix in the full matrix.
     */
    template <class CoeffAlpha, class CoeffBeta>
    void compute_elements(
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    values_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    col_idx_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    nnz_per_row_csr,
            int first_col_offset) const
    {
        IdxRangeBSRTheta idx_range_tensor_basis = m_idxrange_tensor_basis;

        MetricTensorType const metric_tensor = m_metric_tensor;
        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_int_volume;
        TensorBSplineStash bspline_stash = get_tensor_bspline_stash();

        IdxRangeQuadratureRTheta full_quad_idx_range = m_idxrange_quadrature;

        IdxRangeBatch idx_range_batch = m_idx_range_batch;
        const int n_batch = idx_range_batch.size();
        const int n_tensor_basis = idx_range_tensor_basis.size();
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                tensor_nnz_per_row
                = m_tensor_nnz_per_row;
        const std::source_location location = std::source_location::current();

        // Number of quadrature points in the support of a tensor product B-spline
        constexpr int n_quad_support = s_n_quad_support_r * s_n_quad_support_theta;

        // Each team needs scratch memory to store the test weights at the quadrature points
        // in the support of the test function.
        Kokkos::TeamPolicy<>
                policy(Kokkos::DefaultExecutionSpace(), n_batch * n_tensor_basis, Kokkos::AUTO);
        policy.set_scratch_size(0, Kokkos::PerTeam(ScratchView::shmem_size(3 * n_quad_support)));

        // Calculate the matrix elements following a stencil
        // Teams loop over rows
        Kokkos::parallel_for(
                location.function_name(),
                policy,
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const int batch_linear_idx = team.league_rank() / n_tensor_basis;
                    const IdxBatch idx_batch = ddcHelper::
                            get_idx_from_linear_index(idx_range_batch, batch_linear_idx);
                    // The row of the tensor-tensor matrix
                    const int row_idx = team.league_rank() % n_tensor_basis;
                    // Calculate the index of the test b-spline
                    const IdxBSRTheta idx_test(
                            ddcHelper::get_idx_from_linear_index(idx_range_tensor_basis, row_idx));

                    // Calculate the row index and the position of the first tensor-tensor
                    // element in the row. These elements are found at the end of the row.
                    int const csr_idx_first_stencil_element
                            = nnz_per_row_csr[row_idx + 1]
                              - (tensor_nnz_per_row[row_idx + 1] - tensor_nnz_per_row[row_idx]);

                    // Calculate the test weights once for all trial b-splines
                    const IdxRangeQuadratureRTheta idx_range_quad = detail_poisson::
                            get_support_quadrature(idx_test, full_quad_idx_range.front());
                    ScratchView test_weights(team.team_scratch(0), 3 * n_quad_support);
                    QuadField test_mass(test_weights.data(), idx_range_quad);
                    DVectorField<IdxRangeQuadratureRTheta, VectorIndexSet<R, Theta>> test_stiffness(
                            idx_range_quad,
                            test_weights.data() + n_quad_support,
                            test_weights.data() + 2 * n_quad_support);
                    fill_test_weights(
                            team,
                            test_mass,
                            test_stiffness,
                            idx_test,
                            idx_batch,
                            coeff_alpha,
                            coeff_beta,
                            metric_tensor,
                            full_quad_idx_range,
                            int_volume_proxy,
                            bspline_stash);

                    // Loop over the poloidal offset to a trial b-spline
                    for (int csr_idx = csr_idx_first_stencil_element;
                         csr_idx < nnz_per_row_csr[row_idx + 1];
                         ++csr_idx) {
                        int const col_idx = col_idx_csr[csr_idx];
                        const IdxBSRTheta idx_trial(ddcHelper::get_idx_from_linear_index(
                                idx_range_tensor_basis,
                                col_idx - first_col_offset));
                        double element = get_matrix_stencil_element(
                                team,
                                idx_test,
                                idx_trial,
                                test_mass,
                                test_stiffness,
                                full_quad_idx_range,
                                bspline_stash);
                        values_csr(batch_linear_idx, csr_idx) = element;
                    }
                });
    }

    /**
     * @brief Computes the test_stiffness_weight and test_mass_weight for a given quadrature point
     * for a tensor product test basis spline.
     *
     * See detail_poisson::get_test_weights for more details.
     *
     * @param[in] idx_test
     *      The 2D index of the test tensor product B-spline.
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
     * @param[in] bspline_stash
     *      The values and derivatives of the 1D B-splines at the quadrature points.
     * @return A tuple containing test_mass_weight and test_stiffness_weight.
     */
    template <class CoeffAlpha, class CoeffBeta>
    static KOKKOS_FUNCTION std::tuple<double, DVector<R, Theta>> get_test_weights(
            IdxBSRTheta idx_test,
            IdxQuadratureRTheta idx_quad,
            IdxBatch idx_batch,
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            MetricTensorType const& metric_tensor,
            DField<IdxRangeQuadratureRTheta> int_volume,
            TensorBSplineStash const& bspline_stash)
    {
        auto [basis_val_test_space, basis_derivs_test_space]
                = bspline_stash.vals_and_derivs(idx_quad, idx_test);
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
     * @brief Computes the matrix element corresponding to two tensor product splines
     *        with index idx_test and idx_trial
     *
     * @param[in] team
     *      The team of threads from which this function is called.
     * @param[in] idx_test
     *      The 2D index of the tensor product B-spline in the test space.
     * @param[in] idx_trial
     *      The 2D index of the tensor product B-spline in the trial space.
     * @param[in] test_mass
     *      The test_mass_weight at the quadrature points in the support of the test B-spline
     *      (see fill_test_weights).
     * @param[in] test_stiffness
     *      The test_stiffness_weight at the quadrature points in the support of the test
     *      B-spline (see fill_test_weights).
     * @param[in] full_quad_idx_range
     *      The index range of all the quadrature points.
     * @param[in] bspline_stash
     *      The values and derivatives of the 1D B-splines at the quadrature points.
     * @return
     *      The value of the matrix element.
     */
    static KOKKOS_FUNCTION double get_matrix_stencil_element(
            const Kokkos::TeamPolicy<>::member_type& team,
            IdxBSRTheta idx_test,
            IdxBSRTheta idx_trial,
            QuadField test_mass,
            DVectorField<IdxRangeQuadratureRTheta, VectorIndexSet<R, Theta>> test_stiffness,
            IdxRangeQuadratureRTheta const& full_quad_idx_range,
            TensorBSplineStash const& bspline_stash)
    {
        const IdxBSR idx_test_r(idx_test);
        const IdxBSR idx_trial_r(idx_trial);
        const IdxBSTheta idx_test_theta(
                detail_poisson::theta_mod<BSplinesTheta>(IdxBSTheta(idx_test)));
        const IdxBSTheta idx_trial_theta(
                detail_poisson::theta_mod<BSplinesTheta>(IdxBSTheta(idx_trial)));

        auto& bspl_r = ddc::discrete_space<BSplinesR>();
        auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

        const Idx<KnotsR> start_non_zero_r(
                Kokkos::
                        max(bspl_r.break_point_domain().front(),
                            Kokkos::
                                    max(bspl_r.get_first_support_knot(idx_test_r),
                                        bspl_r.get_first_support_knot(idx_trial_r))));
        const Idx<KnotsR> end_non_zero_r(
                Kokkos::
                        min(bspl_r.break_point_domain().back(),
                            Kokkos::
                                    min(bspl_r.get_last_support_knot(idx_test_r),
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
                Kokkos::max(first_support_knot_theta_test, first_support_knot_theta_trial));
        const Idx<KnotsTheta> end_non_zero_theta(
                Kokkos::min(last_support_knot_theta_test, last_support_knot_theta_trial));
        const IdxRangeQuadratureRTheta quad_range = detail_poisson::
                get_quadrature_between_knots<QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta>(
                        start_non_zero_r,
                        end_non_zero_r,
                        start_non_zero_theta,
                        end_non_zero_theta,
                        full_quad_idx_range.front());

        IdxRangeQuadratureRTheta test_support = get_idx_range(test_mass);

        double result = 0;
        Kokkos::parallel_reduce(
                Kokkos::TeamThreadMDRange(
                        team,
                        quad_range.template extent<QDimRMesh>(),
                        quad_range.template extent<QDimThetaMesh>()),
                [&](int r_thread_index, int theta_thread_index, double& sum) {
                    IdxQuadratureRTheta idx_quad
                            = quad_range.front()
                              + IdxStep<
                                      QDimRMesh,
                                      QDimThetaMesh>(r_thread_index, theta_thread_index);
                    // Manage periodicity
                    if (!full_quad_idx_range.contains(idx_quad)) {
                        idx_quad -= full_quad_idx_range.template extent<QDimThetaMesh>();
                    }
                    assert(full_quad_idx_range.contains(idx_quad));

                    auto [basis_val_trial_space, basis_derivs_trial_space]
                            = bspline_stash.vals_and_derivs(idx_quad, idx_trial);

                    // The test support may extend past the end of the poloidal domain
                    IdxQuadratureRTheta idx_quad_test = idx_quad;
                    if (!test_support.contains(idx_quad_test)) {
                        idx_quad_test += full_quad_idx_range.template extent<QDimThetaMesh>();
                    }
                    assert(test_support.contains(idx_quad_test));

                    // Assemble the weak integral element
                    sum += scalar_product(test_stiffness(idx_quad_test), basis_derivs_trial_space)
                           + test_mass(idx_quad_test) * basis_val_trial_space;
                },
                result);
        return result;
    }

    /**
     * @brief Computes the test_mass_weight and test_stiffness_weight at every quadrature point
     *        in the support of a test function and saves them in team scratch memory.
     *
     * @param[in] team
     *      The team of threads from which this function is called.
     * @param[out] test_mass
     *      The field where test_mass_weight is saved. Its index range is the range of
     *      quadrature points in the support of the test B-spline
     *      (see detail_poisson::get_support_quadrature).
     * @param[out] test_stiffness
     *      The field where test_stiffness_weight is saved. It is defined on the same index
     *      range as test_mass.
     * @param[in] idx_test
     *      The 2D index of the test B-spline.
     * @param[in] idx_batch
     *      The index of the batch for which the element is being calculated.
     * @param[in] coeff_alpha
     *      A callable with signature `double operator()(CoordRTheta)` returning @f$ \alpha @f$.
     * @param[in] coeff_beta
     *      A callable with signature `double operator()(CoordRTheta)` returning @f$ \beta @f$.
     * @param[in] metric_tensor
     *      The metric tensor associated with the mapping from the logical domain to the
     *      physical domain where the equation is defined.
     * @param[in] full_quad_idx_range
     *      The index range of all the quadrature points.
     * @param[in] int_volume
     *      The integral volume associated with each point used in the quadrature scheme.
     * @param[in] bspline_stash
     *      The values and derivatives of the 1D B-splines at the quadrature points.
     */
    template <class CoeffAlpha, class CoeffBeta>
    static KOKKOS_FUNCTION void fill_test_weights(
            const Kokkos::TeamPolicy<>::member_type& team,
            QuadField test_mass,
            DVectorField<IdxRangeQuadratureRTheta, VectorIndexSet<R, Theta>> test_stiffness,
            IdxBSRTheta idx_test,
            IdxBatch idx_batch,
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            MetricTensorType const& metric_tensor,
            IdxRangeQuadratureRTheta const& full_quad_idx_range,
            DField<IdxRangeQuadratureRTheta> int_volume,
            TensorBSplineStash const& bspline_stash)
    {
        IdxRangeQuadratureRTheta test_support = get_idx_range(test_mass);
        assert(test_support.template extent<QDimRMesh>().value() <= s_n_quad_support_r);
        assert(test_support.template extent<QDimThetaMesh>().value() <= s_n_quad_support_theta);
        Kokkos::parallel_for(
                Kokkos::TeamThreadMDRange(
                        team,
                        test_support.template extent<QDimRMesh>(),
                        test_support.template extent<QDimThetaMesh>()),
                [&](int r_thread_index, int theta_thread_index) {
                    IdxStep<QDimRMesh, QDimThetaMesh>
                            quad_offset(r_thread_index, theta_thread_index);
                    IdxQuadratureRTheta idx_quad = test_support.front() + quad_offset;
                    IdxQuadratureRTheta idx_quad_in_domain = idx_quad;
                    // Manage periodicity
                    if (!full_quad_idx_range.contains(idx_quad_in_domain)) {
                        idx_quad_in_domain -= full_quad_idx_range.template extent<QDimThetaMesh>();
                    }
                    assert(full_quad_idx_range.contains(idx_quad_in_domain));
                    auto [test_mass_weight, test_stiffness_weight] = get_test_weights(
                            idx_test,
                            idx_quad_in_domain,
                            idx_batch,
                            coeff_alpha,
                            coeff_beta,
                            metric_tensor,
                            int_volume,
                            bspline_stash);
                    test_mass(quad_offset) = test_mass_weight;
                    ddcHelper::assign_vector_field_element(
                            test_stiffness,
                            idx_quad,
                            test_stiffness_weight);
                });
        team.team_barrier();
    }

    /**
     * @brief Fills the nnz data structure of the matrix corresponding to the inner
     * products of tensor basis functions and tensor basis functions.
     *
     * The rows and columns of the matrix are the tensor basis functions, ordered as in
     * the 2D (r,theta) index range of the tensor basis functions. After this function the
     * array will contain:
     * nnz_per_row[0] = unchanged (0 for a standalone matrix).
     * nnz_per_row[1] = number of non-zero elements in line 0 of the matrix.
     * ...
     * nnz_per_row[n] = number of non-zero elements in lines 0-(n-1) of the matrix.
     * The values are added to the values already saved in nnz_per_row.
     *
     * @param[inout] nnz_per_row A 1D Kokkos view of length n+1, where n is the number of tensor
     *                  basis functions, which stores the sum of the non-zeros in the matrix
     *                  on all lines up to the one in.
     */
    void init_tensor_tensor_nnz_per_line(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> nnz_per_row)
            const
    {
        IdxRangeBSRTheta idxrange_tensor_basis = m_idxrange_tensor_basis;
        size_t const mat_size = nnz_per_row.extent(0) - 1;
        assert(mat_size == idxrange_tensor_basis.size());
        // nnz per line
        Field<int, IdxRangeBSRTheta> nnz(nnz_per_row.data() + 1, idxrange_tensor_basis);

        size_t constexpr degree = BSplinesR::degree();
        size_t const nbasis_theta_proxy
                = idxrange_tensor_basis.template extent<BSplinesTheta>().value();
        // The number of radial splines which overlap with a radial spline far from the boundaries
        size_t constexpr stencil_overlap_r = 2 * degree + 1;
        // The number of poloidal splines which overlap with a poloidal spline
        size_t const stencil_overlap_theta
                = std::min(size_t(2 * BSplinesTheta::degree() + 1), nbasis_theta_proxy);

        // The number of radial splines which overlap with a given radial spline to its
        // left or to its right
        IdxStepBSR n_one_side_overlap(BSplinesR::degree());

        int nnz_sum = 0;

        IdxRangeBSR idxrange_bsplines_r(idxrange_tensor_basis);
        IdxRangeBSTheta idxrange_bsplines_theta(idxrange_tensor_basis);
        IdxRangeBSR inner_bsplines_r = idxrange_bsplines_r.take_first(n_one_side_overlap);
        IdxRangeBSRTheta inner_bsplines_2d(inner_bsplines_r, idxrange_bsplines_theta);
        // Approaching the internal boundary the overlapping possibilities between two radial splines decrease
        Kokkos::parallel_for(
                "inner boundary",
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        inner_bsplines_r.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    IdxStepBSR ir(team.league_rank());
                    // number of radial bases which overlap with this element
                    int nr = degree + ir + 1;

                    // nnz to this line = ntheta * n_stencil_theta \sum_{k=degree+1}^{nr-1} k
                    int nnz_sum_to_r = (nr * (nr - 1) - degree * (degree + 1)) * nbasis_theta_proxy
                                       * stencil_overlap_theta / 2;

                    // Loop over poloidal dimensions
                    Kokkos::parallel_for(
                            Kokkos::TeamThreadRange(team, nbasis_theta_proxy),
                            [&](int const& thread_index) {
                                IdxStepBSTheta itheta(thread_index);
                                nnz(ir, itheta)
                                        += nnz_sum_to_r
                                           + (thread_index + 1) * nr * stencil_overlap_theta;
                            });
                });

        {
            int nr = inner_bsplines_r.size();
            nnz_sum += nr * (nr + 2 * degree + 1) * stencil_overlap_theta / 2 * nbasis_theta_proxy;
        }

        IdxRangeBSR idxrange_bsplines_r_stencil
                = idxrange_bsplines_r.remove(n_one_side_overlap, n_one_side_overlap);
        IdxRangeBSRTheta idxrange_stencil(idxrange_bsplines_r_stencil, idxrange_bsplines_theta);
        // Stencil for tensor product bsplines which are not affected by the radial boundaries
        Kokkos::parallel_for(
                "stencil",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, idxrange_stencil.size()),
                KOKKOS_LAMBDA(const int k) {
                    IdxBSRTheta k_2d = idxrange_stencil.front() + IdxStepBSR(k / nbasis_theta_proxy)
                                       + IdxStepBSTheta(k % nbasis_theta_proxy);
                    nnz(k_2d) += nnz_sum + (k + 1) * stencil_overlap_r * stencil_overlap_theta;
                });

        nnz_sum += idxrange_stencil.size() * stencil_overlap_r * stencil_overlap_theta;

        IdxRangeBSR outer_bsplines_r = idxrange_bsplines_r.take_last(n_one_side_overlap);
        IdxRangeBSRTheta outer_bsplines_2d(outer_bsplines_r, idxrange_bsplines_theta);
        // Approaching the external boundary the overlapping possibilities between two radial splines decrease
        Kokkos::parallel_for(
                "outer boundary",
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        outer_bsplines_r.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    IdxStepBSR ir(team.league_rank());
                    int nr = 2 * degree - ir;
                    int nnz_sum_to_r = nnz_sum
                                       + (degree * (2 * degree + 1) - nr * (nr + 1) / 2)
                                                 * stencil_overlap_theta * nbasis_theta_proxy;
                    // Loop over poloidal dimensions
                    Kokkos::parallel_for(
                            Kokkos::TeamThreadRange(team, nbasis_theta_proxy),
                            [&](int const& thread_index) {
                                IdxStepBSTheta itheta(thread_index);
                                IdxBSRTheta k = outer_bsplines_2d.front() + ir + itheta;
                                nnz(k) += nnz_sum_to_r
                                          + (thread_index + 1) * nr * stencil_overlap_theta;
                            });
                });
    }
};
