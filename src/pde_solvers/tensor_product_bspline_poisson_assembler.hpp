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
 * @tparam Grid1 The grid type in the first dimension.
 * @tparam Grid2 The grid type in the second dimension.
 * @tparam BSplines1 The 1D B-splines in the first dimension.
 * @tparam BSplines2 The 1D B-splines in the second dimension.
 * @tparam QDim1Mesh The quadrature grid type in the first dimension.
 * @tparam QDim2Mesh The quadrature grid type in the second dimension.
 * @tparam MetricTensorType The type of the metric tensor associated with the mapping from the
 *          logical domain to the physical domain where the equation is defined.
 * @tparam IdxRangeFull The full index range of @f$ \phi @f$ including any batch dimensions.
 *          One matrix is assembled for each index in the batch dimensions.
 */
template <
        typename Grid1,
        typename Grid2,
        typename BSplines1,
        typename BSplines2,
        typename QDim1Mesh,
        typename QDim2Mesh,
        concepts::MetricTensor MetricTensorType,
        typename IdxRangeFull = IdxRange<Grid1, Grid2>>
class TensorProductBSplinePoissonAssembler
{
public:
    /// The first dimension
    using Dim1 = typename Grid1::continuous_dimension_type;
    /// The second dimension
    using Dim2 = typename Grid2::continuous_dimension_type;

    static_assert(Dim1::IS_CONTRAVARIANT);
    static_assert(Dim2::IS_CONTRAVARIANT);
    static_assert(
            std::is_same_v<typename MetricTensorType::CoordArg, Coord<Dim1, Dim2>>,
            "The metric tensor must be evaluated on (Dim1, Dim2) coordinates.");
    static_assert(
            !Dim1::PERIODIC && Dim2::PERIODIC,
            "TensorProductBSplinePoissonAssembler is currently only implemented for a ring "
            "geometry");

    /// The index range of the tensor product B-splines
    using IdxRangeBS2D = IdxRange<BSplines1, BSplines2>;

    /// The index range of the batch dimensions
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFull, Grid1, Grid2>;

    /// The index range of the quadrature points
    using IdxRangeQuadrature2D = IdxRange<QDim1Mesh, QDim2Mesh>;

    /// The values and derivatives of the 1D B-splines at the quadrature points.
    using TensorBSplineStash
            = detail_poisson::TensorBSplineStash<QDim1Mesh, QDim2Mesh, BSplines1, BSplines2>;

private:
    using Knots1 = ddc::knot_discrete_dimension_t<BSplines1>;
    using Knots2 = ddc::knot_discrete_dimension_t<BSplines2>;

    using IdxRangeBS1 = IdxRange<BSplines1>;
    using IdxRangeBS2 = IdxRange<BSplines2>;

    using IdxBS1 = Idx<BSplines1>;
    using IdxBS2 = Idx<BSplines2>;
    using IdxBS2D = Idx<BSplines1, BSplines2>;

    using IdxStepBS1 = IdxStep<BSplines1>;
    using IdxStepBS2 = IdxStep<BSplines2>;

    using IdxBatch = typename IdxRangeBatch::discrete_element_type;

    using IdxRangeQuadrature1 = IdxRange<QDim1Mesh>;
    using IdxRangeQuadrature2 = IdxRange<QDim2Mesh>;
    using IdxQuadrature2D = Idx<QDim1Mesh, QDim2Mesh>;

    /**
     * @brief The index range of the stash of the non-zero B-splines at the quadrature points
     * in the first dimension. The B-spline index is the offset from the first non-zero B-spline
     * at the quadrature point.
     */
    using IdxRangeBasisStash1 = IdxRange<QDim1Mesh, BSplines1, ddc::Deriv<Dim1>>;
    /**
     * @brief The index range of the stash of the non-zero B-splines at the quadrature points
     * in the second dimension. The B-spline index is the offset from the first non-zero B-spline
     * at the quadrature point.
     */
    using IdxRangeBasisStash2 = IdxRange<QDim2Mesh, BSplines2, ddc::Deriv<Dim2>>;

    using QuadField = DField<IdxRangeQuadrature2D>;

    static constexpr int s_n_gauss_legendre_1 = BSplines1::degree() + 1;
    static constexpr int s_n_gauss_legendre_2 = BSplines2::degree() + 1;

    // Number of cells over which a B-spline in the first dimension has its support
    // This is the case for b-splines which are not affected by the higher knot multiplicity at the boundary.
    static constexpr int s_n_non_zero_bases_1 = BSplines1::degree() + 1;

    // Number of cells over which a B-spline in the second dimension has its support
    static constexpr int s_n_non_zero_bases_2 = BSplines2::degree() + 1;

    // Maximum number of quadrature points in the support of a tensor product B-spline
    static constexpr int s_n_quad_support_1 = s_n_gauss_legendre_1 * s_n_non_zero_bases_1;
    static constexpr int s_n_quad_support_2 = s_n_gauss_legendre_2 * s_n_non_zero_bases_2;

    // A view on team scratch memory
    using ScratchView = Kokkos::View<
            double*,
            Kokkos::LayoutRight,
            Kokkos::DefaultExecutionSpace::scratch_memory_space,
            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    // The tensor basis functions in the FEM basis. The order of the rows and columns of the
    // tensor-tensor matrix is the layout_right order of this index range.
    IdxRangeBS2D m_idxrange_tensor_basis;

    IdxRangeQuadrature1 m_idxrange_quadrature_1;
    IdxRangeQuadrature2 m_idxrange_quadrature_2;
    IdxRangeQuadrature2D m_idxrange_quadrature;

    IdxRangeBatch m_idx_range_batch;

    QuadField m_int_volume;

    MetricTensorType m_metric_tensor;

    // Values and derivatives of the non-zero 1D B-splines at the 1D quadrature points
    DFieldMem<IdxRangeBasisStash1> m_basis_vals_1;
    DFieldMem<IdxRangeBasisStash2> m_basis_vals_2;
    FieldMem<IdxBS1, IdxRangeQuadrature1> m_first_basis_1;
    FieldMem<IdxBS2, IdxRangeQuadrature2> m_first_basis_2;

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
            IdxRangeBS2D idx_range_tensor_basis)
        : m_idxrange_tensor_basis(idx_range_tensor_basis)
        , m_idxrange_quadrature_1(
                  Idx<QDim1Mesh>(0),
                  IdxStep<QDim1Mesh>(
                          s_n_gauss_legendre_1 * ddc::discrete_space<BSplines1>().ncells()))
        , m_idxrange_quadrature_2(
                  Idx<QDim2Mesh>(0),
                  IdxStep<QDim2Mesh>(
                          s_n_gauss_legendre_2 * ddc::discrete_space<BSplines2>().ncells()))
        , m_idxrange_quadrature(m_idxrange_quadrature_1, m_idxrange_quadrature_2)
        , m_idx_range_batch(idx_range_full)
        , m_int_volume(int_volume)
        , m_metric_tensor(metric_tensor)
        , m_basis_vals_1(IdxRangeBasisStash1(
                  m_idxrange_quadrature_1,
                  ddc::discrete_space<BSplines1>().full_domain().take_first(
                          IdxStepBS1(s_n_non_zero_bases_1)),
                  IdxRange<ddc::Deriv<
                          Dim1>>(Idx<ddc::Deriv<Dim1>>(0), IdxStep<ddc::Deriv<Dim1>>(2))))
        , m_basis_vals_2(IdxRangeBasisStash2(
                  m_idxrange_quadrature_2,
                  ddc::discrete_space<BSplines2>().full_domain().take_first(
                          IdxStepBS2(s_n_non_zero_bases_2)),
                  IdxRange<ddc::Deriv<
                          Dim2>>(Idx<ddc::Deriv<Dim2>>(0), IdxStep<ddc::Deriv<Dim2>>(2))))
        , m_first_basis_1(m_idxrange_quadrature_1)
        , m_first_basis_2(m_idxrange_quadrature_2)
        , m_tensor_nnz_per_row(
                  "m_tensor_nnz_per_row "
                  "(TensorProductBSplinePoissonAssembler::TensorProductBSplinePoissonAssembler)",
                  m_idxrange_tensor_basis.size() + 1)
    {
        stash_basis_eval(get_field(m_basis_vals_1), get_field(m_first_basis_1));
        stash_basis_eval(get_field(m_basis_vals_2), get_field(m_first_basis_2));
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
                get_const_field(m_basis_vals_1),
                get_const_field(m_basis_vals_2),
                get_const_field(m_first_basis_1),
                get_const_field(m_first_basis_2)};
    }

    /**
     * @brief Get the index range of the tensor product B-splines in the FEM basis.
     *
     * @return The index range of the tensor product B-splines.
     */
    IdxRangeBS2D get_tensor_basis_idx_range() const
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
     *      A callable object with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)`
     *      returning the value of @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable object with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)`
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
        const int nbasis_1 = m_idxrange_tensor_basis.template extent<BSplines1>().value();
        const int nbasis_2 = m_idxrange_tensor_basis.template extent<BSplines2>().value();
        const int n_stencil_2 = nbasis_2 * std::min(int(1 + 2 * BSplines2::degree()), nbasis_2);
        const int n_stencil_1 = nbasis_1 * (1 + 2 * BSplines1::degree())
                                - (1 + BSplines1::degree()) * BSplines1::degree();
        // These elements have a tensor product structure
        return n_stencil_1 * n_stencil_2;
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
        IdxRangeBS2D idx_range_tensor_basis = m_idxrange_tensor_basis;

        IdxRangeBS1 idx_range_fem_1(idx_range_tensor_basis);
        IdxRangeBS2 full_idx_range_2(idx_range_tensor_basis);

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

                    // Calculate the index of the test b-spline and its components in each dimension
                    const IdxBS2D idx_test(
                            ddcHelper::get_idx_from_linear_index(idx_range_tensor_basis, row_idx));
                    const IdxBS1 idx_test_1(idx_test);
                    const IdxBS2 idx_test_2(idx_test);

                    // Calculate the offsets to locate all trial b-splines which are non-zero
                    // somewhere in the support of the test b-spline.
                    IdxStepBS1 idx_step_trial_1_offset_min(
                            Kokkos::
                                    max(IdxStepBS1(-BSplines1::degree()),
                                        idx_range_fem_1.front() - idx_test_1));
                    IdxStepBS1 idx_step_trial_1_offset_max(
                            Kokkos::
                                    min(IdxStepBS1(BSplines1::degree() + 1),
                                        idx_range_fem_1.back() + 1 - idx_test_1));
                    IdxStepBS2 idx_step_trial_2_offset_min(-BSplines2::degree());
                    IdxStepBS2 idx_step_trial_2_offset_max(BSplines2::degree() + 1);

                    // Calculate the position of the first tensor-tensor element in the row.
                    // These elements are found at the end of the row.
                    int csr_idx = nnz_per_row_csr(row_idx + 1)
                                  - (tensor_nnz_per_row(row_idx + 1) - tensor_nnz_per_row(row_idx));

                    // Loop over the offset to a trial b-spline in the first dimension
                    for (IdxStepBS1 idx_step_trial_1(idx_step_trial_1_offset_min);
                         idx_step_trial_1 < idx_step_trial_1_offset_max;
                         ++idx_step_trial_1) {
                        const IdxBS1 idx_trial_1 = idx_test_1 + idx_step_trial_1;

                        // As the second dimension is periodic its component must be split in 2
                        // to guarantee that col_idx and values (from the CSR format) are
                        // filled in order.
                        // Calculate the start and end of the offset in the second dimension
                        IdxBS2 first_idx_trial_2 = detail_poisson::
                                mod_add(idx_test_2, idx_step_trial_2_offset_min, full_idx_range_2);
                        const IdxBS2 last_idx_trial_2 = detail_poisson::
                                mod_add(idx_test_2, idx_step_trial_2_offset_max, full_idx_range_2);
                        IdxBS2 first_periodic_idx_trial_2 = full_idx_range_2.back() + 1;
                        // If the start is after the end then the periodicity wraps the b-splines
                        // first_periodic_idx_trial_2 is modified to fill the wrapped area (at the
                        // right hand side of the block) last.
                        if (first_idx_trial_2 > last_idx_trial_2) {
                            first_periodic_idx_trial_2 = first_idx_trial_2;
                            first_idx_trial_2 = full_idx_range_2.front();
                        }
                        // Loop over the offset to a trial b-spline in the second dimension
                        for (IdxBS2 idx_trial_2(first_idx_trial_2); idx_trial_2 < last_idx_trial_2;
                             ++idx_trial_2) {
                            const IdxBS2D idx_trial(idx_trial_1, idx_trial_2);
                            int const col_idx = first_col_offset
                                                + ddcHelper::get_linear_index_from_idx(
                                                        idx_range_tensor_basis,
                                                        idx_trial);
                            col_idx_csr(csr_idx) = col_idx;
                            csr_idx++;
                        }
                        // Loop over the offset to a trial b-spline in the second dimension for any elements that
                        // are to the right of zeros following the elements already set
                        for (IdxBS2 idx_trial_2(first_periodic_idx_trial_2);
                             idx_trial_2 < full_idx_range_2.back() + 1;
                             ++idx_trial_2) {
                            const IdxBS2D idx_trial(idx_trial_1, idx_trial_2);
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
     *      A callable with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)` returning
     *      @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)` returning
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
        IdxRangeBS2D idx_range_tensor_basis = m_idxrange_tensor_basis;

        MetricTensorType const metric_tensor = m_metric_tensor;
        DField<IdxRangeQuadrature2D> int_volume_proxy = m_int_volume;
        TensorBSplineStash bspline_stash = get_tensor_bspline_stash();

        IdxRangeQuadrature2D full_quad_idx_range = m_idxrange_quadrature;

        IdxRangeBatch idx_range_batch = m_idx_range_batch;
        const int n_batch = idx_range_batch.size();
        const int n_tensor_basis = idx_range_tensor_basis.size();
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                tensor_nnz_per_row
                = m_tensor_nnz_per_row;
        const std::source_location location = std::source_location::current();

        // Number of quadrature points in the support of a tensor product B-spline
        constexpr int n_quad_support = s_n_quad_support_1 * s_n_quad_support_2;

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
                    const IdxBS2D idx_test(
                            ddcHelper::get_idx_from_linear_index(idx_range_tensor_basis, row_idx));

                    // Calculate the row index and the position of the first tensor-tensor
                    // element in the row. These elements are found at the end of the row.
                    int const csr_idx_first_stencil_element
                            = nnz_per_row_csr[row_idx + 1]
                              - (tensor_nnz_per_row[row_idx + 1] - tensor_nnz_per_row[row_idx]);

                    // Calculate the test weights once for all trial b-splines
                    const IdxRangeQuadrature2D idx_range_quad = detail_poisson::
                            get_support_quadrature(idx_test, full_quad_idx_range.front());
                    ScratchView test_weights(team.team_scratch(0), 3 * n_quad_support);
                    QuadField test_mass(test_weights.data(), idx_range_quad);
                    DVectorField<IdxRangeQuadrature2D, VectorIndexSet<Dim1, Dim2>> test_stiffness(
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

                    // Loop over the trial b-splines
                    for (int csr_idx = csr_idx_first_stencil_element;
                         csr_idx < nnz_per_row_csr[row_idx + 1];
                         ++csr_idx) {
                        int const col_idx = col_idx_csr[csr_idx];
                        const IdxBS2D idx_trial(ddcHelper::get_idx_from_linear_index(
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
     *      A callable with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)` returning
     *      @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)` returning
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
    static KOKKOS_FUNCTION std::tuple<double, DVector<Dim1, Dim2>> get_test_weights(
            IdxBS2D idx_test,
            IdxQuadrature2D idx_quad,
            IdxBatch idx_batch,
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            MetricTensorType const& metric_tensor,
            DField<IdxRangeQuadrature2D> int_volume,
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
            IdxBS2D idx_test,
            IdxBS2D idx_trial,
            QuadField test_mass,
            DVectorField<IdxRangeQuadrature2D, VectorIndexSet<Dim1, Dim2>> test_stiffness,
            IdxRangeQuadrature2D const& full_quad_idx_range,
            TensorBSplineStash const& bspline_stash)
    {
        // Find the cells where both the test and trial functions are non-zero
        const IdxRange<Knots1> overlap_knots_1 = detail_poisson::get_overlapping_knots_1d<
                BSplines1>(IdxBS1(idx_test), IdxBS1(idx_trial));
        const IdxRange<Knots2> overlap_knots_2 = detail_poisson::get_overlapping_knots_1d<
                BSplines2>(IdxBS2(idx_test), IdxBS2(idx_trial));
        const IdxRangeQuadrature2D quad_range = detail_poisson::
                get_quadrature_between_knots<QDim1Mesh, QDim2Mesh, BSplines1, BSplines2>(
                        overlap_knots_1.front(),
                        overlap_knots_1.front() + overlap_knots_1.extents(),
                        overlap_knots_2.front(),
                        overlap_knots_2.front() + overlap_knots_2.extents(),
                        full_quad_idx_range.front());

        IdxRangeQuadrature2D test_support = get_idx_range(test_mass);

        double result = 0;
        Kokkos::parallel_reduce(
                Kokkos::TeamThreadMDRange(
                        team,
                        quad_range.template extent<QDim1Mesh>(),
                        quad_range.template extent<QDim2Mesh>()),
                [&](int thread_index_1, int thread_index_2, double& sum) {
                    IdxQuadrature2D idx_quad
                            = quad_range.front()
                              + IdxStep<QDim1Mesh, QDim2Mesh>(thread_index_1, thread_index_2);
                    // Manage periodicity
                    idx_quad = detail_poisson::periodic_wrap<
                            BSplines1,
                            BSplines2>(idx_quad, full_quad_idx_range, full_quad_idx_range);

                    auto [basis_val_trial_space, basis_derivs_trial_space]
                            = bspline_stash.vals_and_derivs(idx_quad, idx_trial);

                    // The test support may extend past the end of the domain in periodic dimensions
                    IdxQuadrature2D idx_quad_test = detail_poisson::periodic_wrap<
                            BSplines1,
                            BSplines2>(idx_quad, test_support, full_quad_idx_range);

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
     *      A callable with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)` returning
     *      @f$ \alpha @f$ at the given batch and coordinate.
     * @param[in] coeff_beta
     *      A callable with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)` returning
     *      @f$ \beta @f$ at the given batch and coordinate.
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
            DVectorField<IdxRangeQuadrature2D, VectorIndexSet<Dim1, Dim2>> test_stiffness,
            IdxBS2D idx_test,
            IdxBatch idx_batch,
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            MetricTensorType const& metric_tensor,
            IdxRangeQuadrature2D const& full_quad_idx_range,
            DField<IdxRangeQuadrature2D> int_volume,
            TensorBSplineStash const& bspline_stash)
    {
        IdxRangeQuadrature2D test_support = get_idx_range(test_mass);
        assert(test_support.template extent<QDim1Mesh>().value() <= s_n_quad_support_1);
        assert(test_support.template extent<QDim2Mesh>().value() <= s_n_quad_support_2);
        Kokkos::parallel_for(
                Kokkos::TeamThreadMDRange(
                        team,
                        test_support.template extent<QDim1Mesh>(),
                        test_support.template extent<QDim2Mesh>()),
                [&](int thread_index_1, int thread_index_2) {
                    IdxStep<QDim1Mesh, QDim2Mesh> quad_offset(thread_index_1, thread_index_2);
                    IdxQuadrature2D idx_quad = test_support.front() + quad_offset;
                    // Manage periodicity
                    IdxQuadrature2D idx_quad_in_domain = detail_poisson::periodic_wrap<
                            BSplines1,
                            BSplines2>(idx_quad, full_quad_idx_range, full_quad_idx_range);
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
     * the 2D index range of the tensor basis functions. After this function the
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
        IdxRangeBS2D idxrange_tensor_basis = m_idxrange_tensor_basis;
        assert(nnz_per_row.extent(0) == idxrange_tensor_basis.size() + 1);
        // nnz per line
        Field<int, IdxRangeBS2D> nnz(nnz_per_row.data() + 1, idxrange_tensor_basis);

        size_t constexpr degree_1 = BSplines1::degree();
        size_t const nbasis_2_proxy = idxrange_tensor_basis.template extent<BSplines2>().value();
        // The number of splines in the first dimension which overlap with a spline far from the boundaries
        size_t constexpr stencil_overlap_1 = 2 * degree_1 + 1;
        // The number of splines in the second dimension which overlap with a spline
        size_t const stencil_overlap_2
                = std::min(size_t(2 * BSplines2::degree() + 1), nbasis_2_proxy);

        // The number of splines in the first dimension which overlap with a given spline to its
        // left or to its right
        IdxStepBS1 n_one_side_overlap(BSplines1::degree());

        int nnz_sum = 0;

        IdxRangeBS1 idxrange_bsplines_1(idxrange_tensor_basis);
        IdxRangeBS2 idxrange_bsplines_2(idxrange_tensor_basis);
        IdxRangeBS1 inner_bsplines_1 = idxrange_bsplines_1.take_first(n_one_side_overlap);
        IdxRangeBS2D inner_bsplines_2d(inner_bsplines_1, idxrange_bsplines_2);
        // Approaching the lower boundary of the first dimension the overlapping possibilities
        // between two splines decrease
        Kokkos::parallel_for(
                "inner boundary",
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        inner_bsplines_1.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    IdxStepBS1 i1(team.league_rank());
                    // number of bases in the first dimension which overlap with this element
                    int n1 = degree_1 + i1 + 1;

                    // nnz to this line = nbasis_2 * stencil_overlap_2 \sum_{k=degree_1+1}^{n1-1} k
                    int nnz_sum_to_i1 = (n1 * (n1 - 1) - degree_1 * (degree_1 + 1)) * nbasis_2_proxy
                                        * stencil_overlap_2 / 2;

                    // Loop over the second dimension
                    Kokkos::parallel_for(
                            Kokkos::TeamThreadRange(team, nbasis_2_proxy),
                            [&](int const& thread_index) {
                                IdxStepBS2 i2(thread_index);
                                nnz(i1, i2) += nnz_sum_to_i1
                                               + (thread_index + 1) * n1 * stencil_overlap_2;
                            });
                });

        {
            int n1 = inner_bsplines_1.size();
            nnz_sum += n1 * (n1 + 2 * degree_1 + 1) * stencil_overlap_2 / 2 * nbasis_2_proxy;
        }

        IdxRangeBS1 idxrange_bsplines_1_stencil
                = idxrange_bsplines_1.remove(n_one_side_overlap, n_one_side_overlap);
        IdxRangeBS2D idxrange_stencil(idxrange_bsplines_1_stencil, idxrange_bsplines_2);
        // Stencil for tensor product bsplines which are not affected by the boundaries of the
        // first dimension
        Kokkos::parallel_for(
                "stencil",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, idxrange_stencil.size()),
                KOKKOS_LAMBDA(const int k) {
                    IdxBS2D k_2d = idxrange_stencil.front() + IdxStepBS1(k / nbasis_2_proxy)
                                   + IdxStepBS2(k % nbasis_2_proxy);
                    nnz(k_2d) += nnz_sum + (k + 1) * stencil_overlap_1 * stencil_overlap_2;
                });

        nnz_sum += idxrange_stencil.size() * stencil_overlap_1 * stencil_overlap_2;

        IdxRangeBS1 outer_bsplines_1 = idxrange_bsplines_1.take_last(n_one_side_overlap);
        IdxRangeBS2D outer_bsplines_2d(outer_bsplines_1, idxrange_bsplines_2);
        // Approaching the upper boundary of the first dimension the overlapping possibilities
        // between two splines decrease
        Kokkos::parallel_for(
                "outer boundary",
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        outer_bsplines_1.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    IdxStepBS1 i1(team.league_rank());
                    int n1 = 2 * degree_1 - i1;
                    int nnz_sum_to_i1 = nnz_sum
                                        + (degree_1 * (2 * degree_1 + 1) - n1 * (n1 + 1) / 2)
                                                  * stencil_overlap_2 * nbasis_2_proxy;
                    // Loop over the second dimension
                    Kokkos::parallel_for(
                            Kokkos::TeamThreadRange(team, nbasis_2_proxy),
                            [&](int const& thread_index) {
                                IdxStepBS2 i2(thread_index);
                                IdxBS2D k = outer_bsplines_2d.front() + i1 + i2;
                                nnz(k) += nnz_sum_to_i1
                                          + (thread_index + 1) * n1 * stencil_overlap_2;
                            });
                });
    }
};
