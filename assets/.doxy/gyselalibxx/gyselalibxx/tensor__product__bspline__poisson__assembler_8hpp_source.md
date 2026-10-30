

# File tensor\_product\_bspline\_poisson\_assembler.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**tensor\_product\_bspline\_poisson\_assembler.hpp**](tensor__product__bspline__poisson__assembler_8hpp.md)

[Go to the documentation of this file](tensor__product__bspline__poisson__assembler_8hpp.md)


```C++
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
    using R = typename GridR::continuous_dimension_type;
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

    using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFull, GridR, GridTheta>;

    using IdxRangeQuadratureRTheta = IdxRange<QDimRMesh, QDimThetaMesh>;

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

    using IdxRangeBasisStashR = IdxRange<QDimRMesh, BSplinesR, ddc::Deriv<R>>;
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

    TensorBSplineStash get_tensor_bspline_stash() const
    {
        return TensorBSplineStash {
                get_const_field(m_basis_vals_r),
                get_const_field(m_basis_vals_theta),
                get_const_field(m_first_basis_r),
                get_const_field(m_first_basis_theta)};
    }

    IdxRangeBSRTheta get_tensor_basis_idx_range() const
    {
        return m_idxrange_tensor_basis;
    }

    QuadField get_int_volume() const
    {
        return m_int_volume;
    }

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

    void init_tensor_tensor_nnz_per_line(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> nnz_per_row)
            const
    {
        IdxRangeBSRTheta idxrange_tensor_basis = m_idxrange_tensor_basis;
        assert(nnz_per_row.extent(0) == idxrange_tensor_basis.size() + 1);
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
```


