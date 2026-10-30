

# File polar\_spline\_fem\_poisson\_like\_assembler.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**polar\_spline\_fem\_poisson\_like\_assembler.hpp**](polar__spline__fem__poisson__like__assembler_8hpp.md)

[Go to the documentation of this file](polar__spline__fem__poisson__like__assembler_8hpp.md)


```C++
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
    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;

    static_assert(R::IS_CONTRAVARIANT);
    static_assert(Theta::IS_CONTRAVARIANT);
    static_assert(
            std::is_same_v<typename MetricTensorType::CoordArg, Coord<R, Theta>>,
            "The metric tensor must be evaluated on (r, theta) coordinates.");

private:
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

private:
    using BSplinesR = typename PolarBSplinesRTheta::BSplinesR_tag;
    using BSplinesTheta = typename PolarBSplinesRTheta::BSplinesTheta_tag;

    using KnotsR = ddc::knot_discrete_dimension_t<BSplinesR>;
    using KnotsTheta = ddc::knot_discrete_dimension_t<BSplinesTheta>;

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

    using IdxRangeQuadratureR = IdxRange<QDimRMesh>;
    using IdxRangeQuadratureTheta = IdxRange<QDimThetaMesh>;
    using IdxRangeQuadratureRTheta = IdxRange<QDimRMesh, QDimThetaMesh>;
    using IdxQuadratureR = Idx<QDimRMesh>;
    using IdxQuadratureTheta = Idx<QDimThetaMesh>;
    using IdxQuadratureRTheta = Idx<QDimRMesh, QDimThetaMesh>;
    using IdxStepQuadratureR = IdxStep<QDimRMesh>;
    using IdxStepQuadratureTheta = IdxStep<QDimThetaMesh>;

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
    explicit PolarSplineFEMPoissonLikeAssembler(
            Field<double, IdxRangeQuadratureRTheta> int_volume,
            MetricTensorType const& metric_tensor,
            IdxRangeFull idx_range_full)
        : m_nbasis_theta(ddc::discrete_space<BSplinesTheta>().nbasis())
        , m_matrix_size(ddc::discrete_space<PolarBSplinesRTheta>().nbasis() - m_nbasis_theta)
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
                          m_idxrange_bsplines_r.remove_last(IdxStepBSR(1)),
                          m_idxrange_bsplines_theta))
    {
    }

    TensorBSplineStash get_tensor_bspline_stash() const
    {
        return m_tensor_assembler.get_tensor_bspline_stash();
    }

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

    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> get_tensor_nnz_per_row(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> nnz_per_row)
            const
    {
        int const n_tensor_basis = m_tensor_assembler.get_tensor_basis_idx_range().size();
        return Kokkos::subview(
                nnz_per_row,
                std::pair<int, int>(m_tensor_idx_offset, m_tensor_idx_offset + n_tensor_basis + 1));
    }

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

    void init_nnz_per_line(Kokkos::View<int*, Kokkos::LayoutRight> nnz_per_row) const
    {
        Kokkos::deep_copy(nnz_per_row, 0.);
        init_singular_nnz_per_line(nnz_per_row);
        m_tensor_assembler.init_nnz_per_line(get_tensor_nnz_per_row(nnz_per_row));
    }

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

    static KOKKOS_INLINE_FUNCTION IdxBSPolar to_polar(IdxBSRTheta idx)
    {
        return PolarBSplinesRTheta::template get_polar_index<PolarBSplinesRTheta>(idx);
    }
};
```


