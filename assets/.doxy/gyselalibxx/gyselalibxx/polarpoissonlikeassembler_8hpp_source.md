

# File polarpoissonlikeassembler.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**polarpoissonlikeassembler.hpp**](polarpoissonlikeassembler_8hpp.md)

[Go to the documentation of this file](polarpoissonlikeassembler_8hpp.md)


```C++
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


namespace detail_poisson {

template <typename BSplinesTheta>
KOKKOS_FUNCTION IdxStep<BSplinesTheta> theta_mod(IdxStep<BSplinesTheta> idx_theta)
{
    int n_theta = ddc::discrete_space<BSplinesTheta>().nbasis();
    while (idx_theta < 0)
        idx_theta += n_theta;
    while (idx_theta >= n_theta)
        idx_theta -= n_theta;
    return idx_theta;
}

template <typename BSplinesTheta, class IdxType>
KOKKOS_INLINE_FUNCTION IdxType theta_mod(IdxType idx)
{
    static_assert(ddc::is_discrete_element_v<IdxType>);
    static_assert(ddc::in_tags_v<BSplinesTheta, ddc::to_type_seq_t<IdxType>>);
    IdxRange<BSplinesTheta> idx_range_theta
            = ddc::discrete_space<BSplinesTheta>().full_domain().take_first(
                    IdxStep<BSplinesTheta>(ddc::discrete_space<BSplinesTheta>().nbasis()));
    while (ddc::select<BSplinesTheta>(idx) < idx_range_theta.front())
        idx += idx_range_theta.extents();
    while (ddc::select<BSplinesTheta>(idx) > idx_range_theta.back())
        idx -= idx_range_theta.extents();
    assert(idx_range_theta.contains(ddc::select<BSplinesTheta>(idx)));
    return idx;
}
template <typename PolarBSplinesRTheta, bool calculate_derivs = true>
KOKKOS_FUNCTION std::conditional_t<
        calculate_derivs,
        DVector<typename PolarBSplinesRTheta::R::Dual, typename PolarBSplinesRTheta::Theta::Dual>,
        void>
get_polar_bspline_vals_and_derivs(
        double& val,
        Coord<typename PolarBSplinesRTheta::R, typename PolarBSplinesRTheta::Theta> coord,
        Idx<PolarBSplinesRTheta> idx)
{
    using R = PolarBSplinesRTheta::R;
    using Theta = PolarBSplinesRTheta::Theta;
    std::array<double, PolarBSplinesRTheta::n_singular_basis()> singular_data;
    using BSplinesR = PolarBSplinesRTheta::BSplinesR_tag;
    using BSplinesTheta = PolarBSplinesRTheta::BSplinesTheta_tag;
    constexpr int n_non_zero_bases_r = BSplinesR::degree() + 1;
    constexpr int n_non_zero_bases_theta = BSplinesR::degree() + 1;
    std::array<double, n_non_zero_bases_r * n_non_zero_bases_theta> data;
    // Values of the polar basis splines around the singular point
    // at a given coordinate
    DSpan1D singular_vals(singular_data.data(), PolarBSplinesRTheta::n_singular_basis());
    // Values of the polar basis splines, that do not cover the singular point,
    // at a given coordinate
    DSpan2D vals(data.data(), n_non_zero_bases_r, n_non_zero_bases_theta);

    auto& polar_bspl = ddc::discrete_space<PolarBSplinesRTheta>();

    if (idx < Idx<PolarBSplinesRTheta>(PolarBSplinesRTheta::n_singular_basis())) {
        IdxStep<PolarBSplinesRTheta> offset
                = idx
                  - PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>().front();
        polar_bspl.eval_basis(singular_vals, vals, coord);
        val = singular_vals[offset.value()];
        if constexpr (calculate_derivs) {
            polar_bspl.eval_deriv(singular_vals, vals, coord, Idx<ddc::Deriv<R>>(1));
            double r_deriv = singular_vals[offset.value()];
            polar_bspl.eval_deriv(singular_vals, vals, coord, Idx<ddc::Deriv<Theta>>(1));
            double theta_deriv = singular_vals[offset.value()];
            DVector<typename R::Dual, typename Theta::Dual> derivs(r_deriv, theta_deriv);
            return derivs;
        } else {
            return;
        }
    } else {
        Idx<BSplinesR, BSplinesTheta> idx_front = polar_bspl.eval_basis(singular_vals, vals, coord);
        IdxStep<BSplinesR, BSplinesTheta> offset
                = PolarBSplinesRTheta::get_2d_index(idx) - idx_front;
        IdxStep<BSplinesR> ir(offset);
        IdxStep<BSplinesTheta> itheta(
                detail_poisson::theta_mod<BSplinesTheta>(IdxStep<BSplinesTheta>(offset)));

        val = vals(ir, itheta);
        if constexpr (calculate_derivs) {
            polar_bspl.eval_deriv(singular_vals, vals, coord, Idx<ddc::Deriv<R>>(1));
            double r_deriv = vals(ir, itheta);
            polar_bspl.eval_deriv(singular_vals, vals, coord, Idx<ddc::Deriv<Theta>>(1));
            double theta_deriv = vals(ir, itheta);
            DVector<typename R::Dual, typename Theta::Dual> derivs(r_deriv, theta_deriv);
            return derivs;
        } else {
            return;
        }
    }
}

template <typename QDimRMesh, typename QDimThetaMesh, typename BSplinesR, typename BSplinesTheta>
KOKKOS_FUNCTION IdxRange<QDimRMesh, QDimThetaMesh> get_quadrature_between_knots(
        Idx<ddc::knot_discrete_dimension_t<BSplinesR>> start_knot_r,
        Idx<ddc::knot_discrete_dimension_t<BSplinesR>> end_knot_r,
        Idx<ddc::knot_discrete_dimension_t<BSplinesTheta>> start_knot_theta,
        Idx<ddc::knot_discrete_dimension_t<BSplinesTheta>> end_knot_theta,
        Idx<QDimRMesh, QDimThetaMesh> idx_quad_front)
{
    using KnotsR = ddc::knot_discrete_dimension_t<BSplinesR>;
    using KnotsTheta = ddc::knot_discrete_dimension_t<BSplinesTheta>;

    const IdxRange<KnotsR> k_range_r(start_knot_r, end_knot_r - start_knot_r);
    const IdxRange<KnotsTheta> k_range_theta(start_knot_theta, end_knot_theta - start_knot_theta);

    constexpr int s_n_gauss_legendre_r = BSplinesR::degree() + 1;
    IdxStep<KnotsR> k_r_offset
            = k_range_r.front() - ddc::discrete_space<BSplinesR>().break_point_domain().front();
    Idx<QDimRMesh> q_r_offset
            = Idx<QDimRMesh>(idx_quad_front) + k_r_offset.value() * s_n_gauss_legendre_r;
    IdxStep<QDimRMesh> q_r_len(k_range_r.extents().value() * s_n_gauss_legendre_r);
    IdxRange<QDimRMesh> q_range_r(q_r_offset, q_r_len);

    constexpr int s_n_gauss_legendre_theta = BSplinesTheta::degree() + 1;
    IdxStep<KnotsTheta> k_theta_offset
            = k_range_theta.front()
              - ddc::discrete_space<BSplinesTheta>().break_point_domain().front();
    if (k_theta_offset < 0)
        k_theta_offset += ddc::discrete_space<BSplinesTheta>().nbasis();
    Idx<QDimThetaMesh> q_theta_offset = Idx<QDimThetaMesh>(idx_quad_front)
                                        + k_theta_offset.value() * s_n_gauss_legendre_theta;
    IdxStep<QDimThetaMesh> q_theta_len(k_range_theta.extents().value() * s_n_gauss_legendre_theta);
    IdxRange<QDimThetaMesh> q_range_theta(q_theta_offset, q_theta_len);
    assert(q_range_r.extents() > 0);
    assert(q_range_theta.extents() > 0);

    return IdxRange<QDimRMesh, QDimThetaMesh>(q_range_r, q_range_theta);
}

} // namespace detail_poisson


template <
        typename GridR,
        typename GridTheta,
        typename PolarBSplinesRTheta,
        typename SplineRThetaEvaluatorNullBound,
        typename QDimRMesh,
        typename QDimThetaMesh,
        class IdxRangeFull = IdxRange<GridR, GridTheta>>
class PolarSplineFEMPoissonLikeAssembler
{
    // TODO: Add a batch loop to operator()
    static_assert(
            std::is_same_v<IdxRangeFull, IdxRange<GridR, GridTheta>>,
            "PolarSplineFEMPoissonLikeAssembler is not yet batched");

public:
    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;

    static_assert(R::IS_CONTRAVARIANT);
    static_assert(Theta::IS_CONTRAVARIANT);

    struct InternalBatchDim
    {
    };

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

    using IdxRangeBatchedBSRTheta
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<IdxRangeFull>,
                    ddc::detail::TypeSeq<GridR, GridTheta>,
                    ddc::detail::TypeSeq<BSplinesR, BSplinesTheta>>>;

    using IdxRangeQuadratureR = IdxRange<QDimRMesh>;
    using IdxRangeQuadratureTheta = IdxRange<QDimThetaMesh>;
    using IdxRangeQuadratureRTheta = IdxRange<QDimRMesh, QDimThetaMesh>;
    using IdxQuadratureR = Idx<QDimRMesh>;
    using IdxQuadratureTheta = Idx<QDimThetaMesh>;
    using IdxQuadratureRTheta = Idx<QDimRMesh, QDimThetaMesh>;
    using IdxStepQuadratureR = IdxStep<QDimRMesh>;
    using IdxStepQuadratureTheta = IdxStep<QDimThetaMesh>;

    using ConstSpline2D = DConstField<IdxRangeBatchedBSRTheta>;

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
    IdxRangeQuadratureRTheta m_idxrange_quadrature_singular;
    IdxRangeQuadratureRTheta m_idxrange_quadrature;

    const int m_batch_idx {0}; // TODO: Remove when batching is supported

    Field<double, IdxRangeQuadratureRTheta> m_int_volume;

public:
    explicit PolarSplineFEMPoissonLikeAssembler(Field<double, IdxRangeQuadratureRTheta> int_volume)
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
        , m_int_volume(int_volume)
    {
    }
    template <typename Mapping>
    void operator()(
            std::unique_ptr<
                    MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>&
                    gko_matrix,
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
            std::optional<int> max_iter = std::nullopt,
            std::optional<double> res_tol = std::nullopt,
            std::optional<bool> batch_solver_logger = std::nullopt,
            std::optional<int> preconditioner_max_block_size = std::nullopt)
    {
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
                = m_nbasis_theta * std::min(int(1 + 2 * BSplinesTheta::degree()), m_nbasis_theta);
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

        gko_matrix = std::make_unique<
                MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>(
                1,
                m_matrix_size,
                n_matrix_elements,
                max_iter,
                res_tol,
                batch_solver_logger,
                preconditioner_max_block_size);
        auto [values, col_idx, nnz_per_row] = gko_matrix->get_batch_csr();
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
        gko_matrix->setup_solver();
    }

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

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_int_volume;

        Kokkos::Profiling::pushRegion("PolarPoissonFillFemMatrix");
        const std::source_location location = std::source_location::current();
        // Calculate the matrix elements corresponding to the B-splines which cover the singular point
        ddc::host_for_each(idxrange_singular, [&](IdxBSPolar const idx_test) {
            ddc::host_for_each(idxrange_singular, [&](IdxBSPolar const idx_trial) {
                // Calculate the weak integral
                double const element = ddc::parallel_transform_reduce(
                        location.function_name(),
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

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_int_volume;
        IdxRangeQuadratureRTheta
                full_quad_idx_range(m_idxrange_quadrature_r, m_idxrange_quadrature_theta);

        // Calculate the matrix elements where bspline products overlap the B-splines which cover the singular point
        ddc::host_for_each(idxrange_singular, [&](IdxBSPolar const idx_test) {
            ddc::host_for_each(idxrange_non_singular_near_centre, [&](IdxBSRTheta const idx_trial) {
                const IdxBSPolar idx_trial_polar(to_polar(idx_trial));
                const IdxBSR idx_trial_r(idx_trial);
                const IdxBSTheta idx_trial_theta(idx_trial);

                auto& bspl_r = ddc::discrete_space<BSplinesR>();
                auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

                // Find the index range covering the cells where both the test and trial functions are non-zero
                const Idx<KnotsR> start_non_zero_r(
                        std::
                                max(bspl_r.break_point_domain().front(),
                                    bspl_r.get_first_support_knot(idx_trial_r)));
                const Idx<KnotsR> end_non_zero_r(
                        std::
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
                                m_idxrange_quadrature.front());
                assert(quad_range.size() > 0);
                // Calculate the weak integral
                const std::source_location location = std::source_location::current();
                double element = ddc::parallel_transform_reduce(
                        location.function_name(),
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
        ddc::host_for_each(m_idxrange_fem_non_singular, [&](IdxBSPolar const idx_test_polar) {
            const IdxBSRTheta idx_test(PolarBSplinesRTheta::get_2d_index(idx_test_polar));
            const IdxBSR idx_test_r(idx_test);
            const IdxBSTheta idx_test_theta(idx_test);

            // Calculate the index of the elements that are already filled
            IdxRangeBSTheta remaining_theta(
                    idx_test_theta,
                    IdxStep<BSplinesTheta> {BSplinesTheta::degree() + 1});
            ddc::host_for_each(remaining_theta, [&](IdxBSTheta const idx_trial_theta) {
                IdxBSRTheta idx_trial(idx_test_r, idx_trial_theta);
                IdxBSPolar idx_trial_polar(
                        to_polar(detail_poisson::theta_mod<BSplinesTheta>(idx_trial)));
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
            IdxStepBSR n_remaining_r(std::
                                             min(IdxStepBSR(BSplinesR::degree()),
                                                 full_idx_range_r.back() - idx_test_r));
            IdxRangeBSR remaining_r(idx_test_r + 1, n_remaining_r);
            IdxRangeBSTheta relevant_theta(
                    idx_test_theta + ddc::discrete_space<BSplinesTheta>().nbasis()
                            - BSplinesTheta::degree(),
                    IdxStepBSTheta {2 * BSplinesTheta::degree() + 1});

            IdxRangeBSRTheta trial_idx_range(remaining_r, relevant_theta);

            ddc::host_for_each(trial_idx_range, [&](IdxBSRTheta const idx_trial) {
                IdxBSPolar idx_trial_polar(
                        to_polar(detail_poisson::theta_mod<BSplinesTheta>(idx_trial)));
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
        Coord<R, Theta> coord(ddc::coordinate(idx_quad));
        const double alpha = spline_evaluator(coord, coeff_alpha);
        const double beta = spline_evaluator(coord, coeff_beta);

        // Define the value and gradient of the test and trial basis functions
        double basis_val_test_space;
        DVector<R_cov, Theta_cov> basis_derivs_test_space
                = detail_poisson::get_polar_bspline_vals_and_derivs<
                        PolarBSplinesRTheta>(basis_val_test_space, coord, idx_test);
        double basis_val_trial_space;
        DVector<R_cov, Theta_cov> basis_derivs_trial_space
                = detail_poisson::get_polar_bspline_vals_and_derivs<
                        PolarBSplinesRTheta>(basis_val_trial_space, coord, idx_trial);

        MetricTensorEvaluator<Mapping, Coord<R, Theta>> get_metric_tensor(mapping);

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
        const IdxBSTheta idx_test_theta(
                detail_poisson::theta_mod<BSplinesTheta>(IdxBSTheta(idx_test)));
        const IdxBSTheta idx_trial_theta(
                detail_poisson::theta_mod<BSplinesTheta>(IdxBSTheta(idx_trial)));

        auto& bspl_r = ddc::discrete_space<BSplinesR>();
        auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

        const Idx<KnotsR> start_non_zero_r(
                std::
                        max(bspl_r.break_point_domain().front(),
                            std::
                                    max(bspl_r.get_first_support_knot(idx_test_r),
                                        bspl_r.get_first_support_knot(idx_trial_r))));
        const Idx<KnotsR> end_non_zero_r(
                std::
                        min(bspl_r.break_point_domain().back(),
                            std::
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
                std::max(first_support_knot_theta_test, first_support_knot_theta_trial));
        const Idx<KnotsTheta> end_non_zero_theta(
                std::min(last_support_knot_theta_test, last_support_knot_theta_trial));
        const IdxRangeQuadratureRTheta quad_range = detail_poisson::
                get_quadrature_between_knots<QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta>(
                        start_non_zero_r,
                        end_non_zero_r,
                        start_non_zero_theta,
                        end_non_zero_theta,
                        m_idxrange_quadrature.front());

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_int_volume;

        const IdxBSPolar idx_test_polar(to_polar(idx_test));
        const IdxBSPolar idx_trial_polar(to_polar(idx_trial));

        const std::source_location location = std::source_location::current();
        return ddc::parallel_transform_reduce(
                location.function_name(),
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
    static KOKKOS_INLINE_FUNCTION IdxBSPolar to_polar(IdxBSRTheta idx)
    {
        return PolarBSplinesRTheta::template get_polar_index<PolarBSplinesRTheta>(idx);
    }
};
```


