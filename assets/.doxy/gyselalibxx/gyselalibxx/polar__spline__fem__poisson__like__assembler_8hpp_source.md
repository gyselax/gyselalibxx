

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
#include "gauss_legendre_integration.hpp"
#include "math_tools.hpp"
#include "matrix_batch_csr.hpp"
#include "metric_tensor_evaluator.hpp"
#include "polar_spline_evaluator.hpp"
#include "quadrature_coeffs_nd.hpp"
#include "vector_field.hpp"
#include "view.hpp"
#include "volume_quadrature_nd.hpp"


namespace detail_poisson {

template <class... Dims>
Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> as_batched_vector_view(
        DField<IdxRange<Dims...>> field)
{
    static_assert(sizeof...(Dims) > 0);
    using LastDim = ddc::type_seq_element_t<sizeof...(Dims) - 1, ddc::detail::TypeSeq<Dims...>>;
    IdxRange<Dims...> idx_range = get_idx_range(field);
    std::size_t const matrix_size = idx_range.template extent<LastDim>();
    assert(field.size() == idx_range.size());
    return Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>(
            field.data_handle(),
            idx_range.size() / matrix_size,
            matrix_size);
}


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

template <typename DDom, std::enable_if_t<DDom::continuous_dimension_type::PERIODIC, bool> = true>
KOKKOS_FUNCTION Idx<DDom> mod_add(Idx<DDom> idx, IdxStep<DDom> idx_step, IdxRange<DDom> idx_range)
{
    if (idx - idx_range.front() < -idx_step) {
        return idx + idx_range.extents() + idx_step;
    } else if (idx_range.back() - idx < idx_step) {
        return idx + idx_step - idx_range.extents();
    } else {
        return idx + idx_step;
    }
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

template <typename QDimRMesh, typename QDimThetaMesh, typename BSplinesR, typename BSplinesTheta>
struct TensorBSplineStash
{
private:
    using R = typename BSplinesR::continuous_dimension_type;
    using Theta = typename BSplinesTheta::continuous_dimension_type;
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

public:
    DConstField<IdxRange<QDimRMesh, BSplinesR, ddc::Deriv<R>>> vals_r;
    DConstField<IdxRange<QDimThetaMesh, BSplinesTheta, ddc::Deriv<Theta>>> vals_theta;
    ConstField<Idx<BSplinesR>, IdxRange<QDimRMesh>> first_basis_r;
    ConstField<Idx<BSplinesTheta>, IdxRange<QDimThetaMesh>> first_basis_theta;

    KOKKOS_FUNCTION std::tuple<double, DVector<R_cov, Theta_cov>> vals_and_derivs(
            Idx<QDimRMesh, QDimThetaMesh> idx_quad,
            Idx<BSplinesR, BSplinesTheta> idx) const
    {
        Idx<QDimRMesh> idx_quad_r(idx_quad);
        Idx<QDimThetaMesh> idx_quad_theta(idx_quad);
        IdxStep<BSplinesR> offset_r = Idx<BSplinesR>(idx) - first_basis_r(idx_quad_r);
        IdxStep<BSplinesTheta> offset_theta
                = Idx<BSplinesTheta>(idx) - first_basis_theta(idx_quad_theta);
        DConstField<IdxRange<BSplinesR, ddc::Deriv<R>>> vals_r_at_quad = vals_r[idx_quad_r];
        DConstField<IdxRange<BSplinesTheta, ddc::Deriv<Theta>>> vals_theta_at_quad
                = vals_theta[idx_quad_theta];
        const double val_r = vals_r_at_quad(offset_r, IdxStep<ddc::Deriv<R>>(0));
        const double val_theta = vals_theta_at_quad(offset_theta, IdxStep<ddc::Deriv<Theta>>(0));
        DVector<R_cov, Theta_cov> derivs;
        ddcHelper::get<R_cov>(derivs)
                = vals_r_at_quad(offset_r, IdxStep<ddc::Deriv<R>>(1)) * val_theta;
        ddcHelper::get<Theta_cov>(derivs)
                = val_r * vals_theta_at_quad(offset_theta, IdxStep<ddc::Deriv<Theta>>(1));
        return {val_r * val_theta, derivs};
    }

    KOKKOS_FUNCTION double vals(
            Idx<QDimRMesh, QDimThetaMesh> idx_quad,
            Idx<BSplinesR, BSplinesTheta> idx) const
    {
        Idx<QDimRMesh> idx_quad_r(idx_quad);
        Idx<QDimThetaMesh> idx_quad_theta(idx_quad);
        IdxStep<BSplinesR> offset_r = Idx<BSplinesR>(idx) - first_basis_r(idx_quad_r);
        IdxStep<BSplinesTheta> offset_theta
                = Idx<BSplinesTheta>(idx) - first_basis_theta(idx_quad_theta);
        return vals_r[idx_quad_r](offset_r, IdxStep<ddc::Deriv<R>>(0))
               * vals_theta[idx_quad_theta](offset_theta, IdxStep<ddc::Deriv<Theta>>(0));
    }
};

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

template <typename QDimRMesh, typename QDimThetaMesh, typename BSplinesR, typename BSplinesTheta>
KOKKOS_FUNCTION IdxRange<QDimRMesh, QDimThetaMesh> get_support_quadrature(
        Idx<BSplinesR, BSplinesTheta> idx,
        Idx<QDimRMesh, QDimThetaMesh> idx_quad_front)
{
    using KnotsR = ddc::knot_discrete_dimension_t<BSplinesR>;
    using KnotsTheta = ddc::knot_discrete_dimension_t<BSplinesTheta>;

    const Idx<BSplinesR> idx_r(idx);
    const Idx<BSplinesTheta> idx_theta(idx);

    auto& bspl_r = ddc::discrete_space<BSplinesR>();
    auto& bspl_theta = ddc::discrete_space<BSplinesTheta>();

    // Find the cells on which the bspline is non-zero
    const Idx<KnotsR> start_non_zero_r(
            Kokkos::max(bspl_r.break_point_domain().front(), bspl_r.get_first_support_knot(idx_r)));
    const Idx<KnotsR> end_non_zero_r(
            Kokkos::min(bspl_r.break_point_domain().back(), bspl_r.get_last_support_knot(idx_r)));

    const Idx<KnotsTheta> start_non_zero_theta(bspl_theta.get_first_support_knot(idx_theta));
    const Idx<KnotsTheta> end_non_zero_theta(bspl_theta.get_last_support_knot(idx_theta));
    return get_quadrature_between_knots<QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta>(
            start_non_zero_r,
            end_non_zero_r,
            start_non_zero_theta,
            end_non_zero_theta,
            idx_quad_front);
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

    using IdxRangeBasisStashR = IdxRange<QDimRMesh, BSplinesR, ddc::Deriv<R>>;
    using IdxRangeBasisStashTheta = IdxRange<QDimThetaMesh, BSplinesTheta, ddc::Deriv<Theta>>;

public:
    using TensorBSplineStash = detail_poisson::
            TensorBSplineStash<QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta>;

private:
    static constexpr int s_n_gauss_legendre_r = BSplinesR::degree() + 1;
    static constexpr int s_n_gauss_legendre_theta = BSplinesTheta::degree() + 1;
    // The number of cells (in the radial direction) in which both types of basis splines can be found
    static constexpr int s_n_overlap_cells = PolarBSplinesRTheta::continuity + 1;

    // Number of cells over which a radial B-splines has its support
    // This is the case for b-splines which are not affected by the higher knot multiplicity at the boundary.
    static constexpr int s_n_non_zero_bases_r = BSplinesR::degree() + 1;

    // Number of cells over which a poloidal B-splines has its support
    static constexpr int s_n_non_zero_bases_theta = BSplinesTheta::degree() + 1;

    // Maximum number of quadrature points in the support of a tensor product B-spline
    static constexpr int s_n_quad_support_r = s_n_gauss_legendre_r * s_n_non_zero_bases_r;
    static constexpr int s_n_quad_support_theta
            = s_n_gauss_legendre_theta * s_n_non_zero_bases_theta;

    using QuadField = DField<IdxRangeQuadratureRTheta>;

    // A view on team scratch memory
    using ScratchView = Kokkos::View<
            double*,
            Kokkos::LayoutRight,
            Kokkos::DefaultExecutionSpace::scratch_memory_space,
            Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    const int m_nbasis_r;
    const int m_nbasis_theta;

    // Matrix size is equal to the number of Polar bspline
    const int m_matrix_size;

    // Domains
    IdxRangeBSPolar m_idxrange_fem_tensor_basis;
    IdxRangeBSR m_idxrange_bsplines_r;
    IdxRangeBSTheta m_idxrange_bsplines_theta;

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

public:
    explicit PolarSplineFEMPoissonLikeAssembler(
            Field<double, IdxRangeQuadratureRTheta> int_volume,
            MetricTensorType const& metric_tensor,
            IdxRangeFull idx_range_full)
        : m_nbasis_r(ddc::discrete_space<BSplinesR>().nbasis() - s_n_overlap_cells - 1)
        , m_nbasis_theta(ddc::discrete_space<BSplinesTheta>().nbasis())
        , m_matrix_size(ddc::discrete_space<PolarBSplinesRTheta>().nbasis() - m_nbasis_theta)
        , m_idxrange_fem_tensor_basis(
                  ddc::discrete_space<PolarBSplinesRTheta>().tensor_bspline_idx_range().remove_last(
                          IdxStepBSPolar {m_nbasis_theta}))
        , m_idxrange_bsplines_r(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                  IdxStepBSR {s_n_overlap_cells}))
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
    {
        stash_basis_eval(get_field(m_basis_vals_r), get_field(m_first_basis_r));
        stash_basis_eval(get_field(m_basis_vals_theta), get_field(m_first_basis_theta));
    }

    TensorBSplineStash get_tensor_bspline_stash() const
    {
        return TensorBSplineStash {
                get_const_field(m_basis_vals_r),
                get_const_field(m_basis_vals_theta),
                get_const_field(m_first_basis_r),
                get_const_field(m_first_basis_theta)};
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
        // Number of elements in the matrix that correspond to the splines
        // that cover the singular point
        constexpr int n_elements_singular
                = PolarBSplinesRTheta::n_singular_basis() * PolarBSplinesRTheta::n_singular_basis();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // polar bsplines which traverse the singular point and other polar bsplines
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

        const int n_matrix_elements = n_elements_singular + n_elements_overlap + n_elements_stencil;

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
        compute_tensor_tensor_col_idx(col_idx, nnz_per_row);
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
        compute_tensor_tensor_elements(
                coeff_alpha,
                coeff_beta,
                m_metric_tensor,
                values,
                col_idx,
                nnz_per_row);

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

        Kokkos::Profiling::pushRegion("(GSLX) PolarPoissonFillFemMatrix");
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

    void compute_tensor_tensor_col_idx(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    col_idx_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const
                    nnz_per_row_csr)
    {
        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();

        // Get index range for basis elements (last element removed due to homogeneous Dirichlet)
        IdxRangeBSTheta full_idx_range_theta
                = ddc::discrete_space<BSplinesTheta>().full_domain().take_first(
                        IdxStepBSTheta(ddc::discrete_space<BSplinesTheta>().nbasis()));

        IdxRangeBSR central_radial_bspline_idx_range(
                m_idxrange_bsplines_r.take_first(IdxStep<BSplinesR> {BSplinesR::degree()}));

        IdxRangeBSR idx_range_fem_r = ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                IdxStepBSR(PolarBSplinesRTheta::continuity + 1));

        IdxBSPolar idxrange_fem_non_singular_front = m_idxrange_fem_tensor_basis.front();

        const int n_singular = idxrange_singular.size();
        const std::source_location location = std::source_location::current();

        // Calculate the matrix elements following a stencil
        // Teams loop over rows
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        m_idxrange_fem_tensor_basis.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    // Calculate the index of the test b-spline
                    IdxBSPolar const idx_test_polar(
                            idxrange_fem_non_singular_front + IdxStepBSPolar {team.league_rank()});

                    // Calculate the radial and poloidal components of the test b-spline
                    const IdxBSRTheta idx_test(PolarBSplinesRTheta::get_2d_index(idx_test_polar));
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
                                        idx_range_fem_r.back() - idx_test_r));
                    IdxStepBSTheta idx_step_trial_theta_offset_min(-BSplinesTheta::degree());
                    IdxStepBSTheta idx_step_trial_theta_offset_max(BSplinesTheta::degree() + 1);

                    // Calculate the row index and the first possible column index
                    int const row_idx = idx_test_polar - idxrange_singular.front();
                    int col_offset
                            = n_singular * central_radial_bspline_idx_range.contains(idx_test_r);

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
                            int const col_idx = to_polar(idx_trial) - idxrange_singular.front();
                            const int aij_idx = nnz_per_row_csr(row_idx) + col_offset;
                            col_idx_csr(aij_idx) = col_idx;
                            col_offset++;
                        }
                        // Loop over the poloidal offset to a trial b-spline for any elements that
                        // are to the right of zeros following the elements already set
                        for (IdxBSTheta idx_trial_theta(first_periodic_idx_trial_theta);
                             idx_trial_theta < full_idx_range_theta.back() + 1;
                             ++idx_trial_theta) {
                            const IdxBSRTheta idx_trial(idx_trial_r, idx_trial_theta);
                            int const col_idx = to_polar(idx_trial) - idxrange_singular.front();
                            const int aij_idx = nnz_per_row_csr(row_idx) + col_offset;
                            col_idx_csr(aij_idx) = col_idx;
                            col_offset++;
                        }
                    }
                });

        Kokkos::Profiling::popRegion();
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
        IdxRangeQuadratureRTheta idx_range_quad_singular(
                m_idxrange_quadrature_r.take_first(
                        IdxStep<QDimRMesh> {s_n_overlap_cells * s_n_gauss_legendre_r}),
                m_idxrange_quadrature_theta);

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_int_volume;
        TensorBSplineStash bspline_stash = get_tensor_bspline_stash();

        IdxRangeBatch idx_range_batch = m_idx_range_batch;
        const int n_batch = idx_range_batch.size();
        const int n_singular = idxrange_singular.size();

        Kokkos::Profiling::pushRegion("(GSLX) PolarPoissonFillFemMatrix");
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
                                        int_volume_proxy,
                                        bspline_stash);

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

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_int_volume;
        TensorBSplineStash bspline_stash = get_tensor_bspline_stash();
        IdxRangeQuadratureRTheta full_quad_idx_range = m_idxrange_quadrature;

        IdxQuadratureRTheta idxrange_quadrature_front = m_idxrange_quadrature.front();

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
                                        int_volume_proxy,
                                        bspline_stash);

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
    void compute_tensor_tensor_elements(
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

        IdxRangeBSR central_radial_bspline_idx_range(
                m_idxrange_bsplines_r.take_first(IdxStep<BSplinesR> {BSplinesR::degree()}));

        IdxBSPolar idxrange_fem_non_singular_front = m_idxrange_fem_tensor_basis.front();

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_int_volume;
        TensorBSplineStash bspline_stash = get_tensor_bspline_stash();

        IdxRangeQuadratureRTheta full_quad_idx_range = m_idxrange_quadrature;

        IdxRangeBatch idx_range_batch = m_idx_range_batch;
        const int n_batch = idx_range_batch.size();
        const int n_singular = idxrange_singular.size();
        const int n_tensor_basis = m_idxrange_fem_tensor_basis.size();
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
                    // Calculate the index of the test b-spline
                    IdxBSPolar const idx_test_polar(
                            idxrange_fem_non_singular_front
                            + IdxStepBSPolar {team.league_rank() % n_tensor_basis});

                    // Calculate the radial and poloidal components of the test b-spline
                    const IdxBSRTheta idx_test(PolarBSplinesRTheta::get_2d_index(idx_test_polar));
                    const IdxBSR idx_test_r(idx_test);

                    // Calculate the row index and the first possible column index
                    int const row_idx = idx_test_polar - idxrange_singular.front();
                    int const csr_idx_first_stencil_element
                            = n_singular * central_radial_bspline_idx_range.contains(idx_test_r);

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
                    for (int csr_idx = nnz_per_row_csr[row_idx] + csr_idx_first_stencil_element;
                         csr_idx < nnz_per_row_csr[row_idx + 1];
                         ++csr_idx) {
                        int const col_idx = col_idx_csr[csr_idx];
                        IdxBSPolar const idx_trial_polar = idxrange_singular.front() + col_idx;
                        const IdxBSRTheta idx_trial
                                = PolarBSplinesRTheta::get_2d_index(idx_trial_polar);
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

        Kokkos::Profiling::popRegion();
    }

    template <class CoeffAlpha, class CoeffBeta, class IdxTestType>
    static KOKKOS_FUNCTION std::tuple<double, DVector<R, Theta>> get_test_weights(
            IdxTestType idx_test,
            IdxQuadratureRTheta idx_quad,
            IdxBatch idx_batch,
            CoeffAlpha const& coeff_alpha,
            CoeffBeta const& coeff_beta,
            MetricTensorType const& metric_tensor,
            DField<IdxRangeQuadratureRTheta> int_volume,
            TensorBSplineStash const& bspline_stash)
    {
        // Define the value and gradient of the test basis function
        auto [basis_val_test_space, basis_derivs_test_space] = [&]() {
            if constexpr (std::is_same_v<IdxTestType, Idx<PolarBSplinesRTheta>>) {
                return detail_poisson::get_polar_bspline_vals_and_derivs(idx_quad, idx_test);
            } else {
                static_assert(std::is_same_v<IdxTestType, Idx<BSplinesR, BSplinesTheta>>);
                return bspline_stash.vals_and_derivs(idx_quad, idx_test);
            }
        }();


        // Calculate coefficients at quadrature point
        Coord<R, Theta> coord(ddc::coordinate(idx_quad));
        const double alpha = coeff_alpha(idx_batch, coord);
        const double beta = coeff_beta(idx_batch, coord);

        using Spatial2DVectorSpace = VectorIndexSet<R, Theta>;

        Tensor inv_metric_tensor = metric_tensor.inverse(coord);
        DTensor<Spatial2DVectorSpace, Spatial2DVectorSpace> inv_metric_tensor_2d;
        ddcHelper::assign_elements(inv_metric_tensor_2d, inv_metric_tensor);

        double test_mass_weight = int_volume(idx_quad) * beta * basis_val_test_space;
        DVector<R, Theta> test_stiffness_weight = int_volume(idx_quad) * alpha
                                                  * tensor_mul(
                                                          index<'i'>(basis_derivs_test_space),
                                                          index<'i', 'j'>(inv_metric_tensor_2d));
        return {test_mass_weight, test_stiffness_weight};
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

    void init_nnz_per_line(Kokkos::View<int*, Kokkos::LayoutRight> nnz_per_row) const
    {
        Kokkos::Profiling::pushRegion("(GSLX) PolarPoissonInitNnz");
        size_t const mat_size = nnz_per_row.extent(0) - 1;
        IdxStepBSPolar radial_boundary_splines(m_nbasis_theta);
        IdxRangeBSPolar polar_bspl_idx_range
                = ddc::discrete_space<PolarBSplinesRTheta>().full_domain().remove_last(
                        radial_boundary_splines);

        Kokkos::deep_copy(nnz_per_row, 0.);

        assert(mat_size == polar_bspl_idx_range.size());
        // nnz per line
        Field<int, IdxRangeBSPolar>
                nnz(Kokkos::subview(nnz_per_row, std::pair<int, int>(1, mat_size + 1)),
                    polar_bspl_idx_range);

        size_t constexpr n_singular_basis = PolarBSplinesRTheta::n_singular_basis();
        size_t constexpr degree = BSplinesR::degree();
        size_t constexpr stencil_overlap = 2 * degree + 1;
        size_t const nbasis_theta_proxy = m_nbasis_theta;

        // The number of radial splines which overlap with a given radial spline to its
        // left or to its right
        IdxStepBSR n_one_side_overlap(BSplinesR::degree());

        int nnz_for_singular_rows = n_singular_basis + degree * nbasis_theta_proxy;

        // rows representing the bsplines which cover the singular domain
        // These overlap with other singular domain splines and degree radial splines
        Kokkos::parallel_for(
                "overlap singular radial",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(1, n_singular_basis + 1),
                KOKKOS_LAMBDA(const int k) { nnz_per_row(k) = k * nnz_for_singular_rows; });

        int nnz_sum = nnz_for_singular_rows * n_singular_basis;

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
                IdxStepBSPolar(outer_bsplines_2d.size()));
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


