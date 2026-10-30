

# File fem\_poisson\_assembler\_utils.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**fem\_poisson\_assembler\_utils.hpp**](fem__poisson__assembler__utils_8hpp.md)

[Go to the documentation of this file](fem__poisson__assembler__utils_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <tuple>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "indexed_tensor.hpp"
#include "metric_tensor_evaluator.hpp"
#include "tensor.hpp"
#include "vector_index_tools.hpp"


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

template <
        class QDimRMesh,
        class QDimThetaMesh,
        class IdxBatch,
        class CoeffAlpha,
        class CoeffBeta,
        concepts::MetricTensor MetricTensorType>
KOKKOS_FUNCTION std::tuple<
        double,
        DVector<typename QDimRMesh::continuous_dimension_type,
                typename QDimThetaMesh::continuous_dimension_type>>
get_test_weights(
        double basis_val_test,
        DVector<typename QDimRMesh::continuous_dimension_type::Dual,
                typename QDimThetaMesh::continuous_dimension_type::Dual> basis_derivs_test,
        Idx<QDimRMesh, QDimThetaMesh> idx_quad,
        IdxBatch idx_batch,
        CoeffAlpha const& coeff_alpha,
        CoeffBeta const& coeff_beta,
        MetricTensorType const& metric_tensor,
        DField<IdxRange<QDimRMesh, QDimThetaMesh>> int_volume)
{
    using R = typename QDimRMesh::continuous_dimension_type;
    using Theta = typename QDimThetaMesh::continuous_dimension_type;

    // Calculate coefficients at quadrature point
    Coord<R, Theta> coord(ddc::coordinate(idx_quad));
    const double alpha = coeff_alpha(idx_batch, coord);
    const double beta = coeff_beta(idx_batch, coord);

    using Spatial2DVectorSpace = VectorIndexSet<R, Theta>;

    Tensor inv_metric_tensor = metric_tensor.inverse(idx_quad);
    DTensor<Spatial2DVectorSpace, Spatial2DVectorSpace> inv_metric_tensor_2d;
    ddcHelper::assign_elements(inv_metric_tensor_2d, inv_metric_tensor);

    double test_mass_weight = int_volume(idx_quad) * beta * basis_val_test;
    DVector<R, Theta> test_stiffness_weight
            = int_volume(idx_quad) * alpha
              * tensor_mul(index<'i'>(basis_derivs_test), index<'i', 'j'>(inv_metric_tensor_2d));
    return {test_mass_weight, test_stiffness_weight};
}

} // namespace detail_poisson
```


