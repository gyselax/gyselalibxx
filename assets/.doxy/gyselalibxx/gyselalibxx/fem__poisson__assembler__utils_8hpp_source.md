

# File fem\_poisson\_assembler\_utils.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**fem\_poisson\_assembler\_utils.hpp**](fem__poisson__assembler__utils_8hpp.md)

[Go to the documentation of this file](fem__poisson__assembler__utils_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <array>
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
    using LastDim = ddc::type_seq_element_t<sizeof...(Dims) - 1, ddc::TypeSeq<Dims...>>;
    IdxRange<Dims...> idx_range = get_idx_range(field);
    std::size_t const matrix_size = idx_range.template extent<LastDim>();
    assert(field.size() == idx_range.size());
    return Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>(
            field.data_handle(),
            idx_range.size() / matrix_size,
            matrix_size);
}


template <typename BSplines>
KOKKOS_FUNCTION IdxStep<BSplines> periodic_mod(IdxStep<BSplines> idx_step)
{
    int n_basis = ddc::discrete_space<BSplines>().nbasis();
    while (idx_step < 0)
        idx_step += n_basis;
    while (idx_step >= n_basis)
        idx_step -= n_basis;
    return idx_step;
}

template <typename BSplines, class IdxType>
KOKKOS_INLINE_FUNCTION IdxType periodic_mod(IdxType idx)
{
    static_assert(ddc::is_discrete_element_v<IdxType>);
    static_assert(ddc::in_tags_v<BSplines, ddc::to_type_seq_t<IdxType>>);
    IdxRange<BSplines> idx_range_bsplines
            = ddc::discrete_space<BSplines>().full_domain().take_first(
                    IdxStep<BSplines>(ddc::discrete_space<BSplines>().nbasis()));
    while (ddc::select<BSplines>(idx) < idx_range_bsplines.front())
        idx += idx_range_bsplines.extents();
    while (ddc::select<BSplines>(idx) > idx_range_bsplines.back())
        idx -= idx_range_bsplines.extents();
    assert(idx_range_bsplines.contains(ddc::select<BSplines>(idx)));
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

template <typename QDim1Mesh, typename QDim2Mesh, typename BSplines1, typename BSplines2>
struct TensorBSplineStash
{
private:
    using Dim1 = typename BSplines1::continuous_dimension_type;
    using Dim2 = typename BSplines2::continuous_dimension_type;
    using Dim1_cov = typename Dim1::Dual;
    using Dim2_cov = typename Dim2::Dual;

public:
    DConstField<IdxRange<QDim1Mesh, BSplines1, ddc::Deriv<Dim1>>> vals_1;
    DConstField<IdxRange<QDim2Mesh, BSplines2, ddc::Deriv<Dim2>>> vals_2;
    ConstField<Idx<BSplines1>, IdxRange<QDim1Mesh>> first_basis_1;
    ConstField<Idx<BSplines2>, IdxRange<QDim2Mesh>> first_basis_2;

    KOKKOS_FUNCTION std::tuple<double, DVector<Dim1_cov, Dim2_cov>> vals_and_derivs(
            Idx<QDim1Mesh, QDim2Mesh> idx_quad,
            Idx<BSplines1, BSplines2> idx) const
    {
        Idx<QDim1Mesh> idx_quad_1(idx_quad);
        Idx<QDim2Mesh> idx_quad_2(idx_quad);
        IdxStep<BSplines1> offset_1 = Idx<BSplines1>(idx) - first_basis_1(idx_quad_1);
        IdxStep<BSplines2> offset_2 = Idx<BSplines2>(idx) - first_basis_2(idx_quad_2);
        if (BSplines1::is_periodic()) {
            offset_1 = periodic_mod(offset_1);
        }
        if (BSplines2::is_periodic()) {
            offset_2 = periodic_mod(offset_2);
        }
        DConstField<IdxRange<BSplines1, ddc::Deriv<Dim1>>> vals_1_at_quad = vals_1[idx_quad_1];
        DConstField<IdxRange<BSplines2, ddc::Deriv<Dim2>>> vals_2_at_quad = vals_2[idx_quad_2];
        const double val_1 = vals_1_at_quad(offset_1, IdxStep<ddc::Deriv<Dim1>>(0));
        const double val_2 = vals_2_at_quad(offset_2, IdxStep<ddc::Deriv<Dim2>>(0));
        DVector<Dim1_cov, Dim2_cov> derivs;
        ddcHelper::get<Dim1_cov>(derivs)
                = vals_1_at_quad(offset_1, IdxStep<ddc::Deriv<Dim1>>(1)) * val_2;
        ddcHelper::get<Dim2_cov>(derivs)
                = val_1 * vals_2_at_quad(offset_2, IdxStep<ddc::Deriv<Dim2>>(1));
        return {val_1 * val_2, derivs};
    }

    KOKKOS_FUNCTION double vals(Idx<QDim1Mesh, QDim2Mesh> idx_quad, Idx<BSplines1, BSplines2> idx)
            const
    {
        Idx<QDim1Mesh> idx_quad_1(idx_quad);
        Idx<QDim2Mesh> idx_quad_2(idx_quad);
        IdxStep<BSplines1> offset_1 = Idx<BSplines1>(idx) - first_basis_1(idx_quad_1);
        IdxStep<BSplines2> offset_2 = Idx<BSplines2>(idx) - first_basis_2(idx_quad_2);
        if (BSplines1::is_periodic()) {
            offset_1 = periodic_mod(offset_1);
        }
        if (BSplines2::is_periodic()) {
            offset_2 = periodic_mod(offset_2);
        }
        return vals_1[idx_quad_1](offset_1, IdxStep<ddc::Deriv<Dim1>>(0))
               * vals_2[idx_quad_2](offset_2, IdxStep<ddc::Deriv<Dim2>>(0));
    }
};

template <typename QDimMesh, typename BSplines>
KOKKOS_FUNCTION IdxRange<QDimMesh> get_quadrature_between_knots_1d(
        Idx<ddc::knot_discrete_dimension_t<BSplines>> start_knot,
        Idx<ddc::knot_discrete_dimension_t<BSplines>> end_knot,
        Idx<QDimMesh> idx_quad_front)
{
    using Knots = ddc::knot_discrete_dimension_t<BSplines>;

    constexpr int n_gauss_legendre = BSplines::degree() + 1;
    IdxStep<Knots> k_offset
            = start_knot - ddc::discrete_space<BSplines>().break_point_domain().front();
    if constexpr (BSplines::is_periodic()) {
        if (k_offset < 0)
            k_offset += ddc::discrete_space<BSplines>().ncells();
    }
    assert(k_offset >= 0);
    Idx<QDimMesh> q_offset = idx_quad_front + k_offset.value() * n_gauss_legendre;
    IdxStep<QDimMesh> q_len((end_knot - start_knot).value() * n_gauss_legendre);
    assert(q_len > 0);
    return IdxRange<QDimMesh>(q_offset, q_len);
}

template <typename QDim1Mesh, typename QDim2Mesh, typename BSplines1, typename BSplines2>
KOKKOS_FUNCTION IdxRange<QDim1Mesh, QDim2Mesh> get_quadrature_between_knots(
        Idx<ddc::knot_discrete_dimension_t<BSplines1>> start_knot_1,
        Idx<ddc::knot_discrete_dimension_t<BSplines1>> end_knot_1,
        Idx<ddc::knot_discrete_dimension_t<BSplines2>> start_knot_2,
        Idx<ddc::knot_discrete_dimension_t<BSplines2>> end_knot_2,
        Idx<QDim1Mesh, QDim2Mesh> idx_quad_front)
{
    return IdxRange<QDim1Mesh, QDim2Mesh>(
            get_quadrature_between_knots_1d<
                    QDim1Mesh,
                    BSplines1>(start_knot_1, end_knot_1, Idx<QDim1Mesh>(idx_quad_front)),
            get_quadrature_between_knots_1d<
                    QDim2Mesh,
                    BSplines2>(start_knot_2, end_knot_2, Idx<QDim2Mesh>(idx_quad_front)));
}

template <typename BSplines>
KOKKOS_INLINE_FUNCTION IdxRange<ddc::knot_discrete_dimension_t<BSplines>> get_overlapping_knots_1d(
        Idx<BSplines> idx_test,
        Idx<BSplines> idx_trial)
{
    using Knots = ddc::knot_discrete_dimension_t<BSplines>;

    auto& bspl = ddc::discrete_space<BSplines>();

    if constexpr (BSplines::is_periodic()) {
        idx_test = periodic_mod<BSplines>(idx_test);
        idx_trial = periodic_mod<BSplines>(idx_trial);

        IdxStep<Knots> span(BSplines::degree() + 1);

        Idx<Knots> first_support_knot_test = bspl.get_first_support_knot(idx_test);
        Idx<Knots> first_support_knot_trial = bspl.get_first_support_knot(idx_trial);
        Idx<Knots> last_support_knot_test = first_support_knot_test + span;
        Idx<Knots> last_support_knot_trial = first_support_knot_trial + span;

        // Shift one of the supports by a period if they only overlap across the boundary
        if (first_support_knot_test > last_support_knot_trial) {
            first_support_knot_trial += bspl.nbasis();
            last_support_knot_trial += bspl.nbasis();
        } else if (last_support_knot_test < first_support_knot_trial) {
            first_support_knot_test += bspl.nbasis();
            last_support_knot_test += bspl.nbasis();
        }
        const Idx<Knots> start_non_zero(
                Kokkos::max(first_support_knot_test, first_support_knot_trial));
        const Idx<Knots> end_non_zero(Kokkos::min(last_support_knot_test, last_support_knot_trial));
        return IdxRange<Knots>(start_non_zero, end_non_zero - start_non_zero);
    } else {
        const Idx<Knots> start_non_zero(
                Kokkos::
                        max(bspl.break_point_domain().front(),
                            Kokkos::
                                    max(bspl.get_first_support_knot(idx_test),
                                        bspl.get_first_support_knot(idx_trial))));
        const Idx<Knots> end_non_zero(
                Kokkos::
                        min(bspl.break_point_domain().back(),
                            Kokkos::
                                    min(bspl.get_last_support_knot(idx_test),
                                        bspl.get_last_support_knot(idx_trial))));
        return IdxRange<Knots>(start_non_zero, end_non_zero - start_non_zero);
    }
}

template <typename BSplines>
KOKKOS_FUNCTION std::array<IdxRange<BSplines>, 2> get_overlapping_bsplines_1d(
        Idx<BSplines> idx_test,
        IdxRange<BSplines> idx_range_basis)
{
    using IdxRangeBS = IdxRange<BSplines>;
    using IdxStepBS = IdxStep<BSplines>;
    assert(idx_range_basis.contains(idx_test));
    int const degree = BSplines::degree();
    int const nbasis = idx_range_basis.size();
    // The number of B-splines in the basis before and after the test B-spline
    int const n_before = (idx_test - idx_range_basis.front()).value();
    int const n_after = (idx_range_basis.back() - idx_test).value();
    IdxRangeBS const empty_range(idx_range_basis.front(), IdxStepBS(0));

    if constexpr (BSplines::is_periodic()) {
        assert(nbasis == int(ddc::discrete_space<BSplines>().nbasis()));
        if (nbasis < 2 * degree + 1) {
            // All B-splines overlap
            return {idx_range_basis, empty_range};
        } else if (n_before < degree) {
            // The overlap wraps around the lower boundary
            return {IdxRangeBS(idx_range_basis.front(), IdxStepBS(n_before + degree + 1)),
                    idx_range_basis.take_last(IdxStepBS(degree - n_before))};
        } else if (n_after < degree) {
            // The overlap wraps around the upper boundary
            return {idx_range_basis.take_first(IdxStepBS(degree - n_after)),
                    IdxRangeBS(idx_test - IdxStepBS(degree), IdxStepBS(n_after + degree + 1))};
        } else {
            return {IdxRangeBS(idx_test - IdxStepBS(degree), IdxStepBS(2 * degree + 1)),
                    empty_range};
        }
    } else {
        int const n_lower = Kokkos::min(n_before, degree);
        int const n_upper = Kokkos::min(n_after, degree);
        return {IdxRangeBS(idx_test - IdxStepBS(n_lower), IdxStepBS(n_lower + n_upper + 1)),
                empty_range};
    }
}

template <typename QDimMesh, typename BSplines>
KOKKOS_FUNCTION IdxRange<QDimMesh> get_support_quadrature_1d(
        Idx<BSplines> idx,
        Idx<QDimMesh> idx_quad_front)
{
    using Knots = ddc::knot_discrete_dimension_t<BSplines>;

    auto& bspl = ddc::discrete_space<BSplines>();

    // Find the cells on which the bspline is non-zero
    Idx<Knots> start_non_zero(bspl.get_first_support_knot(idx));
    Idx<Knots> end_non_zero(bspl.get_last_support_knot(idx));
    if constexpr (!BSplines::is_periodic()) {
        start_non_zero = Kokkos::max(bspl.break_point_domain().front(), start_non_zero);
        end_non_zero = Kokkos::min(bspl.break_point_domain().back(), end_non_zero);
    }
    return get_quadrature_between_knots_1d<
            QDimMesh,
            BSplines>(start_non_zero, end_non_zero, idx_quad_front);
}

template <typename QDim1Mesh, typename QDim2Mesh, typename BSplines1, typename BSplines2>
KOKKOS_FUNCTION IdxRange<QDim1Mesh, QDim2Mesh> get_support_quadrature(
        Idx<BSplines1, BSplines2> idx,
        Idx<QDim1Mesh, QDim2Mesh> idx_quad_front)
{
    return IdxRange<QDim1Mesh, QDim2Mesh>(
            get_support_quadrature_1d<
                    QDim1Mesh,
                    BSplines1>(Idx<BSplines1>(idx), Idx<QDim1Mesh>(idx_quad_front)),
            get_support_quadrature_1d<
                    QDim2Mesh,
                    BSplines2>(Idx<BSplines2>(idx), Idx<QDim2Mesh>(idx_quad_front)));
}

template <typename BSplines>
KOKKOS_INLINE_FUNCTION int n_nonzeros_in_rows(int i, int nbasis)
{
    int constexpr degree = BSplines::degree();
    assert(0 <= i && i <= nbasis);
    if constexpr (BSplines::is_periodic()) {
        return i * Kokkos::min(2 * degree + 1, nbasis);
    } else {
        assert(nbasis >= 2 * degree);
        if (i <= degree) {
            // \sum_{j=0}^{i-1} (degree + 1 + j)
            return i * (degree + 1) + i * (i - 1) / 2;
        }
        // The number of B-splines which are not summed
        int const n_after = nbasis - i;
        if (n_after <= degree) {
            int const total = nbasis * (2 * degree + 1) - degree * (degree + 1);
            return total - (n_after * (degree + 1) + n_after * (n_after - 1) / 2);
        }
        return degree * (degree + 1) + degree * (degree - 1) / 2 + (i - degree) * (2 * degree + 1);
    }
}

template <typename BSplines1, typename BSplines2, typename QDim1Mesh, typename QDim2Mesh>
KOKKOS_INLINE_FUNCTION Idx<QDim1Mesh, QDim2Mesh> periodic_wrap(
        Idx<QDim1Mesh, QDim2Mesh> idx_quad,
        IdxRange<QDim1Mesh, QDim2Mesh> const& target_range,
        IdxRange<QDim1Mesh, QDim2Mesh> const& full_quad_idx_range)
{
    if constexpr (BSplines1::is_periodic()) {
        IdxRange<QDim1Mesh> target_range_1(target_range);
        if (Idx<QDim1Mesh>(idx_quad) < target_range_1.front()) {
            idx_quad += full_quad_idx_range.template extent<QDim1Mesh>();
        } else if (Idx<QDim1Mesh>(idx_quad) > target_range_1.back()) {
            idx_quad -= full_quad_idx_range.template extent<QDim1Mesh>();
        }
    }
    if constexpr (BSplines2::is_periodic()) {
        IdxRange<QDim2Mesh> target_range_2(target_range);
        if (Idx<QDim2Mesh>(idx_quad) < target_range_2.front()) {
            idx_quad += full_quad_idx_range.template extent<QDim2Mesh>();
        } else if (Idx<QDim2Mesh>(idx_quad) > target_range_2.back()) {
            idx_quad -= full_quad_idx_range.template extent<QDim2Mesh>();
        }
    }
    assert(target_range.contains(idx_quad));
    return idx_quad;
}

template <
        class QDim1Mesh,
        class QDim2Mesh,
        class IdxBatch,
        class CoeffAlpha,
        class CoeffBeta,
        concepts::MetricTensor MetricTensorType>
KOKKOS_FUNCTION std::tuple<
        double,
        DVector<typename QDim1Mesh::continuous_dimension_type,
                typename QDim2Mesh::continuous_dimension_type>>
get_test_weights(
        double basis_val_test,
        DVector<typename QDim1Mesh::continuous_dimension_type::Dual,
                typename QDim2Mesh::continuous_dimension_type::Dual> basis_derivs_test,
        Idx<QDim1Mesh, QDim2Mesh> idx_quad,
        IdxBatch idx_batch,
        CoeffAlpha const& coeff_alpha,
        CoeffBeta const& coeff_beta,
        MetricTensorType const& metric_tensor,
        DField<IdxRange<QDim1Mesh, QDim2Mesh>> int_volume)
{
    using Dim1 = typename QDim1Mesh::continuous_dimension_type;
    using Dim2 = typename QDim2Mesh::continuous_dimension_type;

    // Calculate coefficients at quadrature point
    Coord<Dim1, Dim2> coord(ddc::coordinate(idx_quad));
    const double alpha = coeff_alpha(idx_batch, coord);
    const double beta = coeff_beta(idx_batch, coord);

    using Spatial2DVectorSpace = VectorIndexSet<Dim1, Dim2>;

    Tensor inv_metric_tensor = metric_tensor.inverse(idx_quad);
    DTensor<Spatial2DVectorSpace, Spatial2DVectorSpace> inv_metric_tensor_2d;
    ddcHelper::assign_elements(inv_metric_tensor_2d, inv_metric_tensor);

    double test_mass_weight = int_volume(idx_quad) * beta * basis_val_test;
    DVector<Dim1, Dim2> test_stiffness_weight
            = int_volume(idx_quad) * alpha
              * tensor_mul(index<'i'>(basis_derivs_test), index<'i', 'j'>(inv_metric_tensor_2d));
    return {test_mass_weight, test_stiffness_weight};
}

} // namespace detail_poisson
```


