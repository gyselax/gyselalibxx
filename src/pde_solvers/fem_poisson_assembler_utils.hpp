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

/**
 * @brief Get a 2D Kokkos view of a 1D field with N batch dimensions.
 *
 * The returned view shares the memory of the field. The first dimension of the view
 * iterates over all the batch dimensions (in layout_right order), the second dimension
 * iterates over the last dimension of the field.
 *
 * @param[in] field A contiguous layout_right field whose last dimension is the dimension
 *          of the matrix equation and whose other dimensions are batch dimensions.
 *
 * @return A 2D Kokkos view of shape (number of batches, matrix size).
 */
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


/**
 * @brief Calculates the modulo idx_step in relation to the number of basis functions of
 * periodic B-splines.
 *
 * @tparam BSplines Periodic B-Splines.
 *
 * @param[in] idx_step A step between two B-spline indices.
 *
 * @return The corresponding step modulo the number of basis functions.
 */
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

/**
 * @brief Calculates the index which is inside the domain of periodic B-splines using the
 * periodicity properties.
 *
 * @tparam BSplines Periodic B-Splines.
 * @tparam IdxType The Type of the input index.
 *
 * @param[in] idx A multi-dimensional index including the B-spline index.
 *
 * @return The corresponding index inside the domain.
 */
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

/**
 * @brief Add an idx_step to an idx in such a way that the result remains in the idx_range.
 */
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

/**
 * @brief The values and first derivatives of the non-zero 1D B-splines at the 1D quadrature points.
 *
 * The tensor product B-splines and the quadrature grid both have a tensor product structure.
 * The value of a tensor product B-spline at a quadrature point can therefore be computed
 * from the values of the 1D B-splines at the 1D quadrature points. The memory required
 * to store these values is O(N_1 + N_2).
 *
 * The B-spline index stored in the stash is the offset from the first non-zero B-spline at the
 * quadrature point.
 *
 * @tparam QDim1Mesh Quadrature mesh in the first dimension.
 * @tparam QDim2Mesh Quadrature mesh in the second dimension.
 * @tparam BSplines1 Splines in the first dimension.
 * @tparam BSplines2 Splines in the second dimension.
 */
template <typename QDim1Mesh, typename QDim2Mesh, typename BSplines1, typename BSplines2>
struct TensorBSplineStash
{
private:
    using Dim1 = typename BSplines1::continuous_dimension_type;
    using Dim2 = typename BSplines2::continuous_dimension_type;
    using Dim1_cov = typename Dim1::Dual;
    using Dim2_cov = typename Dim2::Dual;

public:
    /// The values (derivative order 0) and derivatives (derivative order 1) of the B-splines in the first dimension.
    DConstField<IdxRange<QDim1Mesh, BSplines1, ddc::Deriv<Dim1>>> vals_1;
    /// The values (derivative order 0) and derivatives (derivative order 1) of the B-splines in the second dimension.
    DConstField<IdxRange<QDim2Mesh, BSplines2, ddc::Deriv<Dim2>>> vals_2;
    /// The first non-zero B-spline at each quadrature point in the first dimension.
    ConstField<Idx<BSplines1>, IdxRange<QDim1Mesh>> first_basis_1;
    /// The first non-zero B-spline at each quadrature point in the second dimension.
    ConstField<Idx<BSplines2>, IdxRange<QDim2Mesh>> first_basis_2;

    /**
     * @brief Get the value and the derivatives of a tensor product B-spline at a quadrature point.
     *
     * @param[in] idx_quad The index of the quadrature point.
     * @param[in] idx The 2D index of the tensor product B-spline.
     * @return The value and derivatives of the B-spline at the quadrature point.
     */
    KOKKOS_FUNCTION std::tuple<double, DVector<Dim1_cov, Dim2_cov>> vals_and_derivs(
            Idx<QDim1Mesh, QDim2Mesh> idx_quad,
            Idx<BSplines1, BSplines2> idx) const
    {
        Idx<QDim1Mesh> idx_quad_1(idx_quad);
        Idx<QDim2Mesh> idx_quad_2(idx_quad);
        IdxStep<BSplines1> offset_1 = Idx<BSplines1>(idx) - first_basis_1(idx_quad_1);
        IdxStep<BSplines2> offset_2 = Idx<BSplines2>(idx) - first_basis_2(idx_quad_2);
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

    /**
     * @brief Get the value of a tensor product B-spline at a quadrature point.
     *
     * @param[in] idx The 2D index of the tensor product B-spline.
     * @param[in] idx_quad The index of the quadrature point.
     * @return The value of the B-spline at the quadrature point.
     */
    KOKKOS_FUNCTION double vals(Idx<QDim1Mesh, QDim2Mesh> idx_quad, Idx<BSplines1, BSplines2> idx)
            const
    {
        Idx<QDim1Mesh> idx_quad_1(idx_quad);
        Idx<QDim2Mesh> idx_quad_2(idx_quad);
        IdxStep<BSplines1> offset_1 = Idx<BSplines1>(idx) - first_basis_1(idx_quad_1);
        IdxStep<BSplines2> offset_2 = Idx<BSplines2>(idx) - first_basis_2(idx_quad_2);
        return vals_1[idx_quad_1](offset_1, IdxStep<ddc::Deriv<Dim1>>(0))
               * vals_2[idx_quad_2](offset_2, IdxStep<ddc::Deriv<Dim2>>(0));
    }
};

/**
 * @brief Compute the 1D quadrature range between a provided pair of knots.
 *
 * Compute the range of quadrature points which are found between two knots. If the
 * B-splines are periodic then the front of the returned range is always inside the
 * domain. However, in order to return a contiguous range, the result may include indices
 * beyond the end of the domain. A modulo operator should be applied before using the indices.
 *
 * @tparam QDimMesh Quadrature mesh.
 * @tparam BSplines The 1D B-splines.
 *
 * @param[in] start_knot
 *      The index of the knot describing the lower bound of the domain of interest.
 * @param[in] end_knot
 *      The index of the knot describing the upper bound of the domain of interest.
 * @param[in] idx_quad_front
 *      The first index of the index range of the quadrature points.
 * @return
 *      The range of quadrature points in the specified domain.
 */
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

/**
 * @brief Compute the quadrature range between a provided set of knots.
 *
 * Compute the range of quadrature points which are found between a set of knots
 * in both dimensions. The index of the front of the returned range in the second
 * dimension is always inside the domain. However, in order to return a contiguous
 * range, the result may include indices beyond the end of the domain in the second
 * dimension. A modulo operator should be applied before using the indices.
 * See get_quadrature_between_knots_1d for implementation details.
 *
 * @tparam QDim1Mesh Quadrature mesh in the first dimension.
 * @tparam QDim2Mesh Quadrature mesh in the second dimension.
 * @tparam BSplines1 Splines in the first dimension.
 * @tparam BSplines2 Splines in the second dimension.
 *
 * @param[in] start_knot_1
 *      The index of the knot describing the lower bound of the domain of interest
 *      in the first dimension.
 * @param[in] end_knot_1
 *      The index of the knot describing the upper bound of the domain of interest
 *      in the first dimension.
 * @param[in] start_knot_2
 *      The index of the knot describing the lower bound of the domain of interest
 *      in the second dimension.
 * @param[in] end_knot_2
 *      The index of the knot describing the upper bound of the domain of interest
 *      in the second dimension.
 * @param[in] idx_quad_front
 *      The first index of the index range of the quadrature points.
 * @return
 *      The range of quadrature points in the specified domain.
 */
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

/**
 * @brief Compute the range of knots between which two 1D B-splines are both non-zero.
 *
 * If the B-splines are periodic then, in order to return a contiguous range, the result
 * may include knots beyond the end of the domain.
 *
 * @tparam BSplines The 1D B-splines.
 *
 * @param[in] idx_test
 *      The index of the first B-spline (the test function).
 * @param[in] idx_trial
 *      The index of the second B-spline (the trial function).
 * @return
 *      The range of knots describing the cells where both B-splines are non-zero.
 */
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

/**
 * @brief Get the 1D B-splines in a basis whose supports overlap the support of a given B-spline.
 *
 * The result is returned as 2 contiguous index ranges. The indices in the first range are
 * all smaller than the indices in the second range. This allows the B-splines to be
 * iterated over in order, even when the overlap wraps around the boundary of a periodic
 * basis. The second range is empty if the overlap does not wrap around the boundary.
 *
 * @tparam BSplines The 1D B-splines.
 *
 * @param[in] idx_test
 *      The index of the B-spline (the test function).
 * @param[in] idx_range_basis
 *      The index range of the B-splines in the basis. In a periodic dimension this must
 *      contain all the B-splines.
 * @return
 *      The 2 sorted index ranges containing all the B-splines in the basis which overlap
 *      the given B-spline.
 */
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

/**
 * @brief Compute the range of quadrature points in the support of a 1D B-spline.
 *
 * If the B-splines are periodic then, in order to return a contiguous range, the result
 * may include indices which are outside the domain. A modulo operator should be applied
 * before using the indices.
 *
 * @tparam QDimMesh Quadrature mesh.
 * @tparam BSplines The 1D B-splines.
 *
 * @param[in] idx
 *      The index of the B-spline.
 * @param[in] idx_quad_front
 *      The first index of the index range of the quadrature points.
 * @return
 *      The range of quadrature points in the support of the B-spline.
 */
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

/**
 * @brief Compute the range of quadrature points in the support of a tensor product B-spline.
 *
 * In order to return a contiguous range the result may include indices which are outside
 * the domain in periodic dimensions. A modulo operator should be applied before using the
 * indices.
 *
 * @tparam QDim1Mesh Quadrature mesh in the first dimension.
 * @tparam QDim2Mesh Quadrature mesh in the second dimension.
 * @tparam BSplines1 Splines in the first dimension.
 * @tparam BSplines2 Splines in the second dimension.
 *
 * @param[in] idx
 *      The 2D index of the tensor product B-spline.
 * @param[in] idx_quad_front
 *      The first index of the index range of the quadrature points.
 * @return
 *      The range of quadrature points in the support of the B-spline.
 */
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

/**
 * @brief Get the number of non-zero elements in the first i rows of a 1D stiffness matrix.
 *
 * Row j of the 1D stiffness matrix contains one non-zero element for each B-spline
 * whose support overlaps the support of the j-th B-spline (including the j-th B-spline
 * itself). In a periodic basis all rows contain the same number of non-zero elements.
 * In a non-periodic basis this number decreases near the boundaries. The result is
 * therefore the cumulative quantity required to build the nnz data structure of a CSR
 * matrix. When i == nbasis the result is the total number of non-zero elements in the
 * 1D stiffness matrix.
 *
 * @tparam BSplines The B-splines in the dimension of interest.
 *
 * @param[in] i
 *      The number of rows (counted from the start of the basis) whose non-zero elements are counted.
 * @param[in] nbasis
 *      The number of B-splines in the basis.
 * @return
 *      The number of non-zero elements in rows 0 to i-1 of the 1D stiffness matrix.
 */
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

/**
 * @brief Get the periodic image of a quadrature index which is found in a target range.
 *
 * The quadrature ranges returned by get_support_quadrature and get_quadrature_between_knots
 * are contiguous so they may extend past the end of the domain in periodic dimensions. This
 * function shifts an index by a period in each periodic dimension where it is outside the
 * target range. It can be used to move an index from such a range into the domain
 * (target_range = full_quad_idx_range) or from the domain into such a range.
 *
 * @tparam BSplines1 Splines in the first dimension.
 * @tparam BSplines2 Splines in the second dimension.
 * @tparam QDim1Mesh Quadrature mesh in the first dimension.
 * @tparam QDim2Mesh Quadrature mesh in the second dimension.
 *
 * @param[in] idx_quad
 *      The index of a quadrature point.
 * @param[in] target_range
 *      The range in which the returned index should be found.
 * @param[in] full_quad_idx_range
 *      The index range of all the quadrature points (one period).
 * @return
 *      The equivalent index inside target_range.
 */
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

/**
 * @brief Computes the test_stiffness_weight and test_mass_weight for a given quadrature point.
 *
 * The weak integral element is
 * @f$ w (\alpha \nabla B_i \cdot G^{-1} \nabla B_j + \beta B_i B_j) @f$
 * where @f$ w @f$ is the integral volume, @f$ B_i @f$ is the test function and @f$ B_j @f$ is
 * the trial function. This function computes @f$ w \alpha G^{-1} \nabla B_i @f$ and
 * @f$ w \beta B_i @f$ so that the weak integral element can be computed from the values and
 * derivatives of the trial function alone.
 *
 * @tparam QDim1Mesh Quadrature mesh in the first dimension.
 * @tparam QDim2Mesh Quadrature mesh in the second dimension.
 *
 * @param[in] basis_val_test
 *      The value of the test basis function at the quadrature point.
 * @param[in] basis_derivs_test
 *      The gradient of the test basis function at the quadrature point.
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
 * @return A tuple containing
 *          - test_mass_weight, the quantity @f$ w \beta B_i @f$ which
 *            multiplies the value of the trial function,
 *          - test_stiffness_weight, the quantity @f$ w \alpha G^{-1} \nabla B_i @f$ which
 *            multiplies the gradient of the trial function.
 */
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
