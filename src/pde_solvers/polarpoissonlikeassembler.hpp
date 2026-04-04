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

/// The tag for the batch dimension for the equation.
struct BatchDDim
{
};

template <class ElementType, class... GridDims>
auto to_batch_access(Field<ElementType, IdxRange<GridDims...>> field)
{
    constexpr std::size_t ndims = sizeof...(GridDims);
    using BatchTypeSeq = type_seq_range_t<ddc::detail::TypeSeq<GridDims...>, 0, ndims - 2>;
    using PolarTypeSeq = ddc::type_seq_remove_t<ddc::detail::TypeSeq<GridDims...>, BatchTypeSeq>;
    using OutTypeSeq = ddc::type_seq_cat_t<ddc::detail::TypeSeq<BatchDDim>, PolarTypeSeq>;
    using IdxRangeBatch = ddc::detail::convert_type_seq_to_discrete_domain_t<BatchTypeSeq>;
    using IdxRangePolar = ddc::detail::convert_type_seq_to_discrete_domain_t<PolarTypeSeq>;
    using IdxRangeOut = ddc::detail::convert_type_seq_to_discrete_domain_t<OutTypeSeq>;

    IdxRange<GridDims...> full_idx_range(get_idx_range(field));
    IdxRangeBatch batch_idx_range(full_idx_range);
    IdxRangePolar polar_idx_range(full_idx_range);

    Idx<BatchDDim> batch_start(0);
    IdxStep<BatchDDim> batch_len(batch_idx_range.size());
    IdxRange<BatchDDim> solver_batch_idx_range(batch_start, batch_len);

    IdxRangeOut new_idx_range(solver_batch_idx_range, polar_idx_range);

    assert(new_idx_range.size() == full_idx_range.size());
    return Field<ElementType, IdxRangeOut>(field.data_handle(), new_idx_range);
}


/**
    * @brief Calculates the modulo idx_theta in relation to cells number along  @f$ \theta @f$ direction.
    *
    * @tparam BSplinesTheta Periodic B-Splines in poloidal direction.
    * 
    * @param[in] idx_theta @f$ \theta @f$ index.
    *
    * @return The corresponding indice modulo @f$ \theta @f$ direction cells number
    */
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

/**
    * @brief Calculates the index which is inside the poloidal domain using the periodicity properties.
    * 
    * @tparam BSplinesTheta Periodic B-Splines in poloidal direction.
    * @tparam IdxType The Type of the input index.
    *
    * @param[in] idx A multi-dimensional index including the polar bspline index.
    *
    * @return The corresponding index inside the domain.
    */
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
/**
 * @brief Get the value and derivative of the specified polar bspline at the specified quadrature point.
 *
 * This method calculates the value and the derivatives of polar bsplines. It is templated by
 * calculate_derivs to avoid code duplication between get_polar_bspline_vals_and_derivs and
 * get_polar_bspline_vals. The calling method should not need to use the template parameter.
 * 
 * @tparam PolarBSplinesRTheta The Polar B-Splines to be used.
 * @tparam calculate_derivs Returns the derivatives at the coordinate if true. If false, 
 *     only the value is calculated.
 *
 * @param[out] val
 *      The value of the specified polar bspline at the specified point.
 * @param[in] coord
 *      The coordinate where the value of the polar bspline should be calculated.
 * @param[in] idx
 *      The polar bspline of interest.
 * @return The derivative of the polar bspline (only returned if calculate_derivs is true).
 */
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

        assert(0 <= ir);
        assert(ir < n_non_zero_bases_r);
        assert(0 <= itheta);
        assert(itheta < n_non_zero_bases_theta);

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

/**
 * @brief Compute the quadrature range between a provided set of knots.
 *
 * Compute the range of quadrature points which are found between a set of knots
 * in both the radial and poloidal directions. In order to return a contiguous range
 * the result may include indices which are outside the domain. A modulo operator
 * should be applied before using the indices.
 * 
 * @tparam QDimRMesh Quadrature mesh in radial direction.
 * @tparam QDimThetaMesh Quadrature mesh in poloidal direction.
 * @tparam BSplinesR Splines in radial direction.
 * @tparam BSplinesTheta Splines in poloidal direction.
 *
 * @param[in] start_knot_r
 *      The index of the knot describing the lower bound of the domain of interest
 *      in the radial direction.
 * @param[in] end_knot_r
 *      The index of the knot describing the upper bound of the domain of interest
 *      in the radial direction.
 * @param[in] start_knot_theta
 *      The index of the knot describing the lower bound of the domain of interest
 *      in the poloidal direction.
 * @param[in] end_knot_theta
 *      The index of the knot describing the upper bound of the domain of interest
 *      in the poloidal direction.
 * @param[in] idx_quad_front
 *      The first index of the index range of the quadrature points.
 * @return 
 *      The range of quadrature points in the specified domain.
 */
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
 * @tparam SplineRThetaEvaluatorNullBound The type of the 2D (cross-product) spline evaluator.
 * @tparam QDimRMesh The radial quadrature grid type.
 * @tparam QDimThetaMesh The poloidal quadrature grid type.
 * @tparam IdxRangeFull The full index range of @f$ \phi @f$ including any batch dimensions.
 */

template <
        typename GridR,
        typename GridTheta,
        typename PolarBSplinesRTheta,
        typename SplineRThetaEvaluatorNullBound,
        typename QDimRMesh,
        typename QDimThetaMesh>
class PolarSplineFEMPoissonLikeAssembler
{
public:
    /// The radial dimension
    using R = typename GridR::continuous_dimension_type;
    /// The poloidal dimension
    using Theta = typename GridTheta::continuous_dimension_type;

    static_assert(R::IS_CONTRAVARIANT);
    static_assert(Theta::IS_CONTRAVARIANT);

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

    using IdxRangeFull = IdxRange<detail_poisson::BatchDDim, GridR, GridTheta>;

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

    using IdxRangeBatchedBSRTheta = IdxRange<detail_poisson::BatchDDim, BSplinesR, BSplinesTheta>;

    using IdxBatch = Idx<detail_poisson::BatchDDim>;

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

    using ConstSplineBatched2D = DConstField<IdxRangeBatchedBSRTheta>;

    using ConstSpline2D = DConstField<IdxRangeBSRTheta>;

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
    IdxRangeBSPolar m_idxrange_fem_tensor_basis;
    IdxRangeBSR m_idxrange_bsplines_r;
    IdxRangeBSTheta m_idxrange_bsplines_theta;

    IdxRangeQuadratureR m_idxrange_quadrature_r;
    IdxRangeQuadratureTheta m_idxrange_quadrature_theta;
    IdxRangeQuadratureRTheta m_idxrange_quadrature;

    DField<IdxRangeQuadratureRTheta> m_int_volume;

public:
    /**
     * @brief Instantiate the assembler operator.
     *
     * @param int_volume The initialised field of Jacobian values of the mapping. 
     */
    explicit PolarSplineFEMPoissonLikeAssembler(Field<double, IdxRangeQuadratureRTheta> int_volume)
        : m_nbasis_r(ddc::discrete_space<BSplinesR>().nbasis() - m_n_overlap_cells - 1)
        , m_nbasis_theta(ddc::discrete_space<BSplinesTheta>().nbasis())
        , m_matrix_size(ddc::discrete_space<PolarBSplinesRTheta>().nbasis() - m_nbasis_theta)
        , m_idxrange_fem_tensor_basis(
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
        , m_idxrange_quadrature(m_idxrange_quadrature_r, m_idxrange_quadrature_theta)
        , m_int_volume(int_volume)
    {
    }

    /**
     * @brief Compute the sparsity pattern for the stiffness matrix.
     * 
     * @param[out] gko_matrix The pointer to the assembled matrix.
     * @param[in] n_batch The number of matrix equations to solve.
     * @param[in] max_iter
     *      The maximum number of iterations possible for the batched CSR solver.
     * @param[in] res_tol
     *      The residual tolerance for the batched CSR solver. Be careful! the relative residual
     *      provided here, will be used as "implicit residual" in ginkgo solver.
     * @param[in] batch_solver_logger
     *      Indicates whether log information such as the residual and the number of iterations
     *      should be monitored.
     * @param[in] preconditioner_max_block_size
     *      The maximum size of the Jacobi preconditioner used by the batched CSR solver.
     *
     * @tparam Mapping A class describing a mapping from curvilinear coordinates to Cartesian coordinates.
     */
    void setup_sparse_matrix(
            std::unique_ptr<
                    MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>&
                    gko_matrix,
            std::size_t n_batch,
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
                n_batch,
                m_matrix_size,
                n_matrix_elements,
                max_iter,
                res_tol,
                batch_solver_logger,
                preconditioner_max_block_size);
        auto [values, col_idx, nnz_per_row] = gko_matrix->get_batch_csr();
        init_nnz_per_line(nnz_per_row);
        compute_singular_singular_col_idx(col_idx, nnz_per_row);
        compute_singular_tensor_col_idx(col_idx, nnz_per_row);
        compute_tensor_tensor_col_idx(col_idx, nnz_per_row);
    }

    /**
     * @brief Assemble the stiffness matrix.
     * 
     * @param[out] gko_matrix The pointer to the assembled matrix.
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @param[in] spline_evaluator
     *      An evaluator for evaluating 2D splines on @f$(r,\theta)@f$.
     *
     * @tparam Mapping A class describing a mapping from curvilinear coordinates to Cartesian coordinates.
     */
    template <typename Mapping>
    void operator()(
            std::unique_ptr<
                    MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>> const&
                    gko_matrix,
            ConstSplineBatched2D coeff_alpha,
            ConstSplineBatched2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator)
    {
        //CSR data storage
        auto [values, col_idx, nnz_per_row] = gko_matrix->get_batch_csr();

        assert(get_idx_range<detail_poisson::BatchDDim>(coeff_alpha).size() == values.extent(0));

        compute_singular_singular_elements(
                coeff_alpha,
                coeff_beta,
                mapping,
                spline_evaluator,
                values,
                col_idx,
                nnz_per_row);
        compute_singular_tensor_elements(
                coeff_alpha,
                coeff_beta,
                mapping,
                spline_evaluator,
                values,
                col_idx,
                nnz_per_row);
        compute_tensor_tensor_elements(
                coeff_alpha,
                coeff_beta,
                mapping,
                spline_evaluator,
                values,
                col_idx,
                nnz_per_row);

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

        Kokkos::Profiling::pushRegion("PolarPoissonFillFemMatrix");
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
     * @brief Computes the column indices of the matrix elements corresponding to the
     * inner products of tensor basis functions and tensor basis functions.
     *
     * @param[out] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @param[in] nnz_per_row_csr
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
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

    /**
     * @brief Computes the matrix elements corresponding to the inner products of singular
     * basis functions and singular basis functions.
     *
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @param[in] spline_evaluator
     *      An evaluator for evaluating 2D splines on @f$(r,\theta)@f$.
     * @param[out] values_csr
     *             A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[in] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix).
     * @param[in] nnz_per_row_csr
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    template <class Mapping>
    void compute_singular_singular_elements(
            ConstSplineBatched2D coeff_alpha,
            ConstSplineBatched2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
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
                        IdxStep<QDimRMesh> {m_n_overlap_cells * s_n_gauss_legendre_r}),
                m_idxrange_quadrature_theta);

        DField<IdxRangeQuadratureRTheta> int_volume_proxy = m_int_volume;

        const int n_batch = values_csr.extent(0);
        const int n_singular = idxrange_singular.size();

        Kokkos::Profiling::pushRegion("PolarPoissonFillFemMatrix");
        const std::source_location location = std::source_location::current();
        // Calculate the matrix elements corresponding to the B-splines which cover the singular point
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        n_batch * n_singular * n_singular,
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const int idx_batch = team.league_rank() / (n_singular * n_singular);
                    const int row_col_idx = team.league_rank() % (n_singular * n_singular);
                    const int row_idx = row_col_idx / n_singular;
                    const int col_idx = team.league_rank() % n_singular;
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
                                sum += weak_integral_element(
                                        idx_test,
                                        idx_trial,
                                        idx_quad,
                                        coeff_alpha[IdxBatch {idx_batch}],
                                        coeff_beta[IdxBatch {idx_batch}],
                                        spline_evaluator,
                                        mapping,
                                        int_volume_proxy);
                            },
                            element);
                    const int csr_idx_singular_area = nnz_per_row_csr(row_idx) + col_idx;
                    //Fill the dense matrix corresponding to the b-splines on the singular point
                    values_csr(idx_batch, csr_idx_singular_area) = element;
                });
    }

    /**
     * @brief Computes the matrix elements corresponding to the inner products of singular
     * basis functions and tensor basis functions.
     *
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @param[in] spline_evaluator
     *      An evaluator for evaluating 2D splines on @f$(r,\theta)@f$.
     * @param[out] values_csr
     *             A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[in] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @param[in] nnz_per_row_csr
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    template <class Mapping>
    void compute_singular_tensor_elements(
            ConstSplineBatched2D coeff_alpha,
            ConstSplineBatched2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
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
        IdxRangeQuadratureRTheta full_quad_idx_range = m_idxrange_quadrature;

        IdxQuadratureRTheta idxrange_quadrature_front = m_idxrange_quadrature.front();

        const int n_batch = values_csr.extent(0);
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
                    const int idx_batch(team.league_rank() / (n_singular * n_overlapping_singular));
                    const int idx_test_trial(
                            team.league_rank() % (n_singular * n_overlapping_singular));
                    const IdxStepBSPolar idx_step_test(idx_test_trial / n_overlapping_singular);
                    const IdxStepBSPolar idx_step_trial(
                            team.league_rank() % n_overlapping_singular);
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

                                sum += weak_integral_element<Mapping>(
                                        idx_test,
                                        idx_trial_polar,
                                        idx_quad,
                                        coeff_alpha[IdxBatch {idx_batch}],
                                        coeff_beta[IdxBatch {idx_batch}],
                                        spline_evaluator,
                                        mapping,
                                        int_volume_proxy);
                            },
                            element);

                    const int row_idx = idx_step_test.value();
                    const int col_idx = idx_step_trial.value() + n_singular;

                    //a_ij
                    values_csr(idx_batch, nnz_per_row_csr(row_idx) + col_idx) = element;
                    //a_ji
                    values_csr(idx_batch, nnz_per_row_csr(col_idx) + row_idx) = element;
                });
    }

    /**
     * @brief Computes the matrix elements corresponding to the inner products of tensor
     * basis functions and tensor basis functions.
     *
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @param[in] spline_evaluator
     *      An evaluator for evaluating 2D splines on @f$(r,\theta)@f$.
     * @param[out] values_csr
     *             A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[out] col_idx_csr
     *             A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @param[in] nnz_per_row_csr
     *               A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     */
    template <class Mapping>
    void compute_tensor_tensor_elements(
            ConstSplineBatched2D coeff_alpha,
            ConstSplineBatched2D coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
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

        IdxRangeQuadratureRTheta full_quad_idx_range = m_idxrange_quadrature;

        const int n_batch = values_csr.extent(0);
        const int n_singular = idxrange_singular.size();
        const std::source_location location = std::source_location::current();

        // Calculate the matrix elements following a stencil
        // Teams loop over rows
        Kokkos::parallel_for(
                location.function_name(),
                Kokkos::TeamPolicy<>(
                        Kokkos::DefaultExecutionSpace(),
                        n_batch * m_idxrange_fem_tensor_basis.size(),
                        Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const int idx_batch = team.league_rank() / m_idxrange_fem_tensor_basis.size();
                    // Calculate the index of the test b-spline
                    IdxBSPolar const idx_test_polar(
                            idxrange_fem_non_singular_front
                            + IdxStepBSPolar {
                                    team.league_rank() % m_idxrange_fem_tensor_basis.size()});

                    // Calculate the radial and poloidal components of the test b-spline
                    const IdxBSRTheta idx_test(PolarBSplinesRTheta::get_2d_index(idx_test_polar));
                    const IdxBSR idx_test_r(idx_test);

                    // Calculate the row index and the first possible column index
                    int const row_idx = idx_test_polar - idxrange_singular.front();
                    int const csr_idx_first_stencil_element
                            = n_singular * central_radial_bspline_idx_range.contains(idx_test_r);

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
                                coeff_alpha[IdxBatch {idx_batch}],
                                coeff_beta[IdxBatch {idx_batch}],
                                spline_evaluator,
                                mapping,
                                full_quad_idx_range,
                                int_volume_proxy);
                        values_csr(idx_batch, csr_idx) = element;
                    }
                });

        Kokkos::Profiling::popRegion();
    }

    /**
     * @brief Computes a quadrature summand for computing the integral
     *      corresponding to the inner product.
     *
     * Inner product of the test and trial spline is computed using a 
     * quadrature. This function returns one summand of the quadrature for 
     * the quadrature point given by the indices.
     *
     * @param[in] idx_test
     *      The index of the test basis spline.
     * @param[in] idx_trial
     *      The index of the trial basis spline.
     * @param[in] idx_quad
     *      The index for the point in the quadrature scheme.
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] spline_evaluator
     *      An evaluator for evaluating 2D splines on @f$(r,\theta)@f$.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @param[in] int_volume
     *      The integral volume associated with each point used in the quadrature scheme.
     * @return Value of the quadrature summand
     */
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

    /**
     * @brief Computes the matrix element corresponding to two tensor product splines
     *        with index idx_test and idx_trial
     *
     * @param[in] team
     *      The team of threads from which this function is called.
     * @param[in] idx_test
     *      The index for polar B-spline in the test space.
     * @param[in] idx_trial
     *      The index for polar B-spline in the trial space.
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] evaluator
     *      An evaluator for evaluating 2D splines on @f$ (r, \theta) @f$.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     * @param[in] full_quad_idx_range
     *      The index range of all the quadrature points.
     * @param[in] int_volume
     *      The field describing the quadrature coefficients including the
     *      Jacobian multiplication factor responsible for ensuring the correct
     *      volume for the integral.
     * @return 
     *      The value of the matrix element.
     */
    template <class Mapping>
    static KOKKOS_FUNCTION double get_matrix_stencil_element(
            const Kokkos::TeamPolicy<>::member_type& team,
            IdxBSRTheta idx_test,
            IdxBSRTheta idx_trial,
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping,
            IdxRangeQuadratureRTheta const& full_quad_idx_range,
            DField<IdxRangeQuadratureRTheta> int_volume)
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

        const IdxBSPolar idx_test_polar(to_polar(idx_test));
        const IdxBSPolar idx_trial_polar(to_polar(idx_trial));

        double result = 0;
        Kokkos::parallel_reduce(
                Kokkos::TeamThreadMDRange(
                        team,
                        quad_range.template extent<QDimRMesh>(),
                        quad_range.template extent<QDimThetaMesh>()),
                KOKKOS_LAMBDA(int r_thread_index, int theta_thread_index, double& sum) {
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
                    sum += weak_integral_element(
                            idx_test_polar,
                            idx_trial_polar,
                            idx_quad,
                            coeff_alpha,
                            coeff_beta,
                            evaluator,
                            mapping,
                            int_volume);
                },
                result);
        return result;
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
};
