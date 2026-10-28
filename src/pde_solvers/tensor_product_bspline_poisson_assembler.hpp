#pragma once

#include <ginkgo/ginkgo.hpp>

#include "circular_to_cartesian.hpp"
#include "gauss_legendre_integration.hpp"
#include "identity_coordinate_change.hpp"
#include "metric_tensor_evaluator.hpp"
#include "polar_spline_fem_poisson_like_solver.hpp"
#include "spline_interpolation.hpp"
#include "volume_quadrature_nd.hpp"


namespace detail {

/**
 * @brief Computes the value and the gradient of a tensor-product basis function from its 1D factors.
 *
 * @param[out] value
 *      The value of the basis function.
 * @param[out] derivs
 *      The gradient of the basis function (covariant components).
 * @param[in] r_basis
 *      The value and derivative of the radial factor.
 * @param[in] theta_basis
 *      The value and derivative of the poloidal factor.
 */
template <class R_cov, class Theta_cov, class EvalDeriv1DType>
KOKKOS_INLINE_FUNCTION void get_value_and_gradient(
        double& value,
        DVector<R_cov, Theta_cov>& derivs,
        EvalDeriv1DType const& r_basis,
        EvalDeriv1DType const& theta_basis)
{
    value = r_basis.value * theta_basis.value;
    ddcHelper::get<R_cov>(derivs) = r_basis.derivative * theta_basis.value;
    ddcHelper::get<Theta_cov>(derivs) = r_basis.value * theta_basis.derivative;
}

/**
 * @brief Computes a quadrature summand corresponding to the
 *        inner product.
 *
 * @tparam ThisPatch
 *      The patch on which the inner product is computed (given explicitly; the other template
 *      parameters are deduced from the arguments).
 *
 * @param[in] idx_r
 *      The index for radial direction.
 * @param[in] idx_theta
 *      The index for poloidal direction
 * @param[in] test_bspline_val_and_deriv
 *      The data structure containing the derivatives over radial and poloidal directions for test space.
 * @param[in] trial_bspline_val_and_deriv
 *      The data structure containing the derivatives over radial and poloidal directions for trial space.
 * @param[in] test_bspline_val_and_deriv_theta
 *       The data structure containing the value and derivative along poloidal direction for test space.
 * @param[in] trial_bspline_val_and_deriv_theta
 *       The data structure containing the value and derivative along poloidal direction for trial space.
 * @param[in] coeff_alpha
 *      The spline representation of the @f$ \alpha @f$ function in the
 *      definition of the Poisson-like equation.
 * @param[in] coeff_beta
 *      The spline representation of the  @f$ \beta @f$ function in the
 *      definition of the Poisson-like equation.
 * @param[in] spline_evaluator
 *      An evaluator for evaluating 2D splines on @f$(r,\theta)@f$.
 * @param[in] mapping
 *      The mapping from the logical index range to the physical index range where
 *      the equation is defined.
 * @param[in] int_volume
 *      The integral volume associated with each point used in the quadrature scheme.
 * @return
  inner product of the test and trial spline is computed using a
 * quadrature. This function returns one summand of the quadrature for
 * the quadrature point given by the indices.
 */
template <
        class ThisPatch,
        class QDimRMesh,
        class QDimThetaMesh,
        class EvalDeriv1DType,
        class SplineRThetaEvaluatorNullBound,
        class Mapping>
KOKKOS_FUNCTION double weak_integral_element(
        Idx<QDimRMesh> idx_r,
        Idx<QDimThetaMesh> idx_theta,
        EvalDeriv1DType const& test_bspline_val_and_deriv,
        EvalDeriv1DType const& trial_bspline_val_and_deriv,
        EvalDeriv1DType const& test_bspline_val_and_deriv_theta,
        EvalDeriv1DType const& trial_bspline_val_and_deriv_theta,
        DConstField<IdxRange<typename ThisPatch::BSplines1, typename ThisPatch::BSplines2>>
                coeff_alpha,
        DConstField<IdxRange<typename ThisPatch::BSplines1, typename ThisPatch::BSplines2>>
                coeff_beta,
        SplineRThetaEvaluatorNullBound const& spline_evaluator,
        Mapping const& mapping,
        DField<IdxRange<QDimRMesh, QDimThetaMesh>> int_volume)
{
    using R = typename ThisPatch::Dim1;
    using Theta = typename ThisPatch::Dim2;
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;
    using CoordRTheta = Coord<R, Theta>;

    // Calculate coefficients at quadrature point
    CoordRTheta coord(ddc::coordinate(idx_r), ddc::coordinate(idx_theta));
    const double alpha = spline_evaluator(coord, coeff_alpha);
    const double beta = spline_evaluator(coord, coeff_beta);
    // Define the value and gradient of the test and trial basis functions
    double basis_val_test_space;
    double basis_val_trial_space;
    DVector<R_cov, Theta_cov> basis_derivs_test_space;
    DVector<R_cov, Theta_cov> basis_derivs_trial_space;
    get_value_and_gradient(
            basis_val_test_space,
            basis_derivs_test_space,
            test_bspline_val_and_deriv,
            test_bspline_val_and_deriv_theta);
    get_value_and_gradient(
            basis_val_trial_space,
            basis_derivs_trial_space,
            trial_bspline_val_and_deriv,
            trial_bspline_val_and_deriv_theta);

    MetricTensorEvaluator<Mapping, CoordRTheta> get_metric_tensor(mapping);

    Tensor inv_metric_tensor = get_metric_tensor.inverse(coord);

    // Assemble the weak integral element
    return int_volume(idx_r, idx_theta)
           * (alpha
                      * tensor_mul(
                              index<'i'>(basis_derivs_test_space),
                              index<'i', 'j'>(inv_metric_tensor),
                              index<'j'>(basis_derivs_trial_space))
              + beta * basis_val_test_space * basis_val_trial_space);
}

} // namespace detail



template <typename ThisPatch, typename X, typename Y>
class TensorProductBSplinePoissonAssembler
{
    using MatrixType = gko::matrix::Csr<double>;

    using GridR = typename ThisPatch::Grid1;
    using GridTheta = typename ThisPatch::Grid2;
    using BSplinesR = typename ThisPatch::BSplines1;
    using BSplinesTheta = typename ThisPatch::BSplines2; //MARK3
    /// Whether the poloidal direction is periodic; given by the poloidal B-splines.
    static constexpr bool s_theta_is_periodic = BSplinesTheta::is_periodic();
    using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;
    using IdxRangeBSR = IdxRange<BSplinesR>;
    using IdxRangeBSTheta = IdxRange<BSplinesTheta>;

    using IdxBSRTheta = Idx<BSplinesR, BSplinesTheta>;

    using ConstSpline2D = DConstField<IdxRangeBSRTheta>;

    using SplineCoeffMem = DFieldMem<IdxRangeBSRTheta>;
    using SplineCoeff = DField<IdxRangeBSRTheta>;

    using KnotsR = ddc::NonUniformBsplinesKnots<BSplinesR>;
    using KnotsTheta = ddc::NonUniformBsplinesKnots<BSplinesTheta>;

    using R = typename ThisPatch::Dim1;
    using Theta = typename ThisPatch::Dim2;

    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

    using CoordRTheta = Coord<R, Theta>;


    /**
     * @brief Tag the first dimension for the quadrature mesh.
     */
    struct QDimRMesh : NonUniformGridBase<R>
    {
    };
    /**
     * @brief Tag the second dimension for the quadrature mesh.
     */
    struct QDimThetaMesh : NonUniformGridBase<Theta>
    {
    };
    static constexpr int s_n_gauss_legendre_r = BSplinesR::degree() + 1;
    static constexpr int s_n_gauss_legendre_theta = BSplinesTheta::degree() + 1;

    using IdxRangeQuadratureR = IdxRange<QDimRMesh>;
    /**
     * @brief Tag the quadrature index range in the second dimension.
     */
    using IdxRangeQuadratureTheta = IdxRange<QDimThetaMesh>;

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

    using SplineRThetaEvaluatorNullBound = std::conditional_t<
            s_theta_is_periodic,
            ddc::SplineEvaluator2D<
                    Kokkos::DefaultExecutionSpace,
                    typename Kokkos::DefaultExecutionSpace::memory_space,
                    BSplinesR,
                    BSplinesTheta,
                    GridR,
                    GridTheta,
                    ddc::NullExtrapolationRule, // boundary at r=0
                    ddc::NullExtrapolationRule, // boundary at rmax
                    ddc::PeriodicExtrapolationRule<Theta>,
                    ddc::PeriodicExtrapolationRule<Theta>>,
            ddc::SplineEvaluator2D<
                    Kokkos::DefaultExecutionSpace,
                    typename Kokkos::DefaultExecutionSpace::memory_space,
                    BSplinesR,
                    BSplinesTheta, // MARK1
                    GridR,
                    GridTheta,
                    ddc::NullExtrapolationRule, // boundary at r=0
                    ddc::NullExtrapolationRule, // boundary at rmax
                    ddc::NullExtrapolationRule,
                    ddc::NullExtrapolationRule>>;

    const int m_nbasis_r;
    const int m_nbasis_theta;
    const int m_matrix_size;

    // Domains
    IdxRangeBSR m_idxrange_bsplines_r;
    IdxRangeBSTheta m_idxrange_bsplines_theta;
    IdxRangeBSRTheta m_idxrange_bsplines;
    IdxRangeBSRTheta m_idxrange_bsplines_with_reps;

    IdxRangeQuadratureR m_idxrange_quadrature_r;
    IdxRangeQuadratureTheta m_idxrange_quadrature_theta;

    // Gauss-Legendre points and weights
    FieldMem<double, IdxRangeQuadratureR> m_weights_r;
    FieldMem<double, IdxRangeQuadratureTheta> m_weights_theta;


    FieldMem<double, IdxRangeQuadratureRTheta> m_int_volume;


    struct RBasisSubset
    {
    };
    struct ThetaBasisSubset
    {
    };
    struct RCellDim
    {
    };
    struct ThetaCellDim
    {
    };
    using IdxCell = Idx<RCellDim, ThetaCellDim>;

    /**
    * @brief Object storing a value and a value of the derivative
    * of a 1D function.
    */
    struct EvalDeriv1DType
    {
        double value;
        double derivative;
    };

    // Basis Spline values and derivatives at Gauss-Legendre points
    host_t<FieldMem<EvalDeriv1DType, IdxRange<RBasisSubset, QDimRMesh>>> m_r_basis_vals_and_derivs;
    host_t<FieldMem<EvalDeriv1DType, IdxRange<ThetaBasisSubset, QDimThetaMesh>>>
            m_theta_basis_vals_and_derivs;


    // Number of cells over which a radial B-splines has its support
    // This is the case for b-splines which are not affected by the higher knot multiplicity at the boundary.
    static constexpr IdxStep<RBasisSubset> m_n_non_zero_bases_r
            = IdxStep<RBasisSubset>(BSplinesR::degree() + 1);

    // Number of cells over which a poloidal B-splines has its support
    static constexpr IdxStep<ThetaBasisSubset> m_n_non_zero_bases_theta
            = IdxStep<ThetaBasisSubset>(BSplinesTheta::degree() + 1);

    static constexpr IdxRange<RBasisSubset> m_non_zero_bases_r
            = IdxRange<RBasisSubset>(Idx<RBasisSubset> {0}, m_n_non_zero_bases_r);
    static constexpr IdxRange<ThetaBasisSubset> m_non_zero_bases_theta
            = IdxRange<ThetaBasisSubset>(Idx<ThetaBasisSubset> {0}, m_n_non_zero_bases_theta);

    const int m_batch_idx {0}; // TODO: Remove when batching is supported

    using Mapping = std::conditional_t<
            s_theta_is_periodic,
            CircularToCartesian<typename ThisPatch::Dim1, typename ThisPatch::Dim2, X, Y>,
            IdentityCoordinateChange<
                    VectorIndexSet<typename ThisPatch::Dim1, typename ThisPatch::Dim2>,
                    VectorIndexSet<X, Y>>>;
    Mapping m_mapping;

public:
    TensorProductBSplinePoissonAssembler()
        : m_nbasis_r(ddc::discrete_space<BSplinesR>().nbasis())
        , m_nbasis_theta(ddc::discrete_space<BSplinesTheta>().nbasis())
        , m_matrix_size(m_nbasis_r * m_nbasis_theta)
        , m_idxrange_bsplines_r(ddc::discrete_space<BSplinesR>().full_domain())
        , m_idxrange_bsplines_theta(get_idx_range_bsplines_theta())
        , m_idxrange_bsplines(m_idxrange_bsplines_r, m_idxrange_bsplines_theta)
        , m_idxrange_bsplines_with_reps(
                  m_idxrange_bsplines_r,
                  ddc::discrete_space<BSplinesTheta>().full_domain())
        , m_idxrange_quadrature_r(
                  Idx<QDimRMesh>(0),
                  IdxStep<QDimRMesh>(
                          s_n_gauss_legendre_r * ddc::discrete_space<BSplinesR>().ncells()))
        , m_idxrange_quadrature_theta(
                  Idx<QDimThetaMesh>(0),
                  IdxStep<QDimThetaMesh>(
                          s_n_gauss_legendre_theta * ddc::discrete_space<BSplinesTheta>().ncells()))
        , m_weights_r(m_idxrange_quadrature_r)
        , m_weights_theta(m_idxrange_quadrature_theta)
        , m_r_basis_vals_and_derivs(
                  IdxRange<RBasisSubset, QDimRMesh>(m_non_zero_bases_r, m_idxrange_quadrature_r))
        , m_theta_basis_vals_and_derivs(
                  IdxRange<
                          ThetaBasisSubset,
                          QDimThetaMesh>(m_non_zero_bases_theta, m_idxrange_quadrature_theta))
    {
        ////// DUPLICATED FROM CONSTRUCTOR OF POLAR POISSON ////////
        IdxRange<KnotsR> idxrange_r1_break_points
                = ddc::discrete_space<BSplinesR>().break_point_domain();
        IdxRange<KnotsTheta> idxrange_theta_edges
                = ddc::discrete_space<BSplinesTheta>().break_point_domain();
        host_t<FieldMem<Coord<typename ThisPatch::Dim1>, IdxRange<KnotsR>>> breaks_r(
                idxrange_r1_break_points);
        host_t<FieldMem<Coord<typename ThisPatch::Dim2>, IdxRange<KnotsTheta>>> breaks_theta(
                idxrange_theta_edges);

        ddcHelper::dump_coordinates(Kokkos::DefaultHostExecutionSpace(), get_field(breaks_r));
        ddcHelper::dump_coordinates(Kokkos::DefaultHostExecutionSpace(), get_field(breaks_theta));

        // Define quadrature points and weights
        GaussLegendre<QDimRMesh, s_n_gauss_legendre_r> gl_coeffs_r(get_const_field(breaks_r));
        GaussLegendre<QDimThetaMesh, s_n_gauss_legendre_theta> gl_coeffs_theta(
                get_const_field(breaks_theta));
        m_int_volume = compute_coeffs_on_mapping(
                Kokkos::DefaultExecutionSpace(),
                m_mapping,
                gauss_legendre_quadrature_coefficients<
                        Kokkos::DefaultExecutionSpace>(gl_coeffs_r, gl_coeffs_theta));
        // Find value and derivative of 1D B-splines in radial direction
        ddc::host_for_each(m_idxrange_quadrature_r, [&](IdxQuadratureR const idx_r) {
            std::array<double, 2 * m_n_non_zero_bases_r> data;
            DSpan2D vals(data.data(), m_n_non_zero_bases_r, 2);
            ddc::discrete_space<BSplinesR>()
                    .eval_basis_and_n_derivs(vals, ddc::coordinate(idx_r), 1);
            for (auto ib : m_non_zero_bases_r) {
                const int ib_idx = ib - m_non_zero_bases_r.front();
                m_r_basis_vals_and_derivs(ib, idx_r).value = vals(ib_idx, 0);
                m_r_basis_vals_and_derivs(ib, idx_r).derivative = vals(ib_idx, 1);
            }
        });

        // Find value and derivative of 1D B-splines in poloidal direction
        ddc::host_for_each(m_idxrange_quadrature_theta, [&](IdxQuadratureTheta const idx_theta) {
            std::array<double, 2 * m_n_non_zero_bases_theta> data;
            DSpan2D vals(data.data(), m_n_non_zero_bases_theta, 2);
            ddc::discrete_space<BSplinesTheta>()
                    .eval_basis_and_n_derivs(vals, ddc::coordinate(idx_theta), 1);
            for (auto ib : m_non_zero_bases_theta) {
                const int ib_idx = ib - m_non_zero_bases_theta.front();
                m_theta_basis_vals_and_derivs(ib, idx_theta).value = vals(ib_idx, 0);
                m_theta_basis_vals_and_derivs(ib, idx_theta).derivative = vals(ib_idx, 1);
            }
        });
    }

    void operator()(
            Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>& vals_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace>& cols_csr,
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace>& nnz_per_row_csr_host)
    {
        Kokkos::resize(nnz_per_row_csr_host, m_matrix_size + 2);

        Kokkos::Profiling::pushRegion("PolarPoissonInitNnz");

        size_t constexpr degree = BSplinesR::degree();
        std::size_t const n_splines_1d_overlap = 2 * degree + 1;

        // TODO: Add Comments!
        auto compute_nnz_per_row = [&](int j) {
            int min_idx = std::max(j - 3, 0);
            int max_idx = std::min(j + 3, m_nbasis_theta - 1);
            return max_idx + 1 - min_idx;
        };
        auto nnz_per_row_for_thetas = [&](int i) {
            for (int j = 0; j < m_nbasis_theta; j++) {
                const int row_offset = i * m_nbasis_theta;
                nnz_per_row_csr_host(row_offset + j + 2)
                        = (degree + i + 1) * compute_nnz_per_row(j);
            }
        };

        // Fill the nnz_per_row array for the first couple of r values
        for (int i(0); i < degree; i++) {
            if constexpr (s_theta_is_periodic) {
                for (int k(i * m_nbasis_theta); k < (i + 1) * m_nbasis_theta; k++) {
                    nnz_per_row_csr_host(k + 2) = (degree + i + 1) * n_splines_1d_overlap;
                }
            } else {
                for (int j = 0; j < m_nbasis_theta; j++) {
                    const int row_offset = i * m_nbasis_theta;
                    nnz_per_row_csr_host(row_offset + j + 2)
                            = (degree + i + 1) * compute_nnz_per_row(j);
                }
            }
        }
        if constexpr (s_theta_is_periodic) {
            for (int i(degree * m_nbasis_theta); i < m_matrix_size - degree * m_nbasis_theta; i++) {
                nnz_per_row_csr_host(i + 2) = n_splines_1d_overlap * n_splines_1d_overlap;
            }
        } else {
            for (int i = degree; i < m_nbasis_r - degree; i++) {
                for (int j = 0; j < m_nbasis_theta; j++) {
                    const int row_offset = i * m_nbasis_theta;
                    nnz_per_row_csr_host(row_offset + j + 2)
                            = n_splines_1d_overlap * compute_nnz_per_row(j);
                }
            }
        }

        for (int i(1); i < degree + 1; i++) {
            if constexpr (s_theta_is_periodic)
                for (size_t k(m_matrix_size - i * m_nbasis_theta);
                     k < m_matrix_size - (i - 1) * m_nbasis_theta;
                     k++) {
                    nnz_per_row_csr_host(k + 2) = (degree + i) * n_splines_1d_overlap;
                }
            else
                for (int j = 0; j < m_nbasis_theta; j++) {
                    const int row_offset = (m_nbasis_theta - i) * m_nbasis_theta;
                    nnz_per_row_csr_host(row_offset + j + 2)
                            = (degree + i) * compute_nnz_per_row(j);
                }
        }


        // sum non-zero elements count
        for (size_t k = 2; k < m_matrix_size + 2; k++) {
            nnz_per_row_csr_host(k) += nnz_per_row_csr_host(k - 1);
        }
        nnz_per_row_csr_host(0) = 0;
        nnz_per_row_csr_host(1) = 0;

        Kokkos::Profiling::popRegion();

        SplineCoeffMem coeff_alpha_mem(m_idxrange_bsplines_with_reps);
        SplineCoeffMem coeff_beta_mem(m_idxrange_bsplines_with_reps);

        SplineCoeff coeff_alpha = get_field(coeff_alpha_mem);
        SplineCoeff coeff_beta = get_field(coeff_beta_mem);

        // TODO: Let coefficients be given as parameter
        ddc::parallel_fill(coeff_alpha, 1.0);
        ddc::parallel_fill(coeff_beta, 0.0);

        ddc::NullExtrapolationRule bv_r_min;
        ddc::NullExtrapolationRule bv_r_max;
        auto init_spline_evaluator = [&]() {
            if constexpr (s_theta_is_periodic) {
                return SplineRThetaEvaluatorNullBound(
                        bv_r_min,
                        bv_r_max,
                        ddc::PeriodicExtrapolationRule<Theta> {},
                        ddc::PeriodicExtrapolationRule<Theta> {});
            } else
                return SplineRThetaEvaluatorNullBound(bv_r_min, bv_r_max, bv_r_min, bv_r_max);
        };
        SplineRThetaEvaluatorNullBound spline_evaluator = init_spline_evaluator(); // MARK0

        int const nnz_total = nnz_per_row_csr_host(nnz_per_row_csr_host.size() - 1);

        Kokkos::resize(vals_csr, 1, nnz_total);
        Kokkos::resize(cols_csr, nnz_total);

        // Calculate the matrix elements following a stencil
        ddc::host_for_each(m_idxrange_bsplines, [&](IdxBSRTheta const idx_test) {
            const std::size_t idx_test_r(ddc::select<BSplinesR>(idx_test).uid());
            const std::size_t idx_test_theta(ddc::select<BSplinesTheta>(idx_test).uid());
            int n_splines_theta = ddc::get<BSplinesTheta>(m_idxrange_bsplines_theta.extents());
            int int_idx_test = n_splines_theta * idx_test_r + idx_test_theta;
            // The idea is that we only fill the indexes larger or equal than the curent index
            // The other ones are dealt with using the symmetry
            // These are the theta indices that have remain to be filled
            IdxRangeBSTheta remaining_theta;
            if constexpr (s_theta_is_periodic) {
                remaining_theta = IdxRangeBSTheta(
                        Idx<BSplinesTheta> {idx_test_theta},
                        IdxStep<BSplinesTheta> {BSplinesTheta::degree() + 1});
            } else {
                remaining_theta = IdxRangeBSTheta(
                        Idx<BSplinesTheta> {idx_test_theta},
                        IdxStep<BSplinesTheta> {
                                std::
                                        min(static_cast<int>(BSplinesTheta::degree() + 1),
                                            static_cast<int>(m_nbasis_theta)
                                                    - static_cast<int>(idx_test_theta))});
            }
            // Fill the entries with the same r-index first
            ddc::host_for_each(remaining_theta, [&](Idx<BSplinesTheta> const idx_trial_theta) {
                IdxBSRTheta idx_trial(Idx<BSplinesR>(idx_test_r), idx_trial_theta);
                double element = get_matrix_stencil_element(
                        idx_test,
                        idx_trial,
                        get_const_field(coeff_alpha),
                        get_const_field(coeff_beta),
                        spline_evaluator,
                        m_mapping);
                if (idx_test == idx_trial) {
                    //

                    const int idx = nnz_per_row_csr_host(int_idx_test + 1);
                    assert("Idx out of bounds" && idx >= 0 && idx < cols_csr.size());
                    cols_csr(idx) = int_idx_test;
                    assert("Idx out of bounds" && idx >= 0 && idx < vals_csr.size());
                    vals_csr(m_batch_idx, idx) = element;
                    assert("Idx out of bounds" && int_idx_test + 1 >= 0
                           && int_idx_test + 1 < nnz_per_row_csr_host.size());
                    nnz_per_row_csr_host(int_idx_test + 1)++;

                } else {
                    const std::size_t idx_trial_r(ddc::select<BSplinesR>(idx_trial).uid());
                    std::size_t idx_trial_theta_ul;
                    if constexpr (s_theta_is_periodic) {
                        idx_trial_theta_ul = theta_mod(ddc::select<BSplinesTheta>(idx_trial).uid());
                    } else
                        idx_trial_theta_ul = ddc::select<BSplinesTheta>(idx_trial).uid();
                    const int int_idx_trial = idx_trial_r * n_splines_theta + idx_trial_theta_ul;
                    const int aij_idx = nnz_per_row_csr_host(int_idx_test + 1);
                    assert("Idx out of bounds" && aij_idx >= 0 && aij_idx < vals_csr.size());
                    assert("Idx out of bounds" && aij_idx >= 0 && aij_idx < cols_csr.size());
                    assert("Idx out of bounds" && int_idx_test + 1 >= 0
                           && int_idx_test + 1 < nnz_per_row_csr_host.size());

                    cols_csr(aij_idx) = int_idx_trial;
                    vals_csr(m_batch_idx, aij_idx) = element;
                    nnz_per_row_csr_host(int_idx_test + 1)++;


                    const int aji_idx = nnz_per_row_csr_host(int_idx_trial + 1);
                    assert("Idx out of bounds" && aji_idx >= 0 && aji_idx < vals_csr.size());
                    assert("Idx out of bounds" && aji_idx >= 0 && aji_idx < cols_csr.size());
                    assert("Idx out of bounds" && int_idx_trial + 1 >= 0
                           && int_idx_trial + 1 < nnz_per_row_csr_host.size());

                    cols_csr(aji_idx) = int_idx_test;
                    vals_csr(m_batch_idx, aji_idx) = element;
                    nnz_per_row_csr_host(int_idx_trial + 1)++;
                }
            });
            // Due to symmetry we only have to take care of the r indices larger than the test index
            IdxRangeBSR remaining_r(
                    ddc::select<BSplinesR>(idx_test) + 1,
                    IdxStep<BSplinesR> {std::
                                                min(BSplinesR::degree(),
                                                    ddc::discrete_space<BSplinesR>().nbasis() - 1
                                                            - idx_test_r)});
            // Now we need the relevant theta to the left and the right
            // We add the whole N_theta. That can not work in the non-periodic context
            IdxRangeBSTheta relevant_theta;
            if constexpr (s_theta_is_periodic) {
                relevant_theta = IdxRangeBSTheta(
                        Idx<BSplinesTheta> {
                                idx_test_theta + ddc::discrete_space<BSplinesTheta>().nbasis()
                                - BSplinesTheta::degree()},
                        IdxStep<BSplinesTheta> {2 * BSplinesTheta::degree() + 1});
            } else {
                Idx<BSplinesTheta> const idx_min {
                        std::max(0, static_cast<int>(idx_test_theta) - 3)};
                Idx<BSplinesTheta> const idx_max {
                        std::min(m_nbasis_theta - 1, static_cast<int>(idx_test_theta) + 3)};
                relevant_theta
                        = IdxRangeBSTheta(Idx<BSplinesTheta> {idx_min}, idx_max + 1 - idx_min);
            }
            IdxRangeBSRTheta trial_idx_range(remaining_r, relevant_theta);
            ddc::host_for_each(trial_idx_range, [&](IdxBSRTheta const idx_trial) {
                const int idx_trial_r(ddc::select<BSplinesR>(idx_trial).uid());
                int idx_trial_theta;
                if constexpr (s_theta_is_periodic)
                    idx_trial_theta = theta_mod(ddc::select<BSplinesTheta>(idx_trial).uid());
                else
                    idx_trial_theta = ddc::select<BSplinesTheta>(idx_trial).uid();
                int int_idx_trial = n_splines_theta * idx_trial_r + idx_trial_theta;

                double element = get_matrix_stencil_element(
                        idx_test,
                        idx_trial,
                        get_const_field(coeff_alpha),
                        get_const_field(coeff_beta),
                        spline_evaluator,
                        m_mapping);


                if (idx_test == idx_trial) {
                    const int idx = nnz_per_row_csr_host(int_idx_test + 1);
                    cols_csr(idx) = int_idx_test;
                    vals_csr(m_batch_idx, idx) = element;
                    nnz_per_row_csr_host(int_idx_test + 1)++;
                } else {
                    const int aij_idx = nnz_per_row_csr_host(int_idx_test + 1);
                    assert("Idx out of bounds" && aij_idx >= 0 && aij_idx < vals_csr.size());
                    assert("Idx out of bounds" && aij_idx >= 0 && aij_idx < cols_csr.size());
                    assert("Idx out of bounds" && int_idx_test + 1 >= 0
                           && int_idx_test + 1 < nnz_per_row_csr_host.size());

                    cols_csr(aij_idx) = int_idx_trial;
                    vals_csr(m_batch_idx, aij_idx) = element;
                    nnz_per_row_csr_host(int_idx_test + 1)++;
                    const int aji_idx = nnz_per_row_csr_host(int_idx_trial + 1);

                    assert("Idx out of bounds" && aji_idx >= 0 && aji_idx < vals_csr.size());
                    assert("Idx out of bounds" && aji_idx >= 0 && aji_idx < cols_csr.size());
                    assert("Idx out of bounds" && int_idx_trial + 1 >= 0
                           && int_idx_trial + 1 < nnz_per_row_csr_host.size());
                    cols_csr(aji_idx) = int_idx_test;
                    vals_csr(m_batch_idx, aji_idx) = element;
                    nnz_per_row_csr_host(int_idx_trial + 1)
                            = nnz_per_row_csr_host(int_idx_trial + 1) + 1;
                }
            });
        });
    }

    /**
     * @brief Computes the matrix element corresponding to two tensor product splines
     *        with index idx_test and idx_trial
     *
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
     *      The mapping from the logical index range to the physical index range where
     *      the equation is defined.
     * @return 
     *      The value of the matrix element.
     */
    template <class Mapping>
    double get_matrix_stencil_element(
            IdxBSRTheta idx_test,
            IdxBSRTheta idx_trial,
            ConstSpline2D coeff_alpha,
            ConstSpline2D coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        // 0 <= idx_test_r < 8
        // 0 <= idx_trial_r < 8
        // idx_test_r < idx_trial_r
        const int idx_test_r(ddc::select<BSplinesR>(idx_test).uid());
        const int idx_trial_r(ddc::select<BSplinesR>(idx_trial).uid());
        // 0 <= idx_test_theta < 8
        // 0 <= idx_trial_theta < 8
        int idx_test_theta;
        int idx_trial_theta;
        if constexpr (s_theta_is_periodic) {
            idx_test_theta = theta_mod(ddc::select<BSplinesTheta>(idx_test).uid());
            idx_trial_theta = theta_mod(ddc::select<BSplinesTheta>(idx_trial).uid());
        } else {
            idx_test_theta = ddc::select<BSplinesTheta>(idx_test).uid();
            idx_trial_theta = ddc::select<BSplinesTheta>(idx_trial).uid();
        }
        const std::size_t ncells_r = ddc::discrete_space<BSplinesR>().ncells();

        // 0<= r_offset <= degree_r
        // -degree_theta <= theta_offset <= degree_theta
        const int r_offset = idx_trial_r - idx_test_r;
        // Added the comment below. Examples with 10 periodic θ splines and degree 3:
        //      test 9, trial 0: theta_mod(0 - 9) = 1 < 7, so theta_offset = 1 (trial is the right
        //      neighbour across θ = 0); test 0, trial 9: theta_mod(9 - 0) = 9 >= 7, so 9 - 10 = -1.
        // theta_offset = idx_trial_theta - idx_test_theta. With periodicity, the difference is taken
        // modulo nbasis_theta and shifted into [-degree, nbasis_theta - degree - 1], so that two splines
        // that are neighbours across θ = 0 get a small offset instead of one close to ±nbasis_theta.

        // This is needed because the overlap formulas below need the signed distance between the two splines
        // along θ, as it is on the circle: the sign says which spline is to the right (and so where the
        // overlap starts), and |offset| <= p says how many cells they share (n = p + 1 - |offset|).
        // The raw modulo value does not give this across θ = 0. Example with 10 periodic θ splines,
        // p = 3: test 0, trial 9. theta_mod(9 - 0) = 9, which would give n = p + 1 - 9 = -5 cells and the
        // asserts below would fail. Geometrically, trial 9 is the left neighbour of test 0, so the
        // offset must be -1: n = p + 1 - 1 = 3, first = test - p = -3, i.e. cells 7, 8, 9 after the
        // modulo. The threshold nbasis_theta - p is where "far to the right" turns into "close on the
        // left": an offset in [nbasis - p, nbasis - 1] corresponds to [-p, -1].
        int theta_offset;
        if constexpr (s_theta_is_periodic) {
            theta_offset = theta_mod(idx_trial_theta - idx_test_theta);
            if (theta_offset >= int(m_nbasis_theta - BSplinesTheta::degree())) {
                theta_offset -= m_nbasis_theta;
            }
        } else {
            theta_offset = idx_trial_theta - idx_test_theta;
        }
        assert(r_offset >= 0);
        assert(r_offset <= int(BSplinesR::degree()));
        assert(theta_offset >= -int(BSplinesTheta::degree()));
        assert(theta_offset <= int(BSplinesTheta::degree()));
        // Find the index range covering the cells where both the test and trial functions are non-zero
        //
        // Index convention (B-splines of degree p, cells numbered 0, ..., ncells - 1):
        // the 1D B-spline with index i is non-zero on the p + 1 cells i - p, ..., i.
        // Without periodicity, cells outside [0, ncells - 1] do not exist, so the first p and the
        // last p B-splines have shorter supports. With periodicity, cell indices are taken modulo ncells.
        //
        // In each dimension (r and θ separately) the supports of the test spline (index i) and of the
        // trial spline (index j = i + offset, |offset| <= p) intersect in a range of consecutive cell
        // indices, described by
        //   first_overlap_*     : the first cell index of the intersection,
        //   n_overlap_stencil_* : the number of cells in the intersection, before cutting it to the
        //                         cells that exist.
        // With periodic θ these indices are taken modulo ncells, so the range can run past the last
        // cell and continue at cell 0. In the parametric domain the overlap then consists of two pieces,
        // one at each end of the θ interval (e.g. cells 8, 9 and 0 for 10 cells). The quadrature is still
        // correct because each cell is integrated separately: the cell loop below maps every index with
        // `theta_mod` to the cell that exists.
        // For offset >= 0: [i - p, i] and [j - p, j] intersect on [j - p, i]:
        //                  first = j - p,  n = i - (j - p) + 1 = p + 1 - offset.
        // For offset <  0: they intersect on [i - p, j]:
        //                  first = i - p,  n = j - (i - p) + 1 = p + 1 + offset.
        // Afterwards the range is cut to the existing cells (below), which is only needed without
        // periodicity: in r always, in θ only for non-periodic θ.
        int n_overlap_stencil_r(BSplinesR::degree() + 1 - r_offset);
        // (index convention and all first_overlap_* / n_overlap_stencil_* variables) and the line below. In r only trial >= test is assembled
        //      (r_offset >= 0, the rest follows from symmetry), so this is always the "offset >= 0" case.
        // r_offset >= 0, so the intersection starts at the first cell of the trial support.
        int first_overlap_r(idx_trial_r - BSplinesR::degree());

        // Index of left-most cell of the overlap
        int first_overlap_theta;
        // Number of overlap cells
        int n_overlap_stencil_theta;

        if (theta_offset > 0) {
            n_overlap_stencil_theta = BSplinesTheta::degree() + 1 - theta_offset;
            // This is the "offset >= 0" case of the comment block above (here offset > 0):
            //      test i is non-zero on cells [i - p, i], trial j = i + offset on [j - p, j]; since
            //      j > i they intersect on [j - p, i], i.e. first = j - p and n = p + 1 - offset cells.
            //      Example (p = 3): test 5, trial 6: cells [3, 6] cap [2, 5] = [3, 5], first = 3, n = 3.
            //      With periodicity first is taken modulo the number of cells.
            if constexpr (s_theta_is_periodic) {
                first_overlap_theta = theta_mod(idx_trial_theta - BSplinesTheta::degree());
            } else
                first_overlap_theta = idx_trial_theta - BSplinesTheta::degree();
            // The same formula as in the periodic branch, first = j - p, but without the
            //      modulo. It can therefore be negative (j < p) or the range can end beyond the
            //      last cell; both are cut to the existing cells further down (the F1 fix).
        } else {
            n_overlap_stencil_theta = BSplinesTheta::degree() + 1 + theta_offset;
            // first_overlap_* is the first cell where both splines are non-zero, n_overlap_stencil_* the number of such
            //      cells. This branch is the "offset <= 0" case: the trial spline j is left of (or equal
            //      to) the test spline i, so the supports [i - p, i] and [j - p, j] intersect on
            //      [i - p, j]: first = i - p, n = p + 1 + offset.
            //      Example (p = 3): test 6, trial 5 (offset -1): [3, 6] cap [2, 5] = [3, 5], first = 3, n = 3.
            //      offset = 0 (same spline): first = i - p, n = p + 1, the whole support.
            if constexpr (s_theta_is_periodic) {
                first_overlap_theta = theta_mod(idx_test_theta - BSplinesTheta::degree());
            } else
                first_overlap_theta = idx_test_theta - BSplinesTheta::degree();
        }

        if (first_overlap_r < 0) {
            // r is not periodic, so cells with negative index do not exist. If the
            //      intersection starts at first_overlap_r < 0, its first -first_overlap_r cells are
            //      dropped: the count shrinks by that amount (n_compact is negative) and the range starts
            //      at cell 0. Example (p = 3): test 0, trial 1: [-3, 0] cap [-2, 1] = [-2, 0], first = -2,
            //      n = 3; only cell 0 exists, so n = 3 + (-2) = 1, first = 0. The other end (cells beyond
            //      ncells - 1) is cut by `n_to_edge_r` below.
            const int n_compact = first_overlap_r;
            first_overlap_r = 0;
            n_overlap_stencil_r += n_compact;
        }
        if constexpr (!s_theta_is_periodic) {
            // Without periodicity the overlap has to be restricted to the domain, as in r:
            // - start: cells with negative index do not exist; drop them (count shrinks, start at 0);
            // - end: cells beyond ncells - 1 do not exist; keep at most ncells - first_overlap_theta cells.
            // Without this cut, theta_mod in the cell loop below would wrap the non-existing cells onto
            // the other end of the θ interval and add integrals over cells where the splines are zero.
            // Examples with 10 non-periodic θ cells, degree 3:
            //      - test 1, trial 2: [-2, 1] cap [-1, 2] = [-1, 1], first = -1, n = 3. Cell -1 does not
            //        exist; before the fix theta_mod mapped it to cell 9. Now first = 0, n = 2 (cells 0, 1).
            //      - test 11, trial 12: [8, 11] cap [9, 12] = [9, 11], first = 9, n = 3. Cells 10, 11 do
            //        not exist; before the fix they were mapped to cells 0, 1. Now n = min(3, 10 - 9) = 1.
            //      This is exactly why only the rows of θ splines 0-2 and 10-12 were wrong.
            if (first_overlap_theta < 0) {
                n_overlap_stencil_theta += first_overlap_theta;
                first_overlap_theta = 0;
            }
            const int n_to_edge_theta(
                    ddc::discrete_space<BSplinesTheta>().ncells() - first_overlap_theta);
            n_overlap_stencil_theta = std::min(n_overlap_stencil_theta, n_to_edge_theta);
        }

        const int n_to_edge_r(ncells_r - first_overlap_r);
        const IdxStep<RCellDim> n_overlap_r(std::min(n_overlap_stencil_r, n_to_edge_r));
        const IdxStep<ThetaCellDim> n_overlap_theta(n_overlap_stencil_theta);


        const Idx<RCellDim> first_overlap_element_r(first_overlap_r);
        const Idx<ThetaCellDim> first_overlap_element_theta(first_overlap_theta);

        const IdxRange<RCellDim> r_cells(first_overlap_element_r, n_overlap_r);
        const IdxRange<ThetaCellDim> theta_cells(first_overlap_element_theta, n_overlap_theta);
        const IdxRange<RCellDim, ThetaCellDim> non_zero_cells(r_cells, theta_cells);

        auto r_basis_vals_and_derivs_alloc = ddc::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(m_r_basis_vals_and_derivs));
        auto theta_basis_vals_and_derivs_alloc = ddc::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(m_theta_basis_vals_and_derivs));

        Field<EvalDeriv1DType, IdxRange<RBasisSubset, QDimRMesh>> r_basis_vals_and_derivs
                = get_field(r_basis_vals_and_derivs_alloc);
        Field<EvalDeriv1DType, IdxRange<ThetaBasisSubset, QDimThetaMesh>>
                theta_basis_vals_and_derivs = get_field(theta_basis_vals_and_derivs_alloc);
        DField<IdxRangeQuadratureRTheta> int_volume_proxy = get_field(m_int_volume);
        assert(n_overlap_r * n_overlap_theta > 0);
        return ddc::host_transform_reduce(
                non_zero_cells,
                0.0,
                ddc::reducer::sum<double>(),
                [&](IdxCell const cell_idx) {
                    const int cell_idx_r(ddc::select<RCellDim>(cell_idx).uid());
                    const int cell_idx_theta(theta_mod(ddc::select<ThetaCellDim>(cell_idx).uid()));
                    const IdxRangeQuadratureRTheta cell_quad_points(
                            get_quadrature_points_in_cell(cell_idx_r, cell_idx_theta));

                    int ib_test_theta_idx = idx_test_theta - cell_idx_theta;
                    int ib_trial_theta_idx = idx_trial_theta - cell_idx_theta;

                    // Find the column where the non-zero data is stored
                    Idx<RBasisSubset> ib_test_r(idx_test_r - cell_idx_r);
                    Idx<ThetaBasisSubset> ib_test_theta(theta_mod(ib_test_theta_idx));
                    Idx<RBasisSubset> ib_trial_r(idx_trial_r - cell_idx_r);
                    Idx<ThetaBasisSubset> ib_trial_theta(theta_mod(ib_trial_theta_idx));

                    assert(ib_test_r.uid() < BSplinesR::degree() + 1);
                    assert(ib_test_theta.uid() < BSplinesTheta::degree() + 1);
                    assert(ib_trial_r.uid() < BSplinesR::degree() + 1);
                    assert(ib_trial_theta.uid() < BSplinesTheta::degree() + 1);

                    // Calculate the weak integral
                    return ddc::parallel_transform_reduce(
                            Kokkos::DefaultExecutionSpace(),
                            cell_quad_points,
                            0.0,
                            ddc::reducer::sum<double>(),
                            [&](IdxQuadratureRTheta const idx_quad) {
                                IdxQuadratureR const idx_r(idx_quad);
                                IdxQuadratureTheta const idx_theta(idx_quad);
                                return detail::weak_integral_element<ThisPatch>(
                                        idx_r,
                                        idx_theta,
                                        r_basis_vals_and_derivs(ib_test_r, idx_r),
                                        r_basis_vals_and_derivs(ib_trial_r, idx_r),
                                        theta_basis_vals_and_derivs(ib_test_theta, idx_theta),
                                        theta_basis_vals_and_derivs(ib_trial_theta, idx_theta),
                                        coeff_alpha,
                                        coeff_beta,
                                        evaluator,
                                        mapping,
                                        int_volume_proxy);
                            });
                });
    }

    /**
     * @brief compute the quadrature range for a given pair of indices
     *
     * @param[in] cell_idx_r
     *      The index for radial direction
     * @param[in] cell_idx_theta
     *      The index for poloidal direction
     * @return 
     *      The quadrature range corresponding to the  @f$(r,\theta)@f$ indices.
     */
    static KOKKOS_FUNCTION IdxRangeQuadratureRTheta
    get_quadrature_points_in_cell(int cell_idx_r, int cell_idx_theta)
    {
        const IdxQuadratureR first_quad_point_r(cell_idx_r * s_n_gauss_legendre_r);
        const IdxQuadratureTheta first_quad_point_theta(cell_idx_theta * s_n_gauss_legendre_theta);
        constexpr IdxStepQuadratureR n_GL_r(s_n_gauss_legendre_r);
        constexpr IdxStepQuadratureTheta n_GL_theta(s_n_gauss_legendre_theta);
        const IdxRangeQuadratureR quad_points_r(first_quad_point_r, n_GL_r);
        const IdxRangeQuadratureTheta quad_points_theta(first_quad_point_theta, n_GL_theta);
        return IdxRangeQuadratureRTheta(quad_points_r, quad_points_theta);
    }


    /**
     * @brief Calculates the modulo idx_theta in relation to cells number along  @f$ \theta @f$ direction .
     *
     * @param[in] idx_theta @f$ \theta @f$ index.
     *
     * @return The corresponding indice modulo @f$ \theta @f$ direction cells number
     */
    static KOKKOS_FUNCTION int theta_mod(int idx_theta)
    {
        int ncells_theta = ddc::discrete_space<BSplinesTheta>().ncells();
        while (idx_theta < 0)
            idx_theta += ncells_theta;
        while (idx_theta >= ncells_theta)
            idx_theta -= ncells_theta;
        return idx_theta;
    }
    auto get_idx_range_bsplines_theta()
    {
        if constexpr (s_theta_is_periodic)
            return ddc::discrete_space<BSplinesTheta>().full_domain().take_first(
                    IdxStep<BSplinesTheta> {m_nbasis_theta});
        else
            return ddc::discrete_space<BSplinesTheta>().full_domain();
    }
};
