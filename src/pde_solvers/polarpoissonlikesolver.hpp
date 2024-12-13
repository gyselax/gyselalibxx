// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/mapping/mapping_tools.hpp>
#include <sll/mapping/metric_tensor.hpp>
#include <sll/math_tools.hpp>
#include <sll/matrix_batch_csr.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>
#include <sll/view.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"

/**
* @brief Define a polar PDE solver for a Poisson-like equation.
 *
 * Solve the following Partial Differential Equation
 *
 * (1) @f$  L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho @f$, in  @f$ \Omega@f$,
 *
 * @f$  \phi = 0 @f$, on  @f$ \partial \Omega@f$,
 * 
 * As finite element basis functions we will use polar b-splines which are divided into two types:
 * 1) Basis splines that can be written as a tensor product of 1d basis splines
 *    ("non-singular bsplines")
 * 2) Basis splines that cover the centre point and are defined as a linear combination
 *    of basis splines of type 1 ("singular bsplines")
 * 
 * (see in Emily Bourne's thesis "Non-Uniform Numerical Schemes for the Modelling of Turbulence
 * in the 5D GYSELA Code". December 2022.)
 *
 * @tparam GridR The radial grid type.
 * @tparam GridR The poloidal grid type.
 * @tparam PolarBSplinesRTheta The type of the 2D polar bsplines (on the coordinate
 * system (r,theta) including bsplines which traverse the O point).
 * @tparam SplineRThetaEvaluatorNullBound The type of the 2D (cross-product) spline evaluator.
 * @tparam IdxRangeFull The full index range of @f$ \phi @f$ including any batch dimensions.
 */
template <
        class GridR,
        class GridTheta,
        class PolarBSplinesRTheta,
        class SplineRThetaEvaluatorNullBound,
        class IdxRangeFull = IdxRange<GridR, GridTheta>>
class PolarSplineFEMPoissonLikeSolver
{
    // TODO: Add a batch loop to operator()
    static_assert(
            std::is_same_v<IdxRangeFull, IdxRange<GridR, GridTheta>>,
            "PolarSplineFEMPoissonLikeSolver is not yet batched");

public:
    /// The radial dimension
    using R = typename GridR::continuous_dimension_type;
    /// The poloidal dimension
    using Theta = typename GridTheta::continuous_dimension_type;

public:
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


public:
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

private:
    using CoordRTheta = Coord<R, Theta>;

    /// The 1D bsplines in the radial direction
    using BSplinesR = typename PolarBSplinesRTheta::BSplinesR_tag;
    /// The 1D bsplines in the poloidal direction
    using BSplinesTheta = typename PolarBSplinesRTheta::BSplinesTheta_tag;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRTheta = Idx<GridR, GridTheta>;

    /// The type of an index range over the polar bsplines
    using IdxRangeBSPolar = IdxRange<PolarBSplinesRTheta>;
    using IdxBSPolar = Idx<PolarBSplinesRTheta>;

    using IdxRangeBSR = IdxRange<BSplinesR>;
    using IdxRangeBSTheta = IdxRange<BSplinesTheta>;
    using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;
    using IdxBSRTheta = Idx<BSplinesR, BSplinesTheta>;

    using IdxRangeBatchedBSRTheta
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<IdxRangeFull>,
                    ddc::detail::TypeSeq<GridR, GridTheta>,
                    ddc::detail::TypeSeq<BSplinesR, BSplinesTheta>>>;

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFull, IdxRange<GridR>, IdxRange<GridTheta>>;

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

    using KnotsR = ddc::NonUniformBsplinesKnots<BSplinesR>;
    using KnotsTheta = ddc::NonUniformBsplinesKnots<BSplinesTheta>;

    using ConstSpline2D = DConstField<IdxRangeBatchedBSRTheta>;
    using SplinePolar = PolarSpline<PolarBSplinesRTheta>;

    using CoordFieldMemRTheta = FieldMem<CoordRTheta, IdxRangeRTheta>;
    using CoordFieldRTheta = Field<CoordRTheta, IdxRangeRTheta>;
    using DFieldRTheta = DField<IdxRangeRTheta>;

    /**
     * @brief Object storing a value and a value of the derivative
     * of a 1D function.
     */
    struct EvalDeriv1DType
    {
        double value;
        double derivative;
    };
    /**
     * @brief Object storing a value and a value of the derivatives
     * in each direction of a 2D function.
     */
    struct EvalDeriv2DType
    {
        double value;
        double radial_derivative;
        double poloidal_derivative;
    };

    /**
     * @brief Tag an index of cell.
     */
    using IdxCell = Idx<RCellDim, ThetaCellDim>;

private:
    static constexpr int m_n_gauss_legendre_r = BSplinesR::degree() + 1;
    static constexpr int m_n_gauss_legendre_theta = BSplinesTheta::degree() + 1;
    // The number of cells (in the radial direction) in which both types of basis splines can be found
    static constexpr int m_n_overlap_cells = PolarBSplinesRTheta::continuity + 1;

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

    const int m_nbasis_r;
    const int m_nbasis_theta;

    // Domains
    IdxRangeBSPolar m_idxrange_fem_non_singular;
    IdxRangeBSR m_idxrange_bsplines_r;
    IdxRangeBSTheta m_idxrange_bsplines_theta;

    IdxRangeQuadratureR m_idxrange_quadrature_r;
    IdxRangeQuadratureTheta m_idxrange_quadrature_theta;
    IdxRangeQuadratureRTheta m_idxrange_quadrature_singular;

    // Gauss-Legendre points and weights
    host_t<FieldMem<Coord<R>, IdxRangeQuadratureR>> m_points_r;
    host_t<FieldMem<Coord<Theta>, IdxRangeQuadratureTheta>> m_points_theta;
    host_t<FieldMem<double, IdxRangeQuadratureR>> m_weights_r;
    host_t<FieldMem<double, IdxRangeQuadratureTheta>> m_weights_theta;

    // Basis Spline values and derivatives at Gauss-Legendre points
    host_t<FieldMem<EvalDeriv2DType, IdxRange<PolarBSplinesRTheta, QDimRMesh, QDimThetaMesh>>>
            m_singular_basis_vals_and_derivs;
    host_t<FieldMem<EvalDeriv1DType, IdxRange<RBasisSubset, QDimRMesh>>> r_basis_vals_and_derivs;
    host_t<FieldMem<EvalDeriv1DType, IdxRange<ThetaBasisSubset, QDimThetaMesh>>>
            m_theta_basis_vals_and_derivs;

    host_t<FieldMem<double, IdxRangeQuadratureRTheta>> m_int_volume;

    PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule, Kokkos::HostSpace>
            m_polar_spline_evaluator;
    std::unique_ptr<MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>
            m_gko_matrix;
    mutable host_t<SplinePolar> m_phi_spline_coef;
    Kokkos::View<double**, Kokkos::LayoutRight> m_x_init;

    const int m_batch_idx {0}; // TODO: Remove when batching is supported

public:
    /**
     * @brief Instantiate a polar Poisson-like solver using FEM with B-splines.
     *
     * The equation we are studying:
     *
     * (1) @f$  L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho @f$, in  @f$ \Omega@f$,
     *
     *  @f$  \phi = 0 @f$, on  @f$ \partial \Omega@f$.
     *
     * @param[in] coeff_alpha
     *      The spline representation of the @f$ \alpha @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] coeff_beta
     *      The spline representation of the  @f$ \beta @f$ function in the
     *      definition of the Poisson-like equation.
     * @param[in] mapping
     *      The mapping from the logical index range to the physical index range where
     *      the equation is defined.
     * @param[in] spline_evaluator
     *      An evaluator for evaluating 2D splines on (r, theta)
     *
     * @tparam Mapping A class describing a mapping from curvilinear coordinates to cartesian coordinates.
     */
    template <class Mapping>
    PolarSplineFEMPoissonLikeSolver(
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            Mapping const& mapping,
            SplineRThetaEvaluatorNullBound const& spline_evaluator)
        : m_nbasis_r(ddc::discrete_space<BSplinesR>().nbasis() - m_n_overlap_cells - 1)
        , m_nbasis_theta(ddc::discrete_space<BSplinesTheta>().nbasis())
        , m_idxrange_fem_non_singular(
                  ddc::discrete_space<PolarBSplinesRTheta>().tensor_bspline_idx_range().remove_last(
                          IdxStep<PolarBSplinesRTheta> {m_nbasis_theta}))
        , m_idxrange_bsplines_r(ddc::discrete_space<BSplinesR>().full_domain().remove_first(
                  IdxStep<BSplinesR> {m_n_overlap_cells}))
        , m_idxrange_bsplines_theta(ddc::discrete_space<BSplinesTheta>().full_domain().take_first(
                  IdxStep<BSplinesTheta> {m_nbasis_theta}))
        , m_idxrange_quadrature_r(
                  Idx<QDimRMesh>(0),
                  IdxStep<QDimRMesh>(
                          m_n_gauss_legendre_r * ddc::discrete_space<BSplinesR>().ncells()))
        , m_idxrange_quadrature_theta(
                  Idx<QDimThetaMesh>(0),
                  IdxStep<QDimThetaMesh>(
                          m_n_gauss_legendre_theta * ddc::discrete_space<BSplinesTheta>().ncells()))
        , m_idxrange_quadrature_singular(
                  m_idxrange_quadrature_r.take_first(
                          IdxStep<QDimRMesh> {m_n_overlap_cells * m_n_gauss_legendre_r}),
                  m_idxrange_quadrature_theta)
        , m_points_r(m_idxrange_quadrature_r)
        , m_points_theta(m_idxrange_quadrature_theta)
        , m_weights_r(m_idxrange_quadrature_r)
        , m_weights_theta(m_idxrange_quadrature_theta)
        , m_singular_basis_vals_and_derivs(IdxRange<PolarBSplinesRTheta, QDimRMesh, QDimThetaMesh>(
                  PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>(),
                  ddc::select<QDimRMesh>(m_idxrange_quadrature_singular),
                  ddc::select<QDimThetaMesh>(m_idxrange_quadrature_singular)))
        , r_basis_vals_and_derivs(
                  IdxRange<RBasisSubset, QDimRMesh>(m_non_zero_bases_r, m_idxrange_quadrature_r))
        , m_theta_basis_vals_and_derivs(
                  IdxRange<
                          ThetaBasisSubset,
                          QDimThetaMesh>(m_non_zero_bases_theta, m_idxrange_quadrature_theta))
        , m_int_volume(
                  IdxRangeQuadratureRTheta(m_idxrange_quadrature_r, m_idxrange_quadrature_theta))
        , m_polar_spline_evaluator(ddc::NullExtrapolationRule())
        , m_phi_spline_coef(
                  PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>(),
                  IdxRangeBSRTheta(
                          m_idxrange_bsplines_r,
                          ddc::discrete_space<BSplinesTheta>().full_domain()))
        , m_x_init(
                  "x_init",
                  1,
                  ddc::discrete_space<PolarBSplinesRTheta>().nbasis()
                          - ddc::discrete_space<BSplinesTheta>().nbasis())
    {
        static_assert(has_2d_jacobian_v<Mapping, CoordRTheta>);
        //initialize x_init
        Kokkos::deep_copy(m_x_init, 0);
        // Get break points
        IdxRange<KnotsR> idxrange_r_edges = ddc::discrete_space<BSplinesR>().break_point_domain();
        IdxRange<KnotsTheta> idxrange_theta_edges
                = ddc::discrete_space<BSplinesTheta>().break_point_domain();
        host_t<FieldMem<Coord<R>, IdxRange<KnotsR>>> breaks_r(idxrange_r_edges);
        host_t<FieldMem<Coord<Theta>, IdxRange<KnotsTheta>>> breaks_theta(idxrange_theta_edges);

        ddc::for_each(idxrange_r_edges, [&](Idx<KnotsR> i) { breaks_r(i) = ddc::coordinate(i); });
        ddc::for_each(idxrange_theta_edges, [&](Idx<KnotsTheta> i) {
            breaks_theta(i) = ddc::coordinate(i);
        });

        // Define quadrature points and weights
        GaussLegendre<R> gl_coeffs_r(m_n_gauss_legendre_r);
        GaussLegendre<Theta> gl_coeffs_theta(m_n_gauss_legendre_theta);
        gl_coeffs_r.compute_points_and_weights_on_mesh(
                get_field(m_points_r),
                get_field(m_weights_r),
                get_const_field(breaks_r));
        gl_coeffs_theta.compute_points_and_weights_on_mesh(
                get_field(m_points_theta),
                get_field(m_weights_theta),
                get_const_field(breaks_theta));

        std::vector<double> vect_points_r(m_points_r.size());
        for (IdxQuadratureR i : m_idxrange_quadrature_r) {
            vect_points_r[i - m_idxrange_quadrature_r.front()] = m_points_r(i);
        }
        std::vector<double> vect_points_theta(m_points_theta.size());
        for (IdxQuadratureTheta i : m_idxrange_quadrature_theta) {
            vect_points_theta[i - m_idxrange_quadrature_theta.front()] = m_points_theta(i);
        }

        // Create quadrature index range
        ddc::init_discrete_space<QDimRMesh>(vect_points_r);
        ddc::init_discrete_space<QDimThetaMesh>(vect_points_theta);

        // Find value and derivative of 1D bsplines in radial direction
        ddc::for_each(m_idxrange_quadrature_r, [&](IdxQuadratureR const idx_r) {
            std::array<double, 2 * m_n_non_zero_bases_r> data;
            DSpan2D vals(data.data(), m_n_non_zero_bases_r, 2);
            ddc::discrete_space<BSplinesR>()
                    .eval_basis_and_n_derivs(vals, ddc::coordinate(idx_r), 1);
            for (auto ib : m_non_zero_bases_r) {
                const int ib_idx = ib - m_non_zero_bases_r.front();
                r_basis_vals_and_derivs(ib, idx_r).value = vals(ib_idx, 0);
                r_basis_vals_and_derivs(ib, idx_r).derivative = vals(ib_idx, 1);
            }
        });

        // Find value and derivative of 1D bsplines in poloidal direction
        ddc::for_each(m_idxrange_quadrature_theta, [&](IdxQuadratureTheta const idx_theta) {
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

        IdxRangeBSPolar idxrange_singular
                = PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>();

        // Find value and derivative of 2D bsplines covering the singular point
        ddc::for_each(m_idxrange_quadrature_singular, [&](IdxQuadratureRTheta const irp) {
            std::array<double, PolarBSplinesRTheta::n_singular_basis()> singular_data;
            std::array<double, m_n_non_zero_bases_r * m_n_non_zero_bases_theta> data;
            // Values of the polar basis splines around the singular point
            // at a given coordinate
            DSpan1D singular_vals(singular_data.data(), PolarBSplinesRTheta::n_singular_basis());
            // Values of the polar basis splines, that do not cover the singular point,
            // at a given coordinate
            DSpan2D vals(data.data(), m_n_non_zero_bases_r, m_n_non_zero_bases_theta);
            IdxQuadratureR idx_r = ddc::select<QDimRMesh>(irp);
            IdxQuadratureTheta idx_theta = ddc::select<QDimThetaMesh>(irp);

            const CoordRTheta coord(ddc::coordinate(irp));

            // Calculate the value
            ddc::discrete_space<PolarBSplinesRTheta>().eval_basis(singular_vals, vals, coord);
            for (IdxBSPolar ib : idxrange_singular) {
                m_singular_basis_vals_and_derivs(ib, idx_r, idx_theta).value
                        = singular_vals[ib - idxrange_singular.front()];
            }

            // Calculate the radial derivative
            ddc::discrete_space<PolarBSplinesRTheta>().eval_deriv_r(singular_vals, vals, coord);
            for (IdxBSPolar ib : idxrange_singular) {
                m_singular_basis_vals_and_derivs(ib, idx_r, idx_theta).radial_derivative
                        = singular_vals[ib - idxrange_singular.front()];
            }

            // Calculate the poloidal derivative
            ddc::discrete_space<PolarBSplinesRTheta>().eval_deriv_theta(singular_vals, vals, coord);
            for (IdxBSPolar ib : idxrange_singular) {
                m_singular_basis_vals_and_derivs(ib, idx_r, idx_theta).poloidal_derivative
                        = singular_vals[ib - idxrange_singular.front()];
            }
        });

        // Find the integral volume associated with each point used in the quadrature scheme
        IdxRangeQuadratureRTheta
                all_quad_points(m_idxrange_quadrature_r, m_idxrange_quadrature_theta);
        ddc::for_each(all_quad_points, [&](IdxQuadratureRTheta const irp) {
            IdxQuadratureR const idx_r = ddc::select<QDimRMesh>(irp);
            IdxQuadratureTheta const idx_theta = ddc::select<QDimThetaMesh>(irp);
            CoordRTheta coord(ddc::coordinate(irp));
            m_int_volume(idx_r, idx_theta) = abs(mapping.jacobian(coord)) * m_weights_r(idx_r)
                                             * m_weights_theta(idx_theta);
        });

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
                = m_nbasis_theta * min(int(1 + 2 * BSplinesTheta::degree()), m_nbasis_theta);
        const int n_stencil_r = m_nbasis_r * (1 + 2 * BSplinesR::degree())
                                - (1 + BSplinesR::degree()) * BSplinesR::degree();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // non-central splines. These have a tensor product structure
        const int n_elements_stencil = n_stencil_r * n_stencil_theta;

        const int batch_size = 1;
        // Matrix size is equal to the number Polar bspline
        const int matrix_size
                = ddc::discrete_space<PolarBSplinesRTheta>().nbasis() - m_nbasis_theta;
        const int n_matrix_elements = n_elements_singular + n_elements_overlap + n_elements_stencil;

        //CSR data storage
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>
                values_csr_host("values_csr", batch_size, n_matrix_elements);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace>
                col_idx_csr_host("idx_csr", n_matrix_elements);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
                nnz_per_row_csr("nnz_per_row_csr", matrix_size + 1);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace>
                nnz_per_row_csr_host("nnz_per_row_csr", matrix_size + 1);

        init_nnz_per_line(nnz_per_row_csr);
        Kokkos::deep_copy(nnz_per_row_csr_host, nnz_per_row_csr);

        Kokkos::Profiling::pushRegion("PolarPoissonFillFemMatrix");
        // Calculate the matrix elements corresponding to the bsplines which cover the singular point
        ddc::for_each(idxrange_singular, [&](IdxBSPolar const idx_test) {
            ddc::for_each(idxrange_singular, [&](IdxBSPolar const idx_trial) {
                // Calculate the weak integral
                double const element = ddc::transform_reduce(
                        m_idxrange_quadrature_singular,
                        0.0,
                        ddc::reducer::sum<double>(),
                        [&](IdxQuadratureRTheta const idx_quad) {
                            IdxQuadratureR const idx_r = ddc::select<QDimRMesh>(idx_quad);
                            IdxQuadratureTheta const idx_theta
                                    = ddc::select<QDimThetaMesh>(idx_quad);
                            return weak_integral_element(
                                    idx_r,
                                    idx_theta,
                                    m_singular_basis_vals_and_derivs(idx_test, idx_r, idx_theta),
                                    m_singular_basis_vals_and_derivs(idx_trial, idx_r, idx_theta),
                                    coeff_alpha,
                                    coeff_beta,
                                    spline_evaluator,
                                    mapping);
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

        // Create index ranges associated with the 2D splines
        IdxRangeBSR central_radial_bspline_idx_range(
                m_idxrange_bsplines_r.take_first(IdxStep<BSplinesR> {BSplinesR::degree()}));

        IdxRangeBSRTheta idxrange_non_singular_near_centre(
                central_radial_bspline_idx_range,
                m_idxrange_bsplines_theta);

        // Calculate the matrix elements where bspline products overlap the bsplines which cover the singular point
        ddc::for_each(idxrange_singular, [&](IdxBSPolar const idx_test) {
            ddc::for_each(idxrange_non_singular_near_centre, [&](IdxBSRTheta const idx_trial) {
                const IdxBSPolar idx_trial_polar(
                        PolarBSplinesRTheta::template get_polar_index<PolarBSplinesRTheta>(
                                idx_trial));
                const Idx<BSplinesR> idx_trial_r(ddc::select<BSplinesR>(idx_trial));
                const Idx<BSplinesTheta> idx_trial_theta(ddc::select<BSplinesTheta>(idx_trial));

                // Find the index range covering the cells where both the test and trial functions are non-zero
                const Idx<RCellDim> first_overlap_element_r(
                        idx_trial_r.uid() < BSplinesR::degree()
                                ? 0
                                : idx_trial_r.uid() - BSplinesR::degree());
                const Idx<ThetaCellDim> first_overlap_element_theta(
                        theta_mod(idx_trial_theta.uid() - BSplinesTheta::degree()));

                const IdxStep<RCellDim> n_overlap_r(
                        m_n_overlap_cells - first_overlap_element_r.uid());
                const IdxStep<ThetaCellDim> n_overlap_theta(BSplinesTheta::degree() + 1);

                const IdxRange<RCellDim> r_cells(first_overlap_element_r, n_overlap_r);
                const IdxRange<ThetaCellDim>
                        theta_cells(first_overlap_element_theta, n_overlap_theta);
                const IdxRange<RCellDim, ThetaCellDim> non_zero_cells(r_cells, theta_cells);

                if (n_overlap_r > 0) {
                    double element = 0.0;

                    ddc::for_each(non_zero_cells, [&](IdxCell const cell_idx) {
                        const int cell_idx_r(ddc::select<RCellDim>(cell_idx).uid());
                        const int cell_idx_theta(
                                theta_mod(ddc::select<ThetaCellDim>(cell_idx).uid()));

                        const IdxRangeQuadratureRTheta cell_quad_points(
                                get_quadrature_points_in_cell(cell_idx_r, cell_idx_theta));
                        // Find the column where the non-zero data is stored
                        Idx<RBasisSubset> ib_trial_r(idx_trial_r.uid() - cell_idx_r);
                        Idx<ThetaBasisSubset> ib_trial_theta(
                                theta_mod(idx_trial_theta.uid() - cell_idx_theta));
                        // Calculate the weak integral
                        element += ddc::transform_reduce(
                                cell_quad_points,
                                0.0,
                                ddc::reducer::sum<double>(),
                                [&](IdxQuadratureRTheta const idx_quad) {
                                    IdxQuadratureR const idx_r = ddc::select<QDimRMesh>(idx_quad);
                                    IdxQuadratureTheta const idx_theta
                                            = ddc::select<QDimThetaMesh>(idx_quad);
                                    return weak_integral_element<Mapping>(
                                            idx_r,
                                            idx_theta,
                                            m_singular_basis_vals_and_derivs(
                                                    idx_test,
                                                    idx_r,
                                                    idx_theta),
                                            r_basis_vals_and_derivs(ib_trial_r, idx_r),
                                            m_theta_basis_vals_and_derivs(
                                                    ib_trial_theta,
                                                    idx_theta),
                                            coeff_alpha,
                                            coeff_beta,
                                            spline_evaluator,
                                            mapping);
                                });
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
                }
            });
        });

        // Calculate the matrix elements following a stencil
        ddc::for_each(m_idxrange_fem_non_singular, [&](IdxBSPolar const idx_test_polar) {
            const IdxBSRTheta idx_test(PolarBSplinesRTheta::get_2d_index(idx_test_polar));
            const std::size_t idx_test_r(ddc::select<BSplinesR>(idx_test).uid());
            const std::size_t idx_test_theta(ddc::select<BSplinesTheta>(idx_test).uid());

            // Calculate the index of the elements that are already filled
            IdxRangeBSTheta remaining_theta(
                    Idx<BSplinesTheta> {idx_test_theta},
                    IdxStep<BSplinesTheta> {BSplinesTheta::degree() + 1});
            ddc::for_each(remaining_theta, [&](Idx<BSplinesTheta> const idx_trial_theta) {
                IdxBSRTheta idx_trial(Idx<BSplinesR>(idx_test_r), idx_trial_theta);
                IdxBSPolar idx_trial_polar(
                        PolarBSplinesRTheta::template get_polar_index<PolarBSplinesRTheta>(
                                IdxBSRTheta(idx_test_r, theta_mod(idx_trial_theta.uid()))));
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
            IdxRangeBSR remaining_r(
                    ddc::select<BSplinesR>(idx_test) + 1,
                    IdxStep<BSplinesR> {
                            min(BSplinesR::degree(),
                                ddc::discrete_space<BSplinesR>().nbasis() - 2 - idx_test_r)});
            IdxRangeBSTheta relevant_theta(
                    Idx<BSplinesTheta> {
                            idx_test_theta + ddc::discrete_space<BSplinesTheta>().nbasis()
                            - BSplinesTheta::degree()},
                    IdxStep<BSplinesTheta> {2 * BSplinesTheta::degree() + 1});

            IdxRangeBSRTheta trial_idx_range(remaining_r, relevant_theta);

            ddc::for_each(trial_idx_range, [&](IdxBSRTheta const idx_trial) {
                const int idx_trial_r(ddc::select<BSplinesR>(idx_trial).uid());
                const int idx_trial_theta(ddc::select<BSplinesTheta>(idx_trial).uid());
                IdxBSPolar idx_trial_polar(
                        PolarBSplinesRTheta::template get_polar_index<PolarBSplinesRTheta>(
                                IdxBSRTheta(idx_trial_r, theta_mod(idx_trial_theta))));
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
        assert(nnz_per_row_csr_host(matrix_size) == n_matrix_elements);
        m_gko_matrix = std::make_unique<MatrixBatchCsr<
                Kokkos::DefaultExecutionSpace,
                MatrixBatchCsrSolver::CG>>(1, matrix_size, n_matrix_elements);

        auto [values, col_idx, nnz_per_row] = m_gko_matrix->get_batch_csr();
        Kokkos::deep_copy(values, values_csr_host);
        Kokkos::deep_copy(col_idx, col_idx_csr_host);
        Kokkos::deep_copy(nnz_per_row, nnz_per_row_csr_host);
        m_gko_matrix->setup_solver();
        Kokkos::Profiling::popRegion();
    }

    /**
     * @brief Solve the Poisson-like equation.
     *
     * This operator returns the coefficients associated with the B-Splines
     * of the solution @f$\phi@f$.
     *
     * @param[in] rhs
     *      The rhs @f$ \rho@f$ of the Poisson-like equation.
     *      The type is templated but we can use the PoissonLikeRHSFunction
     *      class.
     * @param[inout] spline
     *      The spline representation of the solution @f$\phi@f$, also used as initial data for the iterative solver.
     */
    template <class RHSFunction>
    void operator()(RHSFunction const& rhs, host_t<SplinePolar>& spline) const
    {
        Kokkos::Profiling::pushRegion("PolarPoissonRHS");

        static_assert(
                std::is_invocable_r_v<double, RHSFunction, CoordRTheta>,
                "RHSFunction must have an operator() which takes a coordinate and returns a "
                "double");
        const int b_size = ddc::discrete_space<PolarBSplinesRTheta>().nbasis()
                           - ddc::discrete_space<BSplinesTheta>().nbasis();
        const int batch_size = 1;
        // Create b for rhs
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>
                b_host("b_host", batch_size, b_size);
        //Create an initial guess
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>
                x_init_host("x_init_host", batch_size, b_size);
        // Fill b
        ddc::for_each(
                PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>(),
                [&](IdxBSPolar const idx) {
                    const int bspl_idx = idx
                                         - PolarBSplinesRTheta::template singular_idx_range<
                                                   PolarBSplinesRTheta>()
                                                   .front();
                    b_host(0, bspl_idx) = ddc::transform_reduce(
                            m_idxrange_quadrature_singular,
                            0.0,
                            ddc::reducer::sum<double>(),
                            [&](IdxQuadratureRTheta const idx_quad) {
                                IdxQuadratureR const idx_r = ddc::select<QDimRMesh>(idx_quad);
                                IdxQuadratureTheta const idx_theta
                                        = ddc::select<QDimThetaMesh>(idx_quad);
                                CoordRTheta coord(ddc::coordinate(idx_quad));
                                return rhs(coord)
                                       * m_singular_basis_vals_and_derivs(idx, idx_r, idx_theta)
                                                 .value
                                       * m_int_volume(idx_r, idx_theta);
                            });
                });
        const std::size_t ncells_r = ddc::discrete_space<BSplinesR>().ncells();
        ddc::for_each(m_idxrange_fem_non_singular, [&](IdxBSPolar const idx) {
            const IdxBSRTheta idx_2d(PolarBSplinesRTheta::get_2d_index(idx));
            const std::size_t idx_r(ddc::select<BSplinesR>(idx_2d).uid());
            const std::size_t idx_theta(ddc::select<BSplinesTheta>(idx_2d).uid());

            // Find the cells on which the bspline is non-zero
            int first_cell_r(idx_r - BSplinesR::degree());
            int first_cell_theta(idx_theta - BSplinesTheta::degree());
            std::size_t last_cell_r(idx_r + 1);
            if (first_cell_r < 0)
                first_cell_r = 0;
            if (last_cell_r > ncells_r)
                last_cell_r = ncells_r;
            IdxStep<RCellDim> const r_length(last_cell_r - first_cell_r);
            IdxStep<ThetaCellDim> const theta_length(BSplinesTheta::degree() + 1);


            Idx<RCellDim> const start_r(first_cell_r);
            Idx<ThetaCellDim> const start_theta(theta_mod(first_cell_theta));
            const IdxRange<RCellDim> r_cells(start_r, r_length);
            const IdxRange<ThetaCellDim> theta_cells(start_theta, theta_length);
            const IdxRange<RCellDim, ThetaCellDim> non_zero_cells(r_cells, theta_cells);
            assert(r_length * theta_length > 0);
            double element = 0.0;
            ddc::for_each(non_zero_cells, [&](IdxCell const cell_idx) {
                const int cell_idx_r(ddc::select<RCellDim>(cell_idx).uid());
                const int cell_idx_theta(theta_mod(ddc::select<ThetaCellDim>(cell_idx).uid()));

                const IdxRangeQuadratureRTheta cell_quad_points(
                        get_quadrature_points_in_cell(cell_idx_r, cell_idx_theta));

                // Find the column where the non-zero data is stored
                Idx<RBasisSubset> ib_r(idx_r - cell_idx_r);
                Idx<ThetaBasisSubset> ib_theta(theta_mod(idx_theta - cell_idx_theta));

                // Calculate the weak integral
                element += ddc::transform_reduce(
                        cell_quad_points,
                        0.0,
                        ddc::reducer::sum<double>(),
                        [&](IdxQuadratureRTheta const idx_quad) {
                            IdxQuadratureR const idx_r = ddc::select<QDimRMesh>(idx_quad);
                            IdxQuadratureTheta const idx_theta
                                    = ddc::select<QDimThetaMesh>(idx_quad);
                            CoordRTheta coord(ddc::coordinate(idx_quad));
                            double rb = r_basis_vals_and_derivs(ib_r, idx_r).value;
                            double pb = m_theta_basis_vals_and_derivs(ib_theta, idx_theta).value;
                            return rhs(coord) * rb * pb * m_int_volume(idx_r, idx_theta);
                        });
            });
            const std::size_t singular_index
                    = idx - ddc::discrete_space<PolarBSplinesRTheta>().full_domain().front();
            b_host(0, singular_index) = element;
        });

        Kokkos::View<double**, Kokkos::LayoutRight> b("b", batch_size, b_size);
        Kokkos::deep_copy(b, b_host);
        Kokkos::Profiling::popRegion();

        Kokkos::deep_copy(m_x_init, x_init_host);
        // Solve the matrix equation
        Kokkos::Profiling::pushRegion("PolarPoissonSolve");
        m_gko_matrix->solve(m_x_init, b);
        Kokkos::deep_copy(x_init_host, m_x_init);
        //-----------------
        IdxRangeBSRTheta dirichlet_boundary_idx_range(
                m_idxrange_bsplines_r.take_last(IdxStep<BSplinesR> {1}),
                m_idxrange_bsplines_theta);
        IdxRangeBSTheta idxrange_polar(ddc::discrete_space<BSplinesTheta>().full_domain());


        // Fill the spline
        ddc::for_each(
                PolarBSplinesRTheta::template singular_idx_range<PolarBSplinesRTheta>(),
                [&](IdxBSPolar const idx) {
                    const int bspl_idx = idx
                                         - PolarBSplinesRTheta::template singular_idx_range<
                                                   PolarBSplinesRTheta>()
                                                   .front();
                    spline.singular_spline_coef(idx) = x_init_host(0, bspl_idx);
                });
        ddc::for_each(m_idxrange_fem_non_singular, [&](IdxBSPolar const idx) {
            const IdxBSRTheta idx_2d(PolarBSplinesRTheta::get_2d_index(idx));
            spline.spline_coef(idx_2d) = x_init_host(0, idx.uid());
        });
        ddc::for_each(dirichlet_boundary_idx_range, [&](IdxBSRTheta const idx) {
            spline.spline_coef(idx) = 0.0;
        });

        // Copy the periodic elements
        IdxRangeBSRTheta copy_idx_range(
                m_idxrange_bsplines_r,
                idxrange_polar.remove_first(
                        IdxStep<BSplinesTheta>(ddc::discrete_space<BSplinesTheta>().nbasis())));
        ddc::for_each(copy_idx_range, [&](IdxBSRTheta const idx_2d) {
            spline.spline_coef(ddc::select<BSplinesR>(idx_2d), ddc::select<BSplinesTheta>(idx_2d))
                    = spline.spline_coef(
                            ddc::select<BSplinesR>(idx_2d),
                            ddc::select<BSplinesTheta>(idx_2d)
                                    - ddc::discrete_space<BSplinesTheta>().nbasis());
        });
        Kokkos::Profiling::popRegion();
    }


    /**
     * @brief Solve the Poisson-like equation.
     *
     * This operator uses the other operator () and returns the values on
     * the grid of the solution @f$\phi@f$.
     *
     * @param[in] rhs
     *      The rhs @f$ \rho@f$ of the Poisson-like equation.
     *      The type is templated but we can use the PoissonLikeRHSFunction
     *      class.
     * @param[inout] phi
     *      The values of the solution @f$\phi@f$ on the given coords_eval, also used as initial data for the iterative solver.
     */
    template <class RHSFunction>
    void operator()(RHSFunction const& rhs, host_t<DFieldRTheta> phi) const
    {
        static_assert(
                std::is_invocable_r_v<double, RHSFunction, CoordRTheta>,
                "RHSFunction must have an operator() which takes a coordinate and returns a "
                "double");


        (*this)(rhs, m_phi_spline_coef);
        host_t<CoordFieldMemRTheta> coords_eval_alloc(get_idx_range(phi));
        host_t<CoordFieldRTheta> coords_eval(get_field(coords_eval_alloc));
        ddc::for_each(get_idx_range(phi), [&](IdxRTheta idx) {
            coords_eval(idx) = ddc::coordinate(idx);
        });
        m_polar_spline_evaluator(phi, get_const_field(coords_eval), m_phi_spline_coef);
    }

private:
    static KOKKOS_FUNCTION IdxRangeQuadratureRTheta
    get_quadrature_points_in_cell(int cell_idx_r, int cell_idx_theta)
    {
        const IdxQuadratureR first_quad_point_r(cell_idx_r * m_n_gauss_legendre_r);
        const IdxQuadratureTheta first_quad_point_theta(cell_idx_theta * m_n_gauss_legendre_theta);
        constexpr IdxStepQuadratureR n_GL_r(m_n_gauss_legendre_r);
        constexpr IdxStepQuadratureTheta n_GL_theta(m_n_gauss_legendre_theta);
        const IdxRangeQuadratureR quad_points_r(first_quad_point_r, n_GL_r);
        const IdxRangeQuadratureTheta quad_points_theta(first_quad_point_theta, n_GL_theta);
        return IdxRangeQuadratureRTheta(quad_points_r, quad_points_theta);
    }

    template <class Mapping>
    double weak_integral_element(
            IdxQuadratureR idx_r,
            IdxQuadratureTheta idx_theta,
            EvalDeriv2DType const& test_bspline_val_and_deriv,
            EvalDeriv2DType const& trial_bspline_val_and_deriv,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                idx_r,
                idx_theta,
                test_bspline_val_and_deriv,
                trial_bspline_val_and_deriv,
                test_bspline_val_and_deriv,
                trial_bspline_val_and_deriv,
                coeff_alpha,
                coeff_beta,
                evaluator,
                mapping);
    }

    template <class Mapping>
    double weak_integral_element(
            IdxQuadratureR idx_r,
            IdxQuadratureTheta idx_theta,
            EvalDeriv2DType const& test_bspline_val_and_deriv,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_r,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_theta,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                idx_r,
                idx_theta,
                test_bspline_val_and_deriv,
                trial_bspline_val_and_deriv_r,
                test_bspline_val_and_deriv,
                trial_bspline_val_and_deriv_theta,
                coeff_alpha,
                coeff_beta,
                evaluator,
                mapping);
    }

    template <class Mapping>
    double weak_integral_element(
            IdxQuadratureR idx_r,
            IdxQuadratureTheta idx_theta,
            EvalDeriv1DType const& test_bspline_val_and_deriv_r,
            EvalDeriv2DType const& trial_bspline_val_and_deriv,
            EvalDeriv1DType const& test_bspline_val_and_deriv_theta,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                idx_r,
                idx_theta,
                test_bspline_val_and_deriv_r,
                trial_bspline_val_and_deriv,
                test_bspline_val_and_deriv_theta,
                trial_bspline_val_and_deriv,
                coeff_alpha,
                coeff_beta,
                evaluator,
                mapping);
    }

    template <class Mapping>
    double weak_integral_element(
            IdxQuadratureR idx_r,
            IdxQuadratureTheta idx_theta,
            EvalDeriv1DType const& test_bspline_val_and_deriv_r,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_r,
            EvalDeriv1DType const& test_bspline_val_and_deriv_theta,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_theta,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                idx_r,
                idx_theta,
                test_bspline_val_and_deriv_r,
                trial_bspline_val_and_deriv_r,
                test_bspline_val_and_deriv_theta,
                trial_bspline_val_and_deriv_theta,
                coeff_alpha,
                coeff_beta,
                evaluator,
                mapping);
    }

    inline void get_value_and_gradient(
            double& value,
            std::array<double, 2>& gradient,
            EvalDeriv1DType const& r_basis,
            EvalDeriv1DType const& theta_basis) const
    {
        value = r_basis.value * theta_basis.value;
        gradient = {r_basis.derivative * theta_basis.value, r_basis.value * theta_basis.derivative};
    }

    inline void get_value_and_gradient(
            double& value,
            std::array<double, 2>& gradient,
            EvalDeriv2DType const& basis,
            EvalDeriv2DType const&) const // Last argument is duplicate
    {
        value = basis.value;
        gradient = {basis.radial_derivative, basis.poloidal_derivative};
    }

    /**
     * @brief Computes a quadrature summand corresponding to the 
     *        inner product
     * 
     * The inner product of the test and trial spline is computed using a 
     * quadrature. This function returns one summand of the quadrature for 
     * the quadrature point given by the indices.
     */

    template <class Mapping, class TestValDerivType, class TrialValDerivType>
    double templated_weak_integral_element(
            IdxQuadratureR idx_r,
            IdxQuadratureTheta idx_theta,
            TestValDerivType const& test_bspline_val_and_deriv,
            TrialValDerivType const& trial_bspline_val_and_deriv,
            TestValDerivType const& test_bspline_val_and_deriv_theta,
            TrialValDerivType const& trial_bspline_val_and_deriv_theta,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& spline_evaluator,
            Mapping const& mapping)
    {
        static_assert(
                std::is_same_v<
                        TestValDerivType,
                        EvalDeriv1DType> || std::is_same_v<TestValDerivType, EvalDeriv2DType>);
        static_assert(
                std::is_same_v<
                        TrialValDerivType,
                        EvalDeriv1DType> || std::is_same_v<TrialValDerivType, EvalDeriv2DType>);

        // Calculate coefficients at quadrature point
        CoordRTheta coord(ddc::coordinate(idx_r), ddc::coordinate(idx_theta));
        const double alpha = spline_evaluator(coord, coeff_alpha);
        const double beta = spline_evaluator(coord, coeff_beta);

        // Define the value and gradient of the test and trial basis functions
        double basis_val_test_space;
        double basis_val_trial_space;
        std::array<double, 2> basis_gradient_test_space;
        std::array<double, 2> basis_gradient_trial_space;
        get_value_and_gradient(
                basis_val_test_space,
                basis_gradient_test_space,
                test_bspline_val_and_deriv,
                test_bspline_val_and_deriv_theta);
        get_value_and_gradient(
                basis_val_trial_space,
                basis_gradient_trial_space,
                trial_bspline_val_and_deriv,
                trial_bspline_val_and_deriv_theta);

        MetricTensor<Mapping, CoordRTheta> metric_tensor(mapping);

        // Assemble the weak integral element
        return m_int_volume(idx_r, idx_theta)
               * (alpha
                          * dot_product(
                                  basis_gradient_test_space,
                                  metric_tensor.to_covariant(basis_gradient_trial_space, coord))
                  + beta * basis_val_test_space * basis_val_trial_space);
    }

    /**
     * @brief Computes the matrix element corresponding to two tensor product splines
     *        with index idx_test and idx_trial
     */
    template <class Mapping>
    double get_matrix_stencil_element(
            IdxBSRTheta idx_test,
            IdxBSRTheta idx_trial,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
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
        int idx_test_theta(theta_mod(ddc::select<BSplinesTheta>(idx_test).uid()));
        int idx_trial_theta(theta_mod(ddc::select<BSplinesTheta>(idx_trial).uid()));

        const std::size_t ncells_r = ddc::discrete_space<BSplinesR>().ncells();

        // 0<= r_offset <= degree_r
        // -degree_theta <= theta_offset <= degree_theta
        const int r_offset = idx_trial_r - idx_test_r;
        int theta_offset = theta_mod(idx_trial_theta - idx_test_theta);
        if (theta_offset >= int(m_nbasis_theta - BSplinesTheta::degree())) {
            theta_offset -= m_nbasis_theta;
        }
        assert(r_offset >= 0);
        assert(r_offset <= int(BSplinesR::degree()));
        assert(theta_offset >= -int(BSplinesTheta::degree()));
        assert(theta_offset <= int(BSplinesTheta::degree()));

        // Find the index range covering the cells where both the test and trial functions are non-zero
        int n_overlap_stencil_r(BSplinesR::degree() + 1 - r_offset);
        int first_overlap_r(idx_trial_r - BSplinesR::degree());

        int first_overlap_theta;
        int n_overlap_stencil_theta;
        if (theta_offset > 0) {
            n_overlap_stencil_theta = BSplinesTheta::degree() + 1 - theta_offset;
            first_overlap_theta = theta_mod(idx_trial_theta - BSplinesTheta::degree());
        } else {
            n_overlap_stencil_theta = BSplinesTheta::degree() + 1 + theta_offset;
            first_overlap_theta = theta_mod(idx_test_theta - BSplinesTheta::degree());
        }

        if (first_overlap_r < 0) {
            const int n_compact = first_overlap_r;
            first_overlap_r = 0;
            n_overlap_stencil_r += n_compact;
        }

        const int n_to_edge_r(ncells_r - first_overlap_r);

        const IdxStep<RCellDim> n_overlap_r(min(n_overlap_stencil_r, n_to_edge_r));
        const IdxStep<ThetaCellDim> n_overlap_theta(n_overlap_stencil_theta);

        const Idx<RCellDim> first_overlap_element_r(first_overlap_r);
        const Idx<ThetaCellDim> first_overlap_element_theta(first_overlap_theta);

        const IdxRange<RCellDim> r_cells(first_overlap_element_r, n_overlap_r);
        const IdxRange<ThetaCellDim> theta_cells(first_overlap_element_theta, n_overlap_theta);
        const IdxRange<RCellDim, ThetaCellDim> non_zero_cells(r_cells, theta_cells);

        assert(n_overlap_r * n_overlap_theta > 0);
        return ddc::transform_reduce(
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
                    return ddc::transform_reduce(
                            cell_quad_points,
                            0.0,
                            ddc::reducer::sum<double>(),
                            [&](IdxQuadratureRTheta const idx_quad) {
                                IdxQuadratureR const idx_r = ddc::select<QDimRMesh>(idx_quad);
                                IdxQuadratureTheta const idx_theta
                                        = ddc::select<QDimThetaMesh>(idx_quad);
                                return weak_integral_element(
                                        idx_r,
                                        idx_theta,
                                        r_basis_vals_and_derivs(ib_test_r, idx_r),
                                        r_basis_vals_and_derivs(ib_trial_r, idx_r),
                                        m_theta_basis_vals_and_derivs(ib_test_theta, idx_theta),
                                        m_theta_basis_vals_and_derivs(ib_trial_theta, idx_theta),
                                        coeff_alpha,
                                        coeff_beta,
                                        evaluator,
                                        mapping);
                            });
                });
    }

    static KOKKOS_FUNCTION int theta_mod(int idx_theta)
    {
        int ncells_theta = ddc::discrete_space<BSplinesTheta>().ncells();
        while (idx_theta < 0)
            idx_theta += ncells_theta;
        while (idx_theta >= ncells_theta)
            idx_theta -= ncells_theta;
        return idx_theta;
    }

public:
    /**
    * @brief Fills the nnz data structure by computing the number of non-zero per line.
    * This number is linked to the weak formulation and depends on @f$(r,\theta)@f$ splines.
    * After this function the array will contain:
    * nnz[0] = 0.
    * nnz[1] = 0.
    * nnz[2] = number of non-zero elements in line 0.
    * nnz[3] = number of non-zero elements in lines 0-1.
    * ...
    * nnz[matrix_size] = number of non-zero elements in lines 0-(matrix_size-1).
    *
    * @param[out]  nnz A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
    */
    void init_nnz_per_line(Kokkos::View<int*, Kokkos::LayoutRight> nnz) const
    {
        Kokkos::Profiling::pushRegion("PolarPoissonInitNnz");
        size_t const mat_size = nnz.extent(0) - 1;
        size_t constexpr n_singular_basis = PolarBSplinesRTheta::n_singular_basis();
        size_t constexpr degree = BSplinesR::degree();
        size_t constexpr radial_overlap = 2 * degree + 1;
        size_t const nbasis_theta_proxy = m_nbasis_theta;

        // overlap between singular domain splines and radial splines
        Kokkos::parallel_for(
                "overlap singular radial",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(1, n_singular_basis + 1),
                KOKKOS_LAMBDA(const int k) {
                    nnz(k + 1) = n_singular_basis + degree * nbasis_theta_proxy;
                });

        // going from the internal boundary the overlapping possiblities between two radial splines increase
        Kokkos::parallel_for(
                "inner overlap",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(1, degree + 2),
                KOKKOS_LAMBDA(const int i) {
                    for (size_t k = n_singular_basis + (i - 1) * nbasis_theta_proxy;
                         k < n_singular_basis + i * nbasis_theta_proxy;
                         k++) {
                        nnz(k + 2) = n_singular_basis + (degree + i) * radial_overlap;
                    }
                });

        // Stencil with maximum possible overlap from two sides for radial spline
        Kokkos::parallel_for(
                "Inner Stencil",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                        n_singular_basis + degree * nbasis_theta_proxy,
                        mat_size - degree * nbasis_theta_proxy),
                KOKKOS_LAMBDA(const int k) { nnz(k + 2) = radial_overlap * radial_overlap; });

        // Approaching the external boundary the overlapping possiblities between two radial splines decrease
        Kokkos::parallel_for(
                "outer overlap",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(1, degree + 1),
                KOKKOS_LAMBDA(const int i) {
                    for (size_t k = mat_size - i * nbasis_theta_proxy;
                         k < mat_size - (i - 1) * nbasis_theta_proxy;
                         k++) {
                        nnz(k + 2) = (degree + i) * radial_overlap;
                    }
                });

        // sum non-zero elements count
        Kokkos::parallel_for(
                "Sum over lines",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, 1),
                KOKKOS_LAMBDA(const int idx) {
                    for (size_t k = 1; k < mat_size; k++) {
                        nnz(k + 1) += nnz(k);
                    }
                    nnz(0) = 0;
                    nnz(1) = 0;
                });
        Kokkos::Profiling::popRegion();
    }
};
