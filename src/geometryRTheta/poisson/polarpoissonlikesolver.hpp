// SPDX-License-Identifier: MIT
#pragma once
#include <iomanip>

#include <ddc/ddc.hpp>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/mapping/metric_tensor.hpp>
#include <sll/math_tools.hpp>
#include <sll/matrix_batch_csr.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>
#include <sll/view.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"

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
 */
class PolarSplineFEMPoissonLikeSolver
{
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

    using BSplinesR_Polar = PolarBSplinesRTheta::BSplinesR_tag;
    using BSplinesTheta_Polar = PolarBSplinesRTheta::BSplinesTheta_tag;

    using IdxBSpline2D_Polar = Idx<BSplinesR_Polar, BSplinesTheta_Polar>;

    using IdxRangeBSR_Polar = IdxRange<BSplinesR_Polar>;
    using IdxRangeBSTheta_Polar = IdxRange<BSplinesTheta_Polar>;

    using KnotsR = ddc::NonUniformBsplinesKnots<BSplinesR_Polar>;
    using KnotsTheta = ddc::NonUniformBsplinesKnots<BSplinesTheta_Polar>;
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

    /**
     * @brief Tag type of matrix element.
     */

private:
    static constexpr int n_gauss_legendre_r = BSplinesR_Polar::degree() + 1;
    static constexpr int n_gauss_legendre_theta = BSplinesTheta_Polar::degree() + 1;
    // The number of cells (in the radial direction) in which both types of basis splines can be found
    static constexpr int n_overlap_cells = PolarBSplinesRTheta::continuity + 1;

    // Number of cells over which a radial B-splines has its support
    // This is the case for b-splines which are not affected by the higher knot multiplicity at the boundary.
    static constexpr IdxStep<RBasisSubset> n_non_zero_bases_r
            = IdxStep<RBasisSubset>(BSplinesR_Polar::degree() + 1);

    // Number of cells over which a poloidal B-splines has its support
    static constexpr IdxStep<ThetaBasisSubset> n_non_zero_bases_theta
            = IdxStep<ThetaBasisSubset>(BSplinesTheta_Polar::degree() + 1);

    static constexpr IdxRange<RBasisSubset> non_zero_bases_r
            = IdxRange<RBasisSubset>(Idx<RBasisSubset> {0}, n_non_zero_bases_r);
    static constexpr IdxRange<ThetaBasisSubset> non_zero_bases_theta
            = IdxRange<ThetaBasisSubset>(Idx<ThetaBasisSubset> {0}, n_non_zero_bases_theta);

    const int nbasis_r;
    const int nbasis_theta;

    // Domains
    IdxRangeBSPolar fem_non_singular_idx_range;
    IdxRangeBSR_Polar radial_bsplines;
    IdxRangeBSTheta_Polar polar_bsplines;

    IdxRangeQuadratureR quadrature_idx_range_r;
    IdxRangeQuadratureTheta quadrature_idx_range_theta;
    IdxRangeQuadratureRTheta quadrature_idx_range_singular;

    // Gauss-Legendre points and weights
    host_t<FieldMem<Coord<R>, IdxRangeQuadratureR>> points_r;
    host_t<FieldMem<Coord<Theta>, IdxRangeQuadratureTheta>> points_theta;
    host_t<FieldMem<double, IdxRangeQuadratureR>> weights_r;
    host_t<FieldMem<double, IdxRangeQuadratureTheta>> weights_theta;

    // Basis Spline values and derivatives at Gauss-Legendre points
    host_t<FieldMem<EvalDeriv2DType, IdxRange<PolarBSplinesRTheta, QDimRMesh, QDimThetaMesh>>>
            singular_basis_vals_and_derivs;
    host_t<FieldMem<EvalDeriv1DType, IdxRange<RBasisSubset, QDimRMesh>>> r_basis_vals_and_derivs;
    host_t<FieldMem<EvalDeriv1DType, IdxRange<ThetaBasisSubset, QDimThetaMesh>>>
            theta_basis_vals_and_derivs;

    host_t<FieldMem<double, IdxRangeQuadratureRTheta>> int_volume;

    PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule> m_polar_spline_evaluator;
    std::unique_ptr<MatrixBatchCsr<Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG>>
            m_gko_matrix;
    const int batch_idx {0};

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
     *
     * @tparam Mapping A Curvilinear2DToCartesian class.
     */
    template <class Mapping>
    PolarSplineFEMPoissonLikeSolver(
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            Mapping const& mapping)
        : nbasis_r(ddc::discrete_space<BSplinesR_Polar>().nbasis() - n_overlap_cells - 1)
        , nbasis_theta(ddc::discrete_space<BSplinesTheta_Polar>().nbasis())
        , fem_non_singular_idx_range(
                  ddc::discrete_space<PolarBSplinesRTheta>().tensor_bspline_idx_range().remove_last(
                          IdxStep<PolarBSplinesRTheta> {nbasis_theta}))
        , radial_bsplines(ddc::discrete_space<BSplinesR_Polar>().full_domain().remove_first(
                  IdxStep<BSplinesR_Polar> {n_overlap_cells}))
        , polar_bsplines(ddc::discrete_space<BSplinesTheta_Polar>().full_domain().take_first(
                  IdxStep<BSplinesTheta_Polar> {nbasis_theta}))
        , quadrature_idx_range_r(
                  Idx<QDimRMesh>(0),
                  IdxStep<QDimRMesh>(
                          n_gauss_legendre_r * ddc::discrete_space<BSplinesR_Polar>().ncells()))
        , quadrature_idx_range_theta(
                  Idx<QDimThetaMesh>(0),
                  IdxStep<QDimThetaMesh>(
                          n_gauss_legendre_theta
                          * ddc::discrete_space<BSplinesTheta_Polar>().ncells()))
        , quadrature_idx_range_singular(
                  quadrature_idx_range_r.take_first(
                          IdxStep<QDimRMesh> {n_overlap_cells * n_gauss_legendre_r}),
                  quadrature_idx_range_theta)
        , points_r(quadrature_idx_range_r)
        , points_theta(quadrature_idx_range_theta)
        , weights_r(quadrature_idx_range_r)
        , weights_theta(quadrature_idx_range_theta)
        , singular_basis_vals_and_derivs(IdxRange<PolarBSplinesRTheta, QDimRMesh, QDimThetaMesh>(
                  PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>(),
                  ddc::select<QDimRMesh>(quadrature_idx_range_singular),
                  ddc::select<QDimThetaMesh>(quadrature_idx_range_singular)))
        , r_basis_vals_and_derivs(
                  IdxRange<RBasisSubset, QDimRMesh>(non_zero_bases_r, quadrature_idx_range_r))
        , theta_basis_vals_and_derivs(
                  IdxRange<
                          ThetaBasisSubset,
                          QDimThetaMesh>(non_zero_bases_theta, quadrature_idx_range_theta))
        , int_volume(IdxRangeQuadratureRTheta(quadrature_idx_range_r, quadrature_idx_range_theta))
        , m_polar_spline_evaluator(ddc::NullExtrapolationRule())
    {
        // Get break points
        IdxRange<KnotsR> idx_range_r_edges
                = ddc::discrete_space<BSplinesR_Polar>().break_point_domain();
        IdxRange<KnotsTheta> idx_range_theta_edges
                = ddc::discrete_space<BSplinesTheta_Polar>().break_point_domain();
        host_t<FieldMem<Coord<R>, IdxRange<KnotsR>>> breaks_r(idx_range_r_edges);
        host_t<FieldMem<Coord<Theta>, IdxRange<KnotsTheta>>> breaks_theta(idx_range_theta_edges);

        ddc::for_each(idx_range_r_edges, [&](Idx<KnotsR> i) { breaks_r(i) = ddc::coordinate(i); });
        ddc::for_each(idx_range_theta_edges, [&](Idx<KnotsTheta> i) {
            breaks_theta(i) = ddc::coordinate(i);
        });

        // Define quadrature points and weights
        GaussLegendre<R> gl_coeffs_r(n_gauss_legendre_r);
        GaussLegendre<Theta> gl_coeffs_theta(n_gauss_legendre_theta);
        gl_coeffs_r.compute_points_and_weights_on_mesh(
                get_field(points_r),
                get_field(weights_r),
                get_const_field(breaks_r));
        gl_coeffs_theta.compute_points_and_weights_on_mesh(
                get_field(points_theta),
                get_field(weights_theta),
                get_const_field(breaks_theta));

        std::vector<double> vect_points_r(points_r.size());
        for (IdxQuadratureR i : quadrature_idx_range_r) {
            vect_points_r[i - quadrature_idx_range_r.front()] = points_r(i);
        }
        std::vector<double> vect_points_theta(points_theta.size());
        for (IdxQuadratureTheta i : quadrature_idx_range_theta) {
            vect_points_theta[i - quadrature_idx_range_theta.front()] = points_theta(i);
        }

        // Create quadrature index range
        ddc::init_discrete_space<QDimRMesh>(vect_points_r);
        ddc::init_discrete_space<QDimThetaMesh>(vect_points_theta);

        // Find value and derivative of 1D bsplines in radial direction
        ddc::for_each(quadrature_idx_range_r, [&](IdxQuadratureR const ir) {
            std::array<double, 2 * n_non_zero_bases_r> data;
            DSpan2D vals(data.data(), n_non_zero_bases_r, 2);
            ddc::discrete_space<BSplinesR_Polar>()
                    .eval_basis_and_n_derivs(vals, ddc::coordinate(ir), 1);
            for (auto ib : non_zero_bases_r) {
                const int ib_idx = ib - non_zero_bases_r.front();
                r_basis_vals_and_derivs(ib, ir).value = vals(ib_idx, 0);
                r_basis_vals_and_derivs(ib, ir).derivative = vals(ib_idx, 1);
            }
        });

        // Find value and derivative of 1D bsplines in poloidal direction
        ddc::for_each(quadrature_idx_range_theta, [&](IdxQuadratureTheta const ip) {
            std::array<double, 2 * n_non_zero_bases_theta> data;
            DSpan2D vals(data.data(), n_non_zero_bases_theta, 2);
            ddc::discrete_space<BSplinesTheta_Polar>()
                    .eval_basis_and_n_derivs(vals, ddc::coordinate(ip), 1);
            for (auto ib : non_zero_bases_theta) {
                const int ib_idx = ib - non_zero_bases_theta.front();
                theta_basis_vals_and_derivs(ib, ip).value = vals(ib_idx, 0);
                theta_basis_vals_and_derivs(ib, ip).derivative = vals(ib_idx, 1);
            }
        });

        IdxRangeBSPolar singular_idx_range
                = PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>();

        // Find value and derivative of 2D bsplines covering the singular point
        ddc::for_each(quadrature_idx_range_singular, [&](IdxQuadratureRTheta const irp) {
            std::array<double, PolarBSplinesRTheta::n_singular_basis()> singular_data;
            std::array<double, n_non_zero_bases_r * n_non_zero_bases_theta> data;
            // Values of the polar basis splines around the singular point
            // at a given coordinate
            DSpan1D singular_vals(singular_data.data(), PolarBSplinesRTheta::n_singular_basis());
            // Values of the polar basis splines, that do not cover the singular point,
            // at a given coordinate
            DSpan2D vals(data.data(), n_non_zero_bases_r, n_non_zero_bases_theta);
            IdxQuadratureR ir = ddc::select<QDimRMesh>(irp);
            IdxQuadratureTheta ip = ddc::select<QDimThetaMesh>(irp);

            const CoordRTheta coord(ddc::coordinate(irp));

            // Calculate the value
            ddc::discrete_space<PolarBSplinesRTheta>().eval_basis(singular_vals, vals, coord);
            for (auto ib : singular_idx_range) {
                singular_basis_vals_and_derivs(ib, ir, ip).value
                        = singular_vals[ib - singular_idx_range.front()];
            }

            // Calculate the radial derivative
            ddc::discrete_space<PolarBSplinesRTheta>().eval_deriv_r(singular_vals, vals, coord);
            for (auto ib : singular_idx_range) {
                singular_basis_vals_and_derivs(ib, ir, ip).radial_derivative
                        = singular_vals[ib - singular_idx_range.front()];
            }

            // Calculate the poloidal derivative
            ddc::discrete_space<PolarBSplinesRTheta>().eval_deriv_theta(singular_vals, vals, coord);
            for (auto ib : singular_idx_range) {
                singular_basis_vals_and_derivs(ib, ir, ip).poloidal_derivative
                        = singular_vals[ib - singular_idx_range.front()];
            }
        });

        // Find the integral volume associated with each point used in the quadrature scheme
        IdxRangeQuadratureRTheta
                all_quad_points(quadrature_idx_range_r, quadrature_idx_range_theta);
        ddc::for_each(all_quad_points, [&](IdxQuadratureRTheta const irp) {
            IdxQuadratureR const ir = ddc::select<QDimRMesh>(irp);
            IdxQuadratureTheta const ip = ddc::select<QDimThetaMesh>(irp);
            CoordRTheta coord(ddc::coordinate(irp));
            int_volume(ir, ip) = abs(mapping.jacobian(coord)) * weights_r(ir) * weights_theta(ip);
        });

        ddc::NullExtrapolationRule r_extrapolation_rule;
        ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
        SplineRThetaEvaluatorNullBound spline_evaluator(
                r_extrapolation_rule,
                r_extrapolation_rule,
                theta_extrapolation_rule,
                theta_extrapolation_rule);

        // Number of elements in the matrix that correspond to the splines
        // that cover the singular point
        constexpr int n_elements_singular
                = PolarBSplinesRTheta::n_singular_basis() * PolarBSplinesRTheta::n_singular_basis();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // polar splines at the singular point and the other splines
        const int n_elements_overlap = 2
                                       * (PolarBSplinesRTheta::n_singular_basis()
                                          * BSplinesR_Polar::degree() * nbasis_theta);
        const int n_stencil_theta
                = nbasis_theta * min(int(1 + 2 * BSplinesTheta_Polar::degree()), nbasis_theta);
        const int n_stencil_r = nbasis_r * (1 + 2 * BSplinesR_Polar::degree())
                                - (1 + BSplinesR_Polar::degree()) * BSplinesR_Polar::degree();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // non-central splines. These have a tensor product structure
        const int n_elements_stencil = n_stencil_r * n_stencil_theta;

        const int batch_size = 1;
        // Matrix size is equal to the number Polar bspline
        const int matrix_size = ddc::discrete_space<PolarBSplinesRTheta>().nbasis() - nbasis_theta;
        const int n_matrix_elements = n_elements_singular + n_elements_overlap + n_elements_stencil;

        //CSR data storage
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                values_csr_host("values_csr", batch_size, n_matrix_elements);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                col_idx_csr_host("idx_csr", n_matrix_elements);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                nnz_per_row_csr_host("nnz_per_row_csr", matrix_size + 1);

        init_nnz_per_line(nnz_per_row_csr_host);

        // Calculate the matrix elements corresponding to the bsplines which cover the singular point
        ddc::for_each(singular_idx_range, [&](IdxPolarBspl const idx_test) {
            ddc::for_each(singular_idx_range, [&](IdxPolarBspl const idx_trial) {
                // Calculate the weak integral
                double const element = ddc::transform_reduce(
                        quadrature_idx_range_singular,
                        0.0,
                        ddc::reducer::sum<double>(),
                        [&](IdxQuadratureRTheta const quad_idx) {
                            IdxQuadratureR const ir = ddc::select<QDimRMesh>(quad_idx);
                            IdxQuadratureTheta const ip = ddc::select<QDimThetaMesh>(quad_idx);
                            return weak_integral_element(
                                    ir,
                                    ip,
                                    singular_basis_vals_and_derivs(idx_test, ir, ip),
                                    singular_basis_vals_and_derivs(idx_trial, ir, ip),
                                    coeff_alpha,
                                    coeff_beta,
                                    spline_evaluator,
                                    mapping);
                        });
                const int row_idx = idx_test - singular_idx_range.front();
                const int col_idx = idx_trial - singular_idx_range.front();
                const int csr_idx_singular_area = nnz_per_row_csr_host(row_idx + 1);
                //Fill the C matrix corresponding to the splines on the singular point
                col_idx_csr_host(csr_idx_singular_area) = col_idx;
                values_csr_host(batch_idx, csr_idx_singular_area) = element;
                nnz_per_row_csr_host(row_idx + 1)++;
            });
        });

        // Create index ranges associated with the 2D splines
        IdxRangeBSR central_radial_bspline_idx_range(
                radial_bsplines.take_first(IdxStep<BSplinesR_Polar> {BSplinesR_Polar::degree()}));

        IdxRangeBSRTheta non_singular_idx_range_near_centre(
                central_radial_bspline_idx_range,
                polar_bsplines);

        // Calculate the matrix elements where bspline products overlap the bsplines which cover the singular point
        ddc::for_each(singular_idx_range, [&](IdxPolarBspl const idx_test) {
            ddc::for_each(
                    non_singular_idx_range_near_centre,
                    [&](IdxBSpline2D_Polar const idx_trial) {
                        const IdxPolarBspl polar_idx_trial(
                                PolarBSplinesRTheta::get_polar_index<PolarBSplinesRTheta>(
                                        idx_trial));
                        const Idx<BSplinesR_Polar> r_idx_trial(
                                ddc::select<BSplinesR_Polar>(idx_trial));
                        const Idx<BSplinesTheta_Polar> theta_idx_trial(
                                ddc::select<BSplinesTheta_Polar>(idx_trial));

                        // Find the index range covering the cells where both the test and trial functions are non-zero
                        const Idx<RCellDim> first_overlap_element_r(
                                r_idx_trial.uid() < BSplinesR_Polar::degree()
                                        ? 0
                                        : r_idx_trial.uid() - BSplinesR_Polar::degree());
                        const Idx<ThetaCellDim> first_overlap_element_theta(
                                theta_mod(theta_idx_trial.uid() - BSplinesTheta_Polar::degree()));

                        const IdxStep<RCellDim> n_overlap_r(
                                n_overlap_cells - first_overlap_element_r.uid());
                        const IdxStep<ThetaCellDim> n_overlap_theta(
                                BSplinesTheta_Polar::degree() + 1);

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
                                Idx<RBasisSubset> ib_trial_r(r_idx_trial.uid() - cell_idx_r);
                                Idx<ThetaBasisSubset> ib_trial_theta(
                                        theta_mod(theta_idx_trial.uid() - cell_idx_theta));
                                // Calculate the weak integral
                                element += ddc::transform_reduce(
                                        cell_quad_points,
                                        0.0,
                                        ddc::reducer::sum<double>(),
                                        [&](IdxQuadratureRTheta const quad_idx) {
                                            IdxQuadratureR const ir
                                                    = ddc::select<QDimRMesh>(quad_idx);
                                            IdxQuadratureTheta const ip
                                                    = ddc::select<QDimThetaMesh>(quad_idx);
                                            return weak_integral_element<Mapping>(
                                                    ir,
                                                    ip,
                                                    singular_basis_vals_and_derivs(
                                                            idx_test,
                                                            ir,
                                                            ip),
                                                    r_basis_vals_and_derivs(ib_trial_r, ir),
                                                    theta_basis_vals_and_derivs(ib_trial_theta, ip),
                                                    coeff_alpha,
                                                    coeff_beta,
                                                    spline_evaluator,
                                                    mapping);
                                        });
                            });

                            int const row_idx = idx_test - singular_idx_range.front();
                            int const col_idx = polar_idx_trial - singular_idx_range.front();
                            //a_ij
                            col_idx_csr_host(nnz_per_row_csr_host(row_idx + 1)) = col_idx;
                            values_csr_host(batch_idx, nnz_per_row_csr_host(row_idx + 1)) = element;
                            nnz_per_row_csr_host(row_idx + 1)++;
                            //a_ji
                            col_idx_csr_host(nnz_per_row_csr_host(col_idx + 1)) = row_idx;
                            values_csr_host(batch_idx, nnz_per_row_csr_host(col_idx + 1)) = element;
                            nnz_per_row_csr_host(col_idx + 1)++;
                        }
                    });
        });

        // Calculate the matrix elements following a stencil
        ddc::for_each(fem_non_singular_idx_range, [&](IdxPolarBspl const polar_idx_test) {
            const IdxBSpline2D_Polar idx_test(PolarBSplinesRTheta::get_2d_index(polar_idx_test));
            const std::size_t r_idx_test(
                    ddc::select<PolarBSplinesRTheta::BSplinesR_tag>(idx_test).uid());
            const std::size_t theta_idx_test(
                    ddc::select<PolarBSplinesRTheta::BSplinesTheta_tag>(idx_test).uid());

            // Calculate the index of the elements that are already filled
            IdxRangeBSTheta_Polar remaining_theta(
                    Idx<BSplinesTheta_Polar> {theta_idx_test},
                    IdxStep<BSplinesTheta_Polar> {BSplinesTheta_Polar::degree() + 1});
            ddc::for_each(remaining_theta, [&](auto const theta_idx_trial) {
                IdxBSpline2D_Polar idx_trial(Idx<BSplinesR_Polar>(r_idx_test), theta_idx_trial);
                IdxPolarBspl polar_idx_trial(
                        PolarBSplinesRTheta::get_polar_index<PolarBSplinesRTheta>(
                                IdxBSpline2D_Polar(r_idx_test, theta_mod(theta_idx_trial.uid()))));
                double element = get_matrix_stencil_element(
                        idx_test,
                        idx_trial,
                        coeff_alpha,
                        coeff_beta,
                        spline_evaluator,
                        mapping);
                int const int_polar_idx_test = polar_idx_test - singular_idx_range.front();
                if (polar_idx_test == polar_idx_trial) {
                    const int idx = nnz_per_row_csr_host(int_polar_idx_test + 1);
                    col_idx_csr_host(idx) = int_polar_idx_test;
                    values_csr_host(batch_idx, idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_test + 1)++;
                } else {
                    int const int_polar_idx_trial = polar_idx_trial - singular_idx_range.front();

                    const int aij_idx = nnz_per_row_csr_host(int_polar_idx_test + 1);
                    col_idx_csr_host(aij_idx) = int_polar_idx_trial;
                    values_csr_host(batch_idx, aij_idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_test + 1)++;

                    const int aji_idx = nnz_per_row_csr_host(int_polar_idx_trial + 1);
                    col_idx_csr_host(aji_idx) = int_polar_idx_test;
                    values_csr_host(batch_idx, aji_idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_trial + 1)++;
                }
            });
            IdxRangeBSR remaining_r(
                    ddc::select<BSplinesR_Polar>(idx_test) + 1,
                    IdxStep<BSplinesR_Polar> {
                            min(BSplinesR_Polar::degree(),
                                ddc::discrete_space<BSplinesR_Polar>().nbasis() - 2 - r_idx_test)});
            IdxRangeBSTheta_Polar relevant_theta(
                    Idx<BSplinesTheta_Polar> {
                            theta_idx_test + ddc::discrete_space<BSplinesTheta_Polar>().nbasis()
                            - BSplinesTheta_Polar::degree()},
                    IdxStep<BSplinesTheta_Polar> {2 * BSplinesTheta_Polar::degree() + 1});

            IdxRangeBSRTheta trial_idx_range(remaining_r, relevant_theta);

            ddc::for_each(trial_idx_range, [&](IdxBSpline2D_Polar const idx_trial) {
                const int r_idx_trial(
                        ddc::select<PolarBSplinesRTheta::BSplinesR_tag>(idx_trial).uid());
                const int theta_idx_trial(
                        ddc::select<PolarBSplinesRTheta::BSplinesTheta_tag>(idx_trial).uid());
                IdxPolarBspl polar_idx_trial(
                        PolarBSplinesRTheta::get_polar_index<PolarBSplinesRTheta>(
                                IdxBSpline2D_Polar(r_idx_trial, theta_mod(theta_idx_trial))));
                double element = get_matrix_stencil_element(
                        idx_test,
                        idx_trial,
                        coeff_alpha,
                        coeff_beta,
                        spline_evaluator,
                        mapping);
                int const int_polar_idx_test = polar_idx_test - singular_idx_range.front();
                if (polar_idx_test == polar_idx_trial) {
                    const int idx = nnz_per_row_csr_host(int_polar_idx_test + 1);
                    col_idx_csr_host(idx) = int_polar_idx_test;
                    values_csr_host(batch_idx, idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_test + 1)++;
                } else {
                    int const int_polar_idx_trial = polar_idx_trial - singular_idx_range.front();
                    const int aij_idx = nnz_per_row_csr_host(int_polar_idx_test + 1);
                    col_idx_csr_host(aij_idx) = int_polar_idx_trial;
                    values_csr_host(batch_idx, aij_idx) = element;
                    nnz_per_row_csr_host(int_polar_idx_test + 1)++;

                    const int aji_idx = nnz_per_row_csr_host(int_polar_idx_trial + 1);
                    col_idx_csr_host(aji_idx) = int_polar_idx_test;
                    values_csr_host(batch_idx, aji_idx) = element;
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
     * @param[out] spline
     *      The spline representation of the solution @f$\phi@f$.
     */
    template <class RHSFunction>
    void operator()(RHSFunction const& rhs, SplinePolar& spline) const
    {
        const int b_size = ddc::discrete_space<PolarBSplinesRTheta>().nbasis()
                           - ddc::discrete_space<BSplinesTheta_Polar>().nbasis();
        const int batch_size = 1;
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                b_host("b_host", batch_size, b_size);

        // Fill b
        ddc::for_each(
                PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>(),
                [&](IdxPolarBspl const idx) {
                    b_host(0, idx.uid()) = ddc::transform_reduce(
                            quadrature_idx_range_singular,
                            0.0,
                            ddc::reducer::sum<double>(),
                            [&](IdxQuadratureRTheta const quad_idx) {
                                IdxQuadratureR const ir = ddc::select<QDimRMesh>(quad_idx);
                                IdxQuadratureTheta const ip = ddc::select<QDimThetaMesh>(quad_idx);
                                CoordRTheta coord(ddc::coordinate(quad_idx));
                                return rhs(coord)
                                       * singular_basis_vals_and_derivs(idx, ir, ip).value
                                       * int_volume(ir, ip);
                            });
                });
        const std::size_t ncells_r = ddc::discrete_space<BSplinesR_Polar>().ncells();
        ddc::for_each(fem_non_singular_idx_range, [&](IdxPolarBspl const idx) {
            const IdxBSpline2D_Polar idx_2d(PolarBSplinesRTheta::get_2d_index(idx));
            const std::size_t r_idx(ddc::select<PolarBSplinesRTheta::BSplinesR_tag>(idx_2d).uid());
            const std::size_t theta_idx(
                    ddc::select<PolarBSplinesRTheta::BSplinesTheta_tag>(idx_2d).uid());

            // Find the cells on which the bspline is non-zero
            int first_cell_r(r_idx - BSplinesR_Polar::degree());
            int first_cell_theta(theta_idx - BSplinesTheta_Polar::degree());
            std::size_t last_cell_r(r_idx + 1);
            if (first_cell_r < 0)
                first_cell_r = 0;
            if (last_cell_r > ncells_r)
                last_cell_r = ncells_r;
            IdxStep<RCellDim> const r_length(last_cell_r - first_cell_r);
            IdxStep<ThetaCellDim> const theta_length(BSplinesTheta_Polar::degree() + 1);


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
                Idx<RBasisSubset> ib_r(r_idx - cell_idx_r);
                Idx<ThetaBasisSubset> ib_theta(theta_mod(theta_idx - cell_idx_theta));

                // Calculate the weak integral
                element += ddc::transform_reduce(
                        cell_quad_points,
                        0.0,
                        ddc::reducer::sum<double>(),
                        [&](IdxQuadratureRTheta const quad_idx) {
                            IdxQuadratureR const ir = ddc::select<QDimRMesh>(quad_idx);
                            IdxQuadratureTheta const ip = ddc::select<QDimThetaMesh>(quad_idx);
                            CoordRTheta coord(ddc::coordinate(quad_idx));
                            double rb = r_basis_vals_and_derivs(ib_r, ir).value;
                            double pb = theta_basis_vals_and_derivs(ib_theta, ip).value;
                            return rhs(coord) * rb * pb * int_volume(ir, ip);
                        });
            });
            b_host(0, idx.uid()) = element;
        });

        // Solve the matrix equation
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
                b("b", batch_size, b_size);
        Kokkos::deep_copy(b, b_host);

        m_gko_matrix->solve(b);

        Kokkos::deep_copy(b_host, b);
        //-----------------
        IdxRangeBSRTheta dirichlet_boundary_idx_range(
                radial_bsplines.take_last(IdxStep<BSplinesR_Polar> {1}),
                polar_bsplines);
        IdxRangeBSTheta_Polar polar_idx_range(
                ddc::discrete_space<BSplinesTheta_Polar>().full_domain());


        // Fill the spline
        ddc::for_each(
                PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>(),
                [&](IdxPolarBspl const idx) {
                    const int bspl_idx
                            = idx
                              - PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>()
                                        .front();
                    spline.singular_spline_coef(idx) = b_host(0, bspl_idx);
                });
        ddc::for_each(fem_non_singular_idx_range, [&](IdxPolarBspl const idx) {
            const IdxBSpline2D_Polar idx_2d(PolarBSplinesRTheta::get_2d_index(idx));
            spline.spline_coef(idx_2d) = b_host(0, idx.uid());
        });
        ddc::for_each(dirichlet_boundary_idx_range, [&](IdxBSpline2D_Polar const idx) {
            spline.spline_coef(idx) = 0.0;
        });

        // Copy the periodic elements
        IdxRangeBSRTheta copy_idx_range(
                radial_bsplines,
                polar_idx_range.remove_first(IdxStep<BSplinesTheta_Polar>(
                        ddc::discrete_space<BSplinesTheta_Polar>().nbasis())));
        ddc::for_each(copy_idx_range, [&](IdxBSpline2D_Polar const idx_2d) {
            spline.spline_coef(
                    ddc::select<BSplinesR_Polar>(idx_2d),
                    ddc::select<BSplinesTheta_Polar>(idx_2d))
                    = spline.spline_coef(
                            ddc::select<BSplinesR_Polar>(idx_2d),
                            ddc::select<BSplinesTheta_Polar>(idx_2d)
                                    - ddc::discrete_space<BSplinesTheta_Polar>().nbasis());
        });
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
     * @param[in] coords_eval
     *      A Field of coordinates where we want to compute the solution.
     * @param[out] result
     *      The values of the solution @f$\phi@f$ on the given coords_eval.
     */
    template <class RHSFunction>
    void operator()(
            RHSFunction const& rhs,
            host_t<ConstFieldRTheta<CoordRTheta>> const coords_eval,
            host_t<DFieldRTheta> result) const
    {
        IdxRangeBSTheta_Polar polar_idx_range(
                ddc::discrete_space<BSplinesTheta_Polar>().full_domain());
        SplinePolar
                spline(PolarBSplinesRTheta::singular_idx_range<PolarBSplinesRTheta>(),
                       IdxRangeBSRTheta(radial_bsplines, polar_idx_range));

        (*this)(rhs, spline);
        m_polar_spline_evaluator(result, coords_eval, spline);
    }


private:
    static IdxRangeQuadratureRTheta get_quadrature_points_in_cell(
            int cell_idx_r,
            int cell_idx_theta)
    {
        const IdxQuadratureR first_quad_point_r(cell_idx_r * n_gauss_legendre_r);
        const IdxQuadratureTheta first_quad_point_theta(cell_idx_theta * n_gauss_legendre_theta);
        constexpr IdxStepQuadratureR n_GL_r(n_gauss_legendre_r);
        constexpr IdxStepQuadratureTheta n_GL_theta(n_gauss_legendre_theta);
        const IdxRangeQuadratureR quad_points_r(first_quad_point_r, n_GL_r);
        const IdxRangeQuadratureTheta quad_points_theta(first_quad_point_theta, n_GL_theta);
        return IdxRangeQuadratureRTheta(quad_points_r, quad_points_theta);
    }

    template <class Mapping>
    double weak_integral_element(
            IdxQuadratureR ir,
            IdxQuadratureTheta ip,
            EvalDeriv2DType const& test_bspline_val_and_deriv,
            EvalDeriv2DType const& trial_bspline_val_and_deriv,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                ir,
                ip,
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
            IdxQuadratureR ir,
            IdxQuadratureTheta ip,
            EvalDeriv2DType const& test_bspline_val_and_deriv,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_r,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_theta,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                ir,
                ip,
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
            IdxQuadratureR ir,
            IdxQuadratureTheta ip,
            EvalDeriv1DType const& test_bspline_val_and_deriv_r,
            EvalDeriv2DType const& trial_bspline_val_and_deriv,
            EvalDeriv1DType const& test_bspline_val_and_deriv_theta,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                ir,
                ip,
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
            IdxQuadratureR ir,
            IdxQuadratureTheta ip,
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
                ir,
                ip,
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
            IdxQuadratureR ir,
            IdxQuadratureTheta ip,
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
        CoordRTheta coord(ddc::coordinate(ir), ddc::coordinate(ip));
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
        return int_volume(ir, ip)
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
            IdxBSpline2D_Polar idx_test,
            IdxBSpline2D_Polar idx_trial,
            host_t<ConstSpline2D> coeff_alpha,
            host_t<ConstSpline2D> coeff_beta,
            SplineRThetaEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        // 0 <= r_idx_test < 8
        // 0 <= r_idx_trial < 8
        // r_idx_test < r_idx_trial
        const int r_idx_test(ddc::select<BSplinesR_Polar>(idx_test).uid());
        const int r_idx_trial(ddc::select<BSplinesR_Polar>(idx_trial).uid());
        // 0 <= theta_idx_test < 8
        // 0 <= theta_idx_trial < 8
        int theta_idx_test(theta_mod(ddc::select<BSplinesTheta_Polar>(idx_test).uid()));
        int theta_idx_trial(theta_mod(ddc::select<BSplinesTheta_Polar>(idx_trial).uid()));

        const std::size_t ncells_r = ddc::discrete_space<BSplinesR_Polar>().ncells();

        // 0<= r_offset <= degree_r
        // -degree_theta <= theta_offset <= degree_theta
        const int r_offset = r_idx_trial - r_idx_test;
        int theta_offset = theta_mod(theta_idx_trial - theta_idx_test);
        if (theta_offset >= int(nbasis_theta - BSplinesTheta_Polar::degree())) {
            theta_offset -= nbasis_theta;
        }
        assert(r_offset >= 0);
        assert(r_offset <= int(BSplinesR_Polar::degree()));
        assert(theta_offset >= -int(BSplinesTheta_Polar::degree()));
        assert(theta_offset <= int(BSplinesTheta_Polar::degree()));

        // Find the index range covering the cells where both the test and trial functions are non-zero
        int n_overlap_stencil_r(BSplinesR_Polar::degree() + 1 - r_offset);
        int first_overlap_r(r_idx_trial - BSplinesR_Polar::degree());

        int first_overlap_theta;
        int n_overlap_stencil_theta;
        if (theta_offset > 0) {
            n_overlap_stencil_theta = BSplinesTheta_Polar::degree() + 1 - theta_offset;
            first_overlap_theta = theta_mod(theta_idx_trial - BSplinesTheta_Polar::degree());
        } else {
            n_overlap_stencil_theta = BSplinesTheta_Polar::degree() + 1 + theta_offset;
            first_overlap_theta = theta_mod(theta_idx_test - BSplinesTheta_Polar::degree());
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

                    int ib_test_theta_idx = theta_idx_test - cell_idx_theta;
                    int ib_trial_theta_idx = theta_idx_trial - cell_idx_theta;

                    // Find the column where the non-zero data is stored
                    Idx<RBasisSubset> ib_test_r(r_idx_test - cell_idx_r);
                    Idx<ThetaBasisSubset> ib_test_theta(theta_mod(ib_test_theta_idx));
                    Idx<RBasisSubset> ib_trial_r(r_idx_trial - cell_idx_r);
                    Idx<ThetaBasisSubset> ib_trial_theta(theta_mod(ib_trial_theta_idx));

                    assert(ib_test_r.uid() < BSplinesR_Polar::degree() + 1);
                    assert(ib_test_theta.uid() < BSplinesTheta_Polar::degree() + 1);
                    assert(ib_trial_r.uid() < BSplinesR_Polar::degree() + 1);
                    assert(ib_trial_theta.uid() < BSplinesTheta_Polar::degree() + 1);

                    // Calculate the weak integral
                    return ddc::transform_reduce(
                            cell_quad_points,
                            0.0,
                            ddc::reducer::sum<double>(),
                            [&](IdxQuadratureRTheta const quad_idx) {
                                IdxQuadratureR const r_idx = ddc::select<QDimRMesh>(quad_idx);
                                IdxQuadratureTheta const theta_idx
                                        = ddc::select<QDimThetaMesh>(quad_idx);
                                return weak_integral_element(
                                        r_idx,
                                        theta_idx,
                                        r_basis_vals_and_derivs(ib_test_r, r_idx),
                                        r_basis_vals_and_derivs(ib_trial_r, r_idx),
                                        theta_basis_vals_and_derivs(ib_test_theta, theta_idx),
                                        theta_basis_vals_and_derivs(ib_trial_theta, theta_idx),
                                        coeff_alpha,
                                        coeff_beta,
                                        evaluator,
                                        mapping);
                            });
                });
    }

    static int theta_mod(int idx_theta)
    {
        int ncells_theta = ddc::discrete_space<BSplinesTheta_Polar>().ncells();
        while (idx_theta < 0)
            idx_theta += ncells_theta;
        while (idx_theta >= ncells_theta)
            idx_theta -= ncells_theta;
        return idx_theta;
    }

    void init_nnz_per_line(
            Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> nnz)
    {
        size_t const mat_size = nnz.extent(0) - 1;
        size_t constexpr n_singular_basis = PolarBSplinesRTheta::n_singular_basis();
        size_t constexpr degree = BSplinesR_Polar::degree();
        size_t constexpr radial_overlap = 2 * degree + 1;
        // overlap between singular domain splines and radial splines
        for (size_t k = 1; k <= n_singular_basis; k++) {
            nnz(k + 1) = n_singular_basis + degree * nbasis_theta;
        }
        // going from the internal boundary the overlapping possiblities between two radial splines increase
        for (size_t i = 1; i <= degree + 1; i++) {
            for (size_t k = n_singular_basis + (i - 1) * nbasis_theta;
                 k < n_singular_basis + i * nbasis_theta;
                 k++) {
                nnz(k + 2) = n_singular_basis + (degree + i) * radial_overlap;
            }
        }
        // Stencil with maximum possible overlap from two sides for radial spline
        for (size_t k = n_singular_basis + degree * nbasis_theta;
             k < mat_size - degree * nbasis_theta;
             k++) {
            nnz(k + 2) = radial_overlap * radial_overlap;
        }
        // Approaching the external boundary the overlapping possiblities between two radial splines decrease
        for (size_t i = 1; i <= degree; i++) {
            for (size_t k = mat_size - i * nbasis_theta; k < mat_size - (i - 1) * nbasis_theta;
                 k++) {
                nnz(k + 2) = (degree + i) * radial_overlap;
            }
        }

        // sum non-zero elements count
        for (size_t k = 1; k < mat_size; k++) {
            nnz(k + 1) += nnz(k);
        }
        nnz(0) = 0;
        nnz(1) = 0;
    }
};