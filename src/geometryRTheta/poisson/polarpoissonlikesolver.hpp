#pragma once

#include <iomanip>

#include <ddc/ddc.hpp>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/math_tools.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>

#include <Eigen/Sparse>

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
    struct PBasisSubset
    {
    };
    struct RCellDim
    {
    };
    struct PCellDim
    {
    };


public:
    /**
     * @brief Tag the first dimension for the quadrature mesh.
     */
    struct QDimRMesh : ddc::NonUniformPointSampling<RDimR>
    {
    };
    /**
     * @brief Tag the second dimension for the quadrature mesh.
     */
    struct QDimPMesh : ddc::NonUniformPointSampling<RDimP>
    {
    };

private:
    /**
     * @brief Tag the quadrature domain in the first dimension.
     */
    using QuadratureDomainR = ddc::DiscreteDomain<QDimRMesh>;
    /**
     * @brief Tag the quadrature domain in the second dimension.
     */
    using QuadratureDomainP = ddc::DiscreteDomain<QDimPMesh>;
    /**
     * @brief Tag the quadrature domain.
     */
    using QuadratureDomainRP = ddc::DiscreteDomain<QDimRMesh, QDimPMesh>;
    /**
     * @brief Tag the elements (index) of the quadrature domain in the first dimension.
     */
    using QuadratureIndexR = ddc::DiscreteElement<QDimRMesh>;
    /**
     * @brief Tag the elements (index) of the quadrature domain in the second dimension.
     */
    using QuadratureIndexP = ddc::DiscreteElement<QDimPMesh>;
    /**
     * @brief Tag the elements (index) of the quadrature domain.
     */
    using QuadratureIndexRP = ddc::DiscreteElement<QDimRMesh, QDimPMesh>;
    /**
     * @brief Tag a vector on the first dimension of the quadrature mesh.
     */
    using QuadratureLengthR = ddc::DiscreteVector<QDimRMesh>;
    /**
     * @brief Tag a vector on the second dimension of the quadrature mesh.
     */
    using QuadratureLengthP = ddc::DiscreteVector<QDimPMesh>;

    using BSplinesR_Polar = PolarBSplinesRP::BSplinesR_tag;
    using BSplinesP_Polar = PolarBSplinesRP::BSplinesP_tag;

    using IDimBSpline2D_Polar = ddc::DiscreteElement<BSplinesR_Polar, BSplinesP_Polar>;

    using BSDomainR_Polar = ddc::DiscreteDomain<BSplinesR_Polar>;
    using BSDomainP_Polar = ddc::DiscreteDomain<BSplinesP_Polar>;

    using KnotsR = ddc::NonUniformBsplinesKnots<BSplinesR_Polar>;
    using KnotsP = ddc::NonUniformBsplinesKnots<BSplinesP_Polar>;
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
    using CellIndex = ddc::DiscreteElement<RCellDim, PCellDim>;

    /**
     * @brief Tag type of matrix element.
     */
    using MatrixElement = Eigen::Triplet<double>;

private:
    static constexpr int n_gauss_legendre_r = BSplinesR_Polar::degree() + 1;
    static constexpr int n_gauss_legendre_p = BSplinesP_Polar::degree() + 1;
    // The number of cells (in the radial direction) in which both types of basis splines can be found
    static constexpr int n_overlap_cells = PolarBSplinesRP::continuity + 1;

    // Number of cells over which a radial B-splines has its support
    // This is the case for b-splines which are not affected by the higher knot multiplicity at the boundary.
    static constexpr ddc::DiscreteVector<RBasisSubset> n_non_zero_bases_r
            = ddc::DiscreteVector<RBasisSubset>(BSplinesR_Polar::degree() + 1);

    // Number of cells over which a poloidal B-splines has its support
    static constexpr ddc::DiscreteVector<PBasisSubset> n_non_zero_bases_p
            = ddc::DiscreteVector<PBasisSubset>(BSplinesP_Polar::degree() + 1);

    static constexpr ddc::DiscreteDomain<RBasisSubset> non_zero_bases_r = ddc::DiscreteDomain<
            RBasisSubset>(ddc::DiscreteElement<RBasisSubset> {0}, n_non_zero_bases_r);
    static constexpr ddc::DiscreteDomain<PBasisSubset> non_zero_bases_p = ddc::DiscreteDomain<
            PBasisSubset>(ddc::DiscreteElement<PBasisSubset> {0}, n_non_zero_bases_p);

    const int nbasis_r;
    const int nbasis_p;

    // Domains
    BSDomainPolar fem_non_singular_domain;
    BSDomainR_Polar radial_bsplines;
    BSDomainP_Polar polar_bsplines;

    QuadratureDomainR quadrature_domain_r;
    QuadratureDomainP quadrature_domain_p;
    QuadratureDomainRP quadrature_domain_singular;

    // Gauss-Legendre points and weights
    ddc::Chunk<ddc::Coordinate<RDimR>, QuadratureDomainR> points_r;
    ddc::Chunk<ddc::Coordinate<RDimP>, QuadratureDomainP> points_p;
    ddc::Chunk<double, QuadratureDomainR> weights_r;
    ddc::Chunk<double, QuadratureDomainP> weights_p;

    // Basis Spline values and derivatives at Gauss-Legendre points
    ddc::Chunk<EvalDeriv2DType, ddc::DiscreteDomain<PolarBSplinesRP, QDimRMesh, QDimPMesh>>
            singular_basis_vals_and_derivs;
    ddc::Chunk<EvalDeriv1DType, ddc::DiscreteDomain<RBasisSubset, QDimRMesh>>
            r_basis_vals_and_derivs;
    ddc::Chunk<EvalDeriv1DType, ddc::DiscreteDomain<PBasisSubset, QDimPMesh>>
            p_basis_vals_and_derivs;

    ddc::Chunk<double, QuadratureDomainRP> int_volume;

    PolarSplineEvaluator<PolarBSplinesRP, ddc::NullExtrapolationRule> m_polar_spline_evaluator;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_matrix;

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
     *      The mapping from the logical domain to the physical domain where
     *      the equation is defined.
     *
     * @tparam Mapping A Curvilinear2DToCartesian class.
     */
    template <class Mapping>
    PolarSplineFEMPoissonLikeSolver(
            Spline2DView coeff_alpha,
            Spline2DView coeff_beta,
            Mapping const& mapping)
        : nbasis_r(ddc::discrete_space<BSplinesR_Polar>().nbasis() - n_overlap_cells - 1)
        , nbasis_p(ddc::discrete_space<BSplinesP_Polar>().nbasis())
        , fem_non_singular_domain(
                  ddc::discrete_space<PolarBSplinesRP>().tensor_bspline_domain().remove_last(
                          ddc::DiscreteVector<PolarBSplinesRP> {nbasis_p}))
        , radial_bsplines(ddc::discrete_space<BSplinesR_Polar>().full_domain().remove_first(
                  ddc::DiscreteVector<BSplinesR_Polar> {n_overlap_cells}))
        , polar_bsplines(ddc::discrete_space<BSplinesP_Polar>().full_domain().take_first(
                  ddc::DiscreteVector<BSplinesP_Polar> {nbasis_p}))
        , quadrature_domain_r(
                  ddc::DiscreteElement<QDimRMesh>(0),
                  ddc::DiscreteVector<QDimRMesh>(
                          n_gauss_legendre_r * ddc::discrete_space<BSplinesR_Polar>().ncells()))
        , quadrature_domain_p(
                  ddc::DiscreteElement<QDimPMesh>(0),
                  ddc::DiscreteVector<QDimPMesh>(
                          n_gauss_legendre_p * ddc::discrete_space<BSplinesP_Polar>().ncells()))
        , quadrature_domain_singular(
                  quadrature_domain_r.take_first(
                          ddc::DiscreteVector<QDimRMesh> {n_overlap_cells * n_gauss_legendre_r}),
                  quadrature_domain_p)
        , points_r(quadrature_domain_r)
        , points_p(quadrature_domain_p)
        , weights_r(quadrature_domain_r)
        , weights_p(quadrature_domain_p)
        , singular_basis_vals_and_derivs(ddc::DiscreteDomain<PolarBSplinesRP, QDimRMesh, QDimPMesh>(
                  PolarBSplinesRP::singular_domain<PolarBSplinesRP>(),
                  ddc::select<QDimRMesh>(quadrature_domain_singular),
                  ddc::select<QDimPMesh>(quadrature_domain_singular)))
        , r_basis_vals_and_derivs(ddc::DiscreteDomain<
                                  RBasisSubset,
                                  QDimRMesh>(non_zero_bases_r, quadrature_domain_r))
        , p_basis_vals_and_derivs(ddc::DiscreteDomain<
                                  PBasisSubset,
                                  QDimPMesh>(non_zero_bases_p, quadrature_domain_p))
        , int_volume(QuadratureDomainRP(quadrature_domain_r, quadrature_domain_p))
        , m_polar_spline_evaluator(ddc::NullExtrapolationRule())
    {
        // Get break points
        ddc::DiscreteDomain<KnotsR> r_edges_dom
                = ddc::discrete_space<BSplinesR_Polar>().break_point_domain();
        ddc::DiscreteDomain<KnotsP> p_edges_dom
                = ddc::discrete_space<BSplinesP_Polar>().break_point_domain();
        ddc::Chunk<ddc::Coordinate<RDimR>, ddc::DiscreteDomain<KnotsR>> breaks_r(r_edges_dom);
        ddc::Chunk<ddc::Coordinate<RDimP>, ddc::DiscreteDomain<KnotsP>> breaks_p(p_edges_dom);

        ddc::for_each(r_edges_dom, [&](ddc::DiscreteElement<KnotsR> i) {
            breaks_r(i) = ddc::coordinate(i);
        });
        ddc::for_each(p_edges_dom, [&](ddc::DiscreteElement<KnotsP> i) {
            breaks_p(i) = ddc::coordinate(i);
        });

        // Define quadrature points and weights
        GaussLegendre<RDimR> gl_coeffs_r(n_gauss_legendre_r);
        GaussLegendre<RDimP> gl_coeffs_p(n_gauss_legendre_p);
        gl_coeffs_r.compute_points_and_weights_on_mesh(
                points_r.span_view(),
                weights_r.span_view(),
                breaks_r.span_cview());
        gl_coeffs_p.compute_points_and_weights_on_mesh(
                points_p.span_view(),
                weights_p.span_view(),
                breaks_p.span_cview());

        std::vector<double> vect_points_r(points_r.size());
        for (auto i : quadrature_domain_r) {
            vect_points_r[i.uid()] = points_r(i);
        }
        std::vector<double> vect_points_p(points_p.size());
        for (auto i : quadrature_domain_p) {
            vect_points_p[i.uid()] = points_p(i);
        }

        // Create quadrature domain
        ddc::init_discrete_space<QDimRMesh>(vect_points_r);
        ddc::init_discrete_space<QDimPMesh>(vect_points_p);

        // Find value and derivative of 1D bsplines in radial direction
        ddc::for_each(quadrature_domain_r, [&](QuadratureIndexR const ir) {
            std::array<double, 2 * n_non_zero_bases_r> data;
            DSpan2D vals(data.data(), n_non_zero_bases_r, 2);
            ddc::discrete_space<BSplinesR_Polar>()
                    .eval_basis_and_n_derivs(vals, ddc::coordinate(ir), 1);
            for (auto ib : non_zero_bases_r) {
                r_basis_vals_and_derivs(ib, ir).value = vals(ib.uid(), 0);
                r_basis_vals_and_derivs(ib, ir).derivative = vals(ib.uid(), 1);
            }
        });

        // Find value and derivative of 1D bsplines in poloidal direction
        ddc::for_each(quadrature_domain_p, [&](QuadratureIndexP const ip) {
            std::array<double, 2 * n_non_zero_bases_p> data;
            DSpan2D vals(data.data(), n_non_zero_bases_p, 2);
            ddc::discrete_space<BSplinesP_Polar>()
                    .eval_basis_and_n_derivs(vals, ddc::coordinate(ip), 1);
            for (auto ib : non_zero_bases_p) {
                p_basis_vals_and_derivs(ib, ip).value = vals(ib.uid(), 0);
                p_basis_vals_and_derivs(ib, ip).derivative = vals(ib.uid(), 1);
            }
        });

        auto singular_domain = PolarBSplinesRP::singular_domain<PolarBSplinesRP>();

        // Find value and derivative of 2D bsplines covering the singular point
        ddc::for_each(quadrature_domain_singular, [&](QuadratureIndexRP const irp) {
            std::array<double, PolarBSplinesRP::n_singular_basis()> singular_data;
            std::array<double, n_non_zero_bases_r * n_non_zero_bases_p> data;
            // Values of the polar basis splines around the singular point
            // at a given coordinate
            DSpan1D singular_vals(singular_data.data(), PolarBSplinesRP::n_singular_basis());
            // Values of the polar basis splines, that do not cover the singular point,
            // at a given coordinate
            DSpan2D vals(data.data(), n_non_zero_bases_r, n_non_zero_bases_p);
            QuadratureIndexR ir = ddc::select<QDimRMesh>(irp);
            QuadratureIndexP ip = ddc::select<QDimPMesh>(irp);

            const CoordRP coord(ddc::coordinate(irp));

            // Calculate the value
            ddc::discrete_space<PolarBSplinesRP>().eval_basis(singular_vals, vals, coord);
            for (auto ib : singular_domain) {
                singular_basis_vals_and_derivs(ib, ir, ip).value = singular_vals[ib.uid()];
            }

            // Calculate the radial derivative
            ddc::discrete_space<PolarBSplinesRP>().eval_deriv_r(singular_vals, vals, coord);
            for (auto ib : singular_domain) {
                singular_basis_vals_and_derivs(ib, ir, ip).radial_derivative
                        = singular_vals[ib.uid()];
            }

            // Calculate the poloidal derivative
            ddc::discrete_space<PolarBSplinesRP>().eval_deriv_p(singular_vals, vals, coord);
            for (auto ib : singular_domain) {
                singular_basis_vals_and_derivs(ib, ir, ip).poloidal_derivative
                        = singular_vals[ib.uid()];
            }
        });

        // Find the integral volume associated with each point used in the quadrature scheme
        QuadratureDomainRP all_quad_points(quadrature_domain_r, quadrature_domain_p);
        ddc::for_each(all_quad_points, [&](QuadratureIndexRP const irp) {
            QuadratureIndexR const ir = ddc::select<QDimRMesh>(irp);
            QuadratureIndexP const ip = ddc::select<QDimPMesh>(irp);
            CoordRP coord(ddc::coordinate(irp));
            int_volume(ir, ip) = abs(mapping.jacobian(coord)) * weights_r(ir) * weights_p(ip);
        });

        ddc::NullExtrapolationRule r_extrapolation_rule;
        ddc::PeriodicExtrapolationRule<RDimP> p_extrapolation_rule;
        SplineRPEvaluatorNullBound spline_evaluator(
                r_extrapolation_rule,
                r_extrapolation_rule,
                p_extrapolation_rule,
                p_extrapolation_rule);

        // Number of elements in the matrix that correspond to the splines
        // that cover the singular point
        constexpr int n_elements_singular
                = PolarBSplinesRP::n_singular_basis() * PolarBSplinesRP::n_singular_basis();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // polar splines at the singular point and the other splines
        const int n_elements_overlap
                = 2 * (PolarBSplinesRP::n_singular_basis() * BSplinesR_Polar::degree() * nbasis_p);
        const int n_stencil_p = nbasis_p * min(int(1 + 2 * BSplinesP_Polar::degree()), nbasis_p);
        const int n_stencil_r = nbasis_r * (1 + 2 * BSplinesR_Polar::degree())
                                - (1 + BSplinesR_Polar::degree()) * BSplinesR_Polar::degree();
        // Number of non-zero elements in the matrix corresponding to the inner product of
        // non-central splines. These have a tensor product structure
        const int n_elements_stencil = n_stencil_r * n_stencil_p;

        // Matrix size is equal to the number Polar bspline
        const int n_matrix_size = ddc::discrete_space<PolarBSplinesRP>().nbasis() - nbasis_p;
        Eigen::SparseMatrix<double> matrix(n_matrix_size, n_matrix_size);
        const int n_matrix_elements = n_elements_singular + n_elements_overlap + n_elements_stencil;
        std::vector<MatrixElement> matrix_elements(n_matrix_elements);
        int matrix_idx(0);
        // Calculate the matrix elements corresponding to the bsplines which cover the singular point
        ddc::for_each(singular_domain, [&](IndexPolarBspl const idx_test) {
            ddc::for_each(singular_domain, [&](IndexPolarBspl const idx_trial) {
                // Calculate the weak integral
                matrix_elements[matrix_idx++] = MatrixElement(
                        idx_test.uid(),
                        idx_trial.uid(),
                        ddc::transform_reduce(
                                quadrature_domain_singular,
                                0.0,
                                ddc::reducer::sum<double>(),
                                [&](QuadratureIndexRP const quad_idx) {
                                    QuadratureIndexR const ir = ddc::select<QDimRMesh>(quad_idx);
                                    QuadratureIndexP const ip = ddc::select<QDimPMesh>(quad_idx);
                                    return weak_integral_element(
                                            ir,
                                            ip,
                                            singular_basis_vals_and_derivs(idx_test, ir, ip),
                                            singular_basis_vals_and_derivs(idx_trial, ir, ip),
                                            coeff_alpha,
                                            coeff_beta,
                                            spline_evaluator,
                                            mapping);
                                }));
            });
        });
        assert(matrix_idx == n_elements_singular);

        // Create domains associated with the 2D splines
        BSDomainR central_radial_bspline_domain(radial_bsplines.take_first(
                ddc::DiscreteVector<BSplinesR_Polar> {BSplinesR_Polar::degree()}));

        BSDomainRP non_singular_domain_near_centre(central_radial_bspline_domain, polar_bsplines);

        // Calculate the matrix elements where bspline products overlap the bsplines which cover the singular point
        ddc::for_each(singular_domain, [&](IndexPolarBspl const idx_test) {
            ddc::for_each(
                    non_singular_domain_near_centre,
                    [&](IDimBSpline2D_Polar const idx_trial) {
                        const IndexPolarBspl polar_idx_trial(
                                PolarBSplinesRP::get_polar_index<PolarBSplinesRP>(idx_trial));
                        const ddc::DiscreteElement<BSplinesR_Polar> r_idx_trial(
                                ddc::select<BSplinesR_Polar>(idx_trial));
                        const ddc::DiscreteElement<BSplinesP_Polar> p_idx_trial(
                                ddc::select<BSplinesP_Polar>(idx_trial));

                        // Find the domain covering the cells where both the test and trial functions are non-zero
                        const ddc::DiscreteElement<RCellDim> first_overlap_element_r(
                                r_idx_trial.uid() < BSplinesR_Polar::degree()
                                        ? 0
                                        : r_idx_trial.uid() - BSplinesR_Polar::degree());
                        const ddc::DiscreteElement<PCellDim> first_overlap_element_p(
                                pmod(p_idx_trial.uid() - BSplinesP_Polar::degree()));

                        const ddc::DiscreteVector<RCellDim> n_overlap_r(
                                n_overlap_cells - first_overlap_element_r.uid());
                        const ddc::DiscreteVector<PCellDim> n_overlap_p(
                                BSplinesP_Polar::degree() + 1);

                        const ddc::DiscreteDomain<RCellDim>
                                r_cells(first_overlap_element_r, n_overlap_r);
                        const ddc::DiscreteDomain<PCellDim>
                                p_cells(first_overlap_element_p, n_overlap_p);
                        const ddc::DiscreteDomain<RCellDim, PCellDim>
                                non_zero_cells(r_cells, p_cells);

                        if (n_overlap_r > 0) {
                            double element = 0.0;

                            ddc::for_each(non_zero_cells, [&](CellIndex const cell_idx) {
                                const int cell_idx_r(ddc::select<RCellDim>(cell_idx).uid());
                                const int cell_idx_p(pmod(ddc::select<PCellDim>(cell_idx).uid()));

                                const QuadratureDomainRP cell_quad_points(
                                        get_quadrature_points_in_cell(cell_idx_r, cell_idx_p));
                                // Find the column where the non-zero data is stored
                                ddc::DiscreteElement<RBasisSubset> ib_trial_r(
                                        r_idx_trial.uid() - cell_idx_r);
                                ddc::DiscreteElement<PBasisSubset> ib_trial_p(
                                        pmod(p_idx_trial.uid() - cell_idx_p));
                                // Calculate the weak integral
                                element += ddc::transform_reduce(
                                        cell_quad_points,
                                        0.0,
                                        ddc::reducer::sum<double>(),
                                        [&](QuadratureIndexRP const quad_idx) {
                                            QuadratureIndexR const ir
                                                    = ddc::select<QDimRMesh>(quad_idx);
                                            QuadratureIndexP const ip
                                                    = ddc::select<QDimPMesh>(quad_idx);
                                            return weak_integral_element<Mapping>(
                                                    ir,
                                                    ip,
                                                    singular_basis_vals_and_derivs(
                                                            idx_test,
                                                            ir,
                                                            ip),
                                                    r_basis_vals_and_derivs(ib_trial_r, ir),
                                                    p_basis_vals_and_derivs(ib_trial_p, ip),
                                                    coeff_alpha,
                                                    coeff_beta,
                                                    spline_evaluator,
                                                    mapping);
                                        });
                            });
                            matrix_elements[matrix_idx++]
                                    = MatrixElement(idx_test.uid(), polar_idx_trial.uid(), element);
                            matrix_elements[matrix_idx++]
                                    = MatrixElement(polar_idx_trial.uid(), idx_test.uid(), element);
                        }
                    });
        });
        assert(matrix_idx == n_elements_singular + n_elements_overlap);

        // Calculate the matrix elements following a stencil
        ddc::for_each(fem_non_singular_domain, [&](IndexPolarBspl const polar_idx_test) {
            const IDimBSpline2D_Polar idx_test(PolarBSplinesRP::get_2d_index(polar_idx_test));
            const std::size_t r_idx_test(
                    ddc::select<PolarBSplinesRP::BSplinesR_tag>(idx_test).uid());
            const std::size_t p_idx_test(
                    ddc::select<PolarBSplinesRP::BSplinesP_tag>(idx_test).uid());

            // Calculate the index of the elements that are already filled
            BSDomainP_Polar remaining_p(
                    ddc::DiscreteElement<BSplinesP_Polar> {p_idx_test},
                    ddc::DiscreteVector<BSplinesP_Polar> {BSplinesP_Polar::degree() + 1});
            ddc::for_each(remaining_p, [&](auto const p_idx_trial) {
                IDimBSpline2D_Polar
                        idx_trial(ddc::DiscreteElement<BSplinesR_Polar>(r_idx_test), p_idx_trial);
                IndexPolarBspl polar_idx_trial(PolarBSplinesRP::get_polar_index<PolarBSplinesRP>(
                        IDimBSpline2D_Polar(r_idx_test, pmod(p_idx_trial.uid()))));
                double element = get_matrix_stencil_element(
                        idx_test,
                        idx_trial,
                        coeff_alpha,
                        coeff_beta,
                        spline_evaluator,
                        mapping);

                if (polar_idx_test.uid() == polar_idx_trial.uid()) {
                    matrix_elements[matrix_idx++]
                            = MatrixElement(polar_idx_test.uid(), polar_idx_trial.uid(), element);
                } else {
                    matrix_elements[matrix_idx++]
                            = MatrixElement(polar_idx_test.uid(), polar_idx_trial.uid(), element);
                    matrix_elements[matrix_idx++]
                            = MatrixElement(polar_idx_trial.uid(), polar_idx_test.uid(), element);
                }
            });
            BSDomainR remaining_r(
                    ddc::select<BSplinesR_Polar>(idx_test) + 1,
                    ddc::DiscreteVector<BSplinesR_Polar> {
                            min(BSplinesR_Polar::degree(),
                                ddc::discrete_space<BSplinesR_Polar>().nbasis() - 2 - r_idx_test)});
            BSDomainP_Polar relevant_p(
                    ddc::DiscreteElement<BSplinesP_Polar> {
                            p_idx_test + ddc::discrete_space<BSplinesP_Polar>().nbasis()
                            - BSplinesP_Polar::degree()},
                    ddc::DiscreteVector<BSplinesP_Polar> {2 * BSplinesP_Polar::degree() + 1});

            BSDomainRP trial_domain(remaining_r, relevant_p);

            ddc::for_each(trial_domain, [&](IDimBSpline2D_Polar const idx_trial) {
                const int r_idx_trial(ddc::select<PolarBSplinesRP::BSplinesR_tag>(idx_trial).uid());
                const int p_idx_trial(ddc::select<PolarBSplinesRP::BSplinesP_tag>(idx_trial).uid());
                IndexPolarBspl polar_idx_trial(PolarBSplinesRP::get_polar_index<PolarBSplinesRP>(
                        IDimBSpline2D_Polar(r_idx_trial, pmod(p_idx_trial))));
                double element = get_matrix_stencil_element(
                        idx_test,
                        idx_trial,
                        coeff_alpha,
                        coeff_beta,
                        spline_evaluator,
                        mapping);
                if (polar_idx_test.uid() == polar_idx_trial.uid()) {
                    matrix_elements[matrix_idx++]
                            = MatrixElement(polar_idx_test.uid(), polar_idx_trial.uid(), element);
                } else {
                    matrix_elements[matrix_idx++]
                            = MatrixElement(polar_idx_test.uid(), polar_idx_trial.uid(), element);
                    matrix_elements[matrix_idx++]
                            = MatrixElement(polar_idx_trial.uid(), polar_idx_test.uid(), element);
                }
            });
        });
        matrix.setFromTriplets(matrix_elements.begin(), matrix_elements.end());
        assert(matrix_idx == n_elements_singular + n_elements_overlap + n_elements_stencil);
        m_matrix.compute(matrix);
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
        Eigen::VectorXd b(
                ddc::discrete_space<PolarBSplinesRP>().nbasis()
                - ddc::discrete_space<BSplinesP_Polar>().nbasis());

        // Fill b
        ddc::for_each(
                PolarBSplinesRP::singular_domain<PolarBSplinesRP>(),
                [&](IndexPolarBspl const idx) {
                    b(idx.uid()) = ddc::transform_reduce(
                            quadrature_domain_singular,
                            0.0,
                            ddc::reducer::sum<double>(),
                            [&](QuadratureIndexRP const quad_idx) {
                                QuadratureIndexR const ir = ddc::select<QDimRMesh>(quad_idx);
                                QuadratureIndexP const ip = ddc::select<QDimPMesh>(quad_idx);
                                CoordRP coord(ddc::coordinate(quad_idx));
                                return rhs(coord)
                                       * singular_basis_vals_and_derivs(idx, ir, ip).value
                                       * int_volume(ir, ip);
                            });
                });
        const std::size_t ncells_r = ddc::discrete_space<BSplinesR_Polar>().ncells();
        ddc::for_each(fem_non_singular_domain, [&](IndexPolarBspl const idx) {
            const IDimBSpline2D_Polar idx_2d(PolarBSplinesRP::get_2d_index(idx));
            const std::size_t r_idx(ddc::select<PolarBSplinesRP::BSplinesR_tag>(idx_2d).uid());
            const std::size_t p_idx(ddc::select<PolarBSplinesRP::BSplinesP_tag>(idx_2d).uid());

            // Find the cells on which the bspline is non-zero
            int first_cell_r(r_idx - BSplinesR_Polar::degree());
            int first_cell_p(p_idx - BSplinesP_Polar::degree());
            std::size_t last_cell_r(r_idx + 1);
            if (first_cell_r < 0)
                first_cell_r = 0;
            if (last_cell_r > ncells_r)
                last_cell_r = ncells_r;
            ddc::DiscreteVector<RCellDim> const r_length(last_cell_r - first_cell_r);
            ddc::DiscreteVector<PCellDim> const p_length(BSplinesP_Polar::degree() + 1);


            ddc::DiscreteElement<RCellDim> const start_r(first_cell_r);
            ddc::DiscreteElement<PCellDim> const start_p(pmod(first_cell_p));
            const ddc::DiscreteDomain<RCellDim> r_cells(start_r, r_length);
            const ddc::DiscreteDomain<PCellDim> p_cells(start_p, p_length);
            const ddc::DiscreteDomain<RCellDim, PCellDim> non_zero_cells(r_cells, p_cells);
            assert(r_length * p_length > 0);
            double element = 0.0;
            ddc::for_each(non_zero_cells, [&](CellIndex const cell_idx) {
                const int cell_idx_r(ddc::select<RCellDim>(cell_idx).uid());
                const int cell_idx_p(pmod(ddc::select<PCellDim>(cell_idx).uid()));

                const QuadratureDomainRP cell_quad_points(
                        get_quadrature_points_in_cell(cell_idx_r, cell_idx_p));

                // Find the column where the non-zero data is stored
                ddc::DiscreteElement<RBasisSubset> ib_r(r_idx - cell_idx_r);
                ddc::DiscreteElement<PBasisSubset> ib_p(pmod(p_idx - cell_idx_p));

                // Calculate the weak integral
                element += ddc::transform_reduce(
                        cell_quad_points,
                        0.0,
                        ddc::reducer::sum<double>(),
                        [&](QuadratureIndexRP const quad_idx) {
                            QuadratureIndexR const ir = ddc::select<QDimRMesh>(quad_idx);
                            QuadratureIndexP const ip = ddc::select<QDimPMesh>(quad_idx);
                            CoordRP coord(ddc::coordinate(quad_idx));
                            double rb = r_basis_vals_and_derivs(ib_r, ir).value;
                            double pb = p_basis_vals_and_derivs(ib_p, ip).value;
                            return rhs(coord) * rb * pb * int_volume(ir, ip);
                        });
            });
            b(idx.uid()) = element;
        });

        // Solve the matrix equation
        Eigen::VectorXd x = m_matrix.solve(b);

        BSDomainRP dirichlet_boundary_domain(
                radial_bsplines.take_last(ddc::DiscreteVector<BSplinesR_Polar> {1}),
                polar_bsplines);
        BSDomainP_Polar polar_domain(ddc::discrete_space<BSplinesP_Polar>().full_domain());


        // Fill the spline
        ddc::for_each(
                PolarBSplinesRP::singular_domain<PolarBSplinesRP>(),
                [&](IndexPolarBspl const idx) { spline.singular_spline_coef(idx) = x(idx.uid()); });
        ddc::for_each(fem_non_singular_domain, [&](IndexPolarBspl const idx) {
            const IDimBSpline2D_Polar idx_2d(PolarBSplinesRP::get_2d_index(idx));
            spline.spline_coef(idx_2d) = x(idx.uid());
        });
        ddc::for_each(dirichlet_boundary_domain, [&](IDimBSpline2D_Polar const idx) {
            spline.spline_coef(idx) = 0.0;
        });

        // Copy the periodic elements
        BSDomainRP copy_domain(
                radial_bsplines,
                polar_domain.remove_first(ddc::DiscreteVector<BSplinesP_Polar>(
                        ddc::discrete_space<BSplinesP_Polar>().nbasis())));
        ddc::for_each(copy_domain, [&](IDimBSpline2D_Polar const idx_2d) {
            spline.spline_coef(
                    ddc::select<BSplinesR_Polar>(idx_2d),
                    ddc::select<BSplinesP_Polar>(idx_2d))
                    = spline.spline_coef(
                            ddc::select<BSplinesR_Polar>(idx_2d),
                            ddc::select<BSplinesP_Polar>(idx_2d)
                                    - ddc::discrete_space<BSplinesP_Polar>().nbasis());
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
     *      A Chunk of coordinates where we want to compute the solution.
     * @param[out] result
     *      The values of the solution @f$\phi@f$ on the given coords_eval.
     */
    template <class RHSFunction>
    void operator()(RHSFunction const& rhs, ViewRP<CoordRP> const coords_eval, DSpanRP result) const
    {
        BSDomainP_Polar polar_domain(ddc::discrete_space<BSplinesP_Polar>().full_domain());
        SplinePolar
                spline(PolarBSplinesRP::singular_domain<PolarBSplinesRP>(),
                       BSDomainRP(radial_bsplines, polar_domain));

        (*this)(rhs, spline);
        m_polar_spline_evaluator(result, coords_eval, spline);
    }


private:
    static QuadratureDomainRP get_quadrature_points_in_cell(int cell_idx_r, int cell_idx_p)
    {
        const QuadratureIndexR first_quad_point_r(cell_idx_r * n_gauss_legendre_r);
        const QuadratureIndexP first_quad_point_p(cell_idx_p * n_gauss_legendre_p);
        constexpr QuadratureLengthR n_GL_r(n_gauss_legendre_r);
        constexpr QuadratureLengthP n_GL_p(n_gauss_legendre_p);
        const QuadratureDomainR quad_points_r(first_quad_point_r, n_GL_r);
        const QuadratureDomainP quad_points_p(first_quad_point_p, n_GL_p);
        return QuadratureDomainRP(quad_points_r, quad_points_p);
    }

    template <class Mapping>
    double weak_integral_element(
            QuadratureIndexR ir,
            QuadratureIndexP ip,
            EvalDeriv2DType const& test_bspline_val_and_deriv,
            EvalDeriv2DType const& trial_bspline_val_and_deriv,
            Spline2DView coeff_alpha,
            Spline2DView coeff_beta,
            SplineRPEvaluatorNullBound const& evaluator,
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
            QuadratureIndexR ir,
            QuadratureIndexP ip,
            EvalDeriv2DType const& test_bspline_val_and_deriv,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_r,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_p,
            Spline2DView coeff_alpha,
            Spline2DView coeff_beta,
            SplineRPEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                ir,
                ip,
                test_bspline_val_and_deriv,
                trial_bspline_val_and_deriv_r,
                test_bspline_val_and_deriv,
                trial_bspline_val_and_deriv_p,
                coeff_alpha,
                coeff_beta,
                evaluator,
                mapping);
    }

    template <class Mapping>
    double weak_integral_element(
            QuadratureIndexR ir,
            QuadratureIndexP ip,
            EvalDeriv1DType const& test_bspline_val_and_deriv_r,
            EvalDeriv2DType const& trial_bspline_val_and_deriv,
            EvalDeriv1DType const& test_bspline_val_and_deriv_p,
            Spline2DView coeff_alpha,
            Spline2DView coeff_beta,
            SplineRPEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                ir,
                ip,
                test_bspline_val_and_deriv_r,
                trial_bspline_val_and_deriv,
                test_bspline_val_and_deriv_p,
                trial_bspline_val_and_deriv,
                coeff_alpha,
                coeff_beta,
                evaluator,
                mapping);
    }

    template <class Mapping>
    double weak_integral_element(
            QuadratureIndexR ir,
            QuadratureIndexP ip,
            EvalDeriv1DType const& test_bspline_val_and_deriv_r,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_r,
            EvalDeriv1DType const& test_bspline_val_and_deriv_p,
            EvalDeriv1DType const& trial_bspline_val_and_deriv_p,
            Spline2DView coeff_alpha,
            Spline2DView coeff_beta,
            SplineRPEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        return templated_weak_integral_element(
                ir,
                ip,
                test_bspline_val_and_deriv_r,
                trial_bspline_val_and_deriv_r,
                test_bspline_val_and_deriv_p,
                trial_bspline_val_and_deriv_p,
                coeff_alpha,
                coeff_beta,
                evaluator,
                mapping);
    }

    inline void get_value_and_gradient(
            double& value,
            std::array<double, 2>& gradient,
            EvalDeriv1DType const& r_basis,
            EvalDeriv1DType const& p_basis) const
    {
        value = r_basis.value * p_basis.value;
        gradient = {r_basis.derivative * p_basis.value, r_basis.value * p_basis.derivative};
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
            QuadratureIndexR ir,
            QuadratureIndexP ip,
            TestValDerivType const& test_bspline_val_and_deriv,
            TrialValDerivType const& trial_bspline_val_and_deriv,
            TestValDerivType const& test_bspline_val_and_deriv_p,
            TrialValDerivType const& trial_bspline_val_and_deriv_p,
            Spline2DView coeff_alpha,
            Spline2DView coeff_beta,
            SplineRPEvaluatorNullBound const& spline_evaluator,
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
        CoordRP coord(ddc::coordinate(ir), ddc::coordinate(ip));
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
                test_bspline_val_and_deriv_p);
        get_value_and_gradient(
                basis_val_trial_space,
                basis_gradient_trial_space,
                trial_bspline_val_and_deriv,
                trial_bspline_val_and_deriv_p);

        // Assemble the weak integral element
        return int_volume(ir, ip)
               * (alpha
                          * dot_product(
                                  basis_gradient_test_space,
                                  mapping.to_covariant(basis_gradient_trial_space, coord))
                  + beta * basis_val_test_space * basis_val_trial_space);
    }

    /**
     * @brief Computes the matrix element corresponding to two tensor product splines
     *        with index idx_test and idx_trial
     */
    template <class Mapping>
    double get_matrix_stencil_element(
            IDimBSpline2D_Polar idx_test,
            IDimBSpline2D_Polar idx_trial,
            Spline2DView coeff_alpha,
            Spline2DView coeff_beta,
            SplineRPEvaluatorNullBound const& evaluator,
            Mapping const& mapping)
    {
        // 0 <= r_idx_test < 8
        // 0 <= r_idx_trial < 8
        // r_idx_test < r_idx_trial
        const int r_idx_test(ddc::select<BSplinesR_Polar>(idx_test).uid());
        const int r_idx_trial(ddc::select<BSplinesR_Polar>(idx_trial).uid());
        // 0 <= p_idx_test < 8
        // 0 <= p_idx_trial < 8
        int p_idx_test(pmod(ddc::select<BSplinesP_Polar>(idx_test).uid()));
        int p_idx_trial(pmod(ddc::select<BSplinesP_Polar>(idx_trial).uid()));

        const std::size_t ncells_r = ddc::discrete_space<BSplinesR_Polar>().ncells();
        const std::size_t nbasis_p = ddc::discrete_space<BSplinesP_Polar>().nbasis();

        // 0<= r_offset <= degree_r
        // -degree_p <= p_offset <= degree_p
        const int r_offset = r_idx_trial - r_idx_test;
        int p_offset = pmod(p_idx_trial - p_idx_test);
        if (p_offset >= int(nbasis_p - BSplinesP_Polar::degree())) {
            p_offset -= nbasis_p;
        }
        assert(r_offset >= 0);
        assert(r_offset <= int(BSplinesR_Polar::degree()));
        assert(p_offset >= -int(BSplinesP_Polar::degree()));
        assert(p_offset <= int(BSplinesP_Polar::degree()));

        // Find the domain covering the cells where both the test and trial functions are non-zero
        int n_overlap_stencil_r(BSplinesR_Polar::degree() + 1 - r_offset);
        int first_overlap_r(r_idx_trial - BSplinesR_Polar::degree());

        int first_overlap_p;
        int n_overlap_stencil_p;
        if (p_offset > 0) {
            n_overlap_stencil_p = BSplinesP_Polar::degree() + 1 - p_offset;
            first_overlap_p = pmod(p_idx_trial - BSplinesP_Polar::degree());
        } else {
            n_overlap_stencil_p = BSplinesP_Polar::degree() + 1 + p_offset;
            first_overlap_p = pmod(p_idx_test - BSplinesP_Polar::degree());
        }

        if (first_overlap_r < 0) {
            const int n_compact = first_overlap_r;
            first_overlap_r = 0;
            n_overlap_stencil_r += n_compact;
        }

        const int n_to_edge_r(ncells_r - first_overlap_r);

        const ddc::DiscreteVector<RCellDim> n_overlap_r(min(n_overlap_stencil_r, n_to_edge_r));
        const ddc::DiscreteVector<PCellDim> n_overlap_p(n_overlap_stencil_p);

        const ddc::DiscreteElement<RCellDim> first_overlap_element_r(first_overlap_r);
        const ddc::DiscreteElement<PCellDim> first_overlap_element_p(first_overlap_p);

        const ddc::DiscreteDomain<RCellDim> r_cells(first_overlap_element_r, n_overlap_r);
        const ddc::DiscreteDomain<PCellDim> p_cells(first_overlap_element_p, n_overlap_p);
        const ddc::DiscreteDomain<RCellDim, PCellDim> non_zero_cells(r_cells, p_cells);

        assert(n_overlap_r * n_overlap_p > 0);
        return ddc::transform_reduce(
                non_zero_cells,
                0.0,
                ddc::reducer::sum<double>(),
                [&](CellIndex const cell_idx) {
                    const int cell_idx_r(ddc::select<RCellDim>(cell_idx).uid());
                    const int cell_idx_p(pmod(ddc::select<PCellDim>(cell_idx).uid()));

                    const QuadratureDomainRP cell_quad_points(
                            get_quadrature_points_in_cell(cell_idx_r, cell_idx_p));

                    int ib_test_p_idx = p_idx_test - cell_idx_p;
                    int ib_trial_p_idx = p_idx_trial - cell_idx_p;

                    // Find the column where the non-zero data is stored
                    ddc::DiscreteElement<RBasisSubset> ib_test_r(r_idx_test - cell_idx_r);
                    ddc::DiscreteElement<PBasisSubset> ib_test_p(pmod(ib_test_p_idx));
                    ddc::DiscreteElement<RBasisSubset> ib_trial_r(r_idx_trial - cell_idx_r);
                    ddc::DiscreteElement<PBasisSubset> ib_trial_p(pmod(ib_trial_p_idx));

                    assert(ib_test_r.uid() < BSplinesR_Polar::degree() + 1);
                    assert(ib_test_p.uid() < BSplinesP_Polar::degree() + 1);
                    assert(ib_trial_r.uid() < BSplinesR_Polar::degree() + 1);
                    assert(ib_trial_p.uid() < BSplinesP_Polar::degree() + 1);

                    // Calculate the weak integral
                    return ddc::transform_reduce(
                            cell_quad_points,
                            0.0,
                            ddc::reducer::sum<double>(),
                            [&](QuadratureIndexRP const quad_idx) {
                                QuadratureIndexR const r_idx = ddc::select<QDimRMesh>(quad_idx);
                                QuadratureIndexP const p_idx = ddc::select<QDimPMesh>(quad_idx);
                                return weak_integral_element(
                                        r_idx,
                                        p_idx,
                                        r_basis_vals_and_derivs(ib_test_r, r_idx),
                                        r_basis_vals_and_derivs(ib_trial_r, r_idx),
                                        p_basis_vals_and_derivs(ib_test_p, p_idx),
                                        p_basis_vals_and_derivs(ib_trial_p, p_idx),
                                        coeff_alpha,
                                        coeff_beta,
                                        evaluator,
                                        mapping);
                            });
                });
    }

    static int pmod(int idx_p)
    {
        int ncells_p = ddc::discrete_space<BSplinesP_Polar>().ncells();
        while (idx_p < 0)
            idx_p += ncells_p;
        while (idx_p >= ncells_p)
            idx_p -= ncells_p;
        return idx_p;
    }
};
