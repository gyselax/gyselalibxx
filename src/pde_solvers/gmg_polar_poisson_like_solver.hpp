// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <vector>

#include <GMGPolar/gmgpolar.h>

#include "ddc_alias_inline_functions.hpp"
#include "ipolar_poisson_like_solver.hpp"

namespace GMGPolarTools {

/**
 * @brief Wraps a gyselalibxx coordinate mapping to satisfy the GMGPolar DomainGeometry concept.
 * @tparam ToPhysicalMapping A mapping from (r, theta) curvilinear coordinates to (x, y) Cartesian.
 */
template <class ToPhysicalMapping>
class MappingToDomainGeometry
{
    using R = typename ToPhysicalMapping::curvilinear_tag_r;
    using Theta = typename ToPhysicalMapping::curvilinear_tag_theta;

    using X = typename ToPhysicalMapping::cartesian_tag_x;
    using Y = typename ToPhysicalMapping::cartesian_tag_y;

    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

private:
    ToPhysicalMapping m_to_physical;

public:
    /// Construct the wrapper class
    explicit MappingToDomainGeometry(ToPhysicalMapping to_physical) : m_to_physical(to_physical) {}

    /// X(r, theta)
    KOKKOS_INLINE_FUNCTION double Fx(const double& r, const double& theta) const
    {
        return Coord<X>(m_to_physical(Coord<R, Theta>(r, theta)));
    }
    /// Y(r, theta)
    KOKKOS_INLINE_FUNCTION double Fy(const double& r, const double& theta) const
    {
        return Coord<Y>(m_to_physical(Coord<R, Theta>(r, theta)));
    }
    /// d/dr X(r, theta)
    KOKKOS_INLINE_FUNCTION double dFx_dr(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<X, R_cov>(Coord<R, Theta>(r, theta));
    }
    /// d/dr Y(r, theta)
    KOKKOS_INLINE_FUNCTION double dFy_dr(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<Y, R_cov>(Coord<R, Theta>(r, theta));
    }
    /// d/(d theta) X(r, theta)
    KOKKOS_INLINE_FUNCTION double dFx_dtheta(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<X, Theta_cov>(Coord<R, Theta>(r, theta));
    }
    /// d/(d theta) Y(r, theta)
    KOKKOS_INLINE_FUNCTION double dFy_dtheta(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<Y, Theta_cov>(Coord<R, Theta>(r, theta));
    }
};

/**
 * @brief Homogeneous Dirichlet boundary conditions satisfying the GMGPolar BoundaryConditions concept.
 */
class HomogeneousDirichletBoundaryConditions
{
public:
    /// The value of the solution on the boundary
    static KOKKOS_INLINE_FUNCTION double u_D(const double& r, const double& theta)
    {
        return 0.0;
    }
    /// The value of the solution on the inner boundary (at r=rmin). Required for the concept, not needed here.
    static KOKKOS_INLINE_FUNCTION double u_D_Interior(const double& r, const double& theta)
    {
        // Only used if DirBC_Interior = true
        assert(false);
        return 0.0;
    }
};

/**
 * @brief Wraps gyselalibxx spline-represented coefficients to satisfy the GMGPolar
 *        DensityProfileCoefficients concept.
 * @tparam SplineEvaluator A 2D spline evaluator for (BSplinesR, BSplinesTheta).
 * @tparam BSplinesR The radial B-spline type.
 * @tparam BSplinesTheta The poloidal B-spline type.
 */
template <class SplineEvaluator, class BSplinesR, class BSplinesTheta>
class PolarPoissonLikeCoefficients
{
    using R = typename BSplinesR::continuous_dimension_type;
    using Theta = typename BSplinesTheta::continuous_dimension_type;

    using DConstSplineRTheta = DConstField<IdxRange<BSplinesR, BSplinesTheta>>;

private:
    SplineEvaluator m_evaluator;
    DConstSplineRTheta m_coeff_alpha;
    DConstSplineRTheta m_coeff_beta;

public:
    /// Build the class instance
    PolarPoissonLikeCoefficients(
            SplineEvaluator evaluator,
            DConstSplineRTheta coeff_alpha,
            DConstSplineRTheta coeff_beta)
        : m_evaluator(evaluator)
        , m_coeff_alpha(coeff_alpha)
        , m_coeff_beta(coeff_beta)
    {
    }

    /// The coefficient alpha in the Poisson-like equation
    KOKKOS_INLINE_FUNCTION double alpha(const double& r, const double& theta) const
    {
        return m_evaluator(Coord<R, Theta>(r, theta), m_coeff_alpha);
    }
    /// The coefficient beta in the Poisson-like equation
    KOKKOS_INLINE_FUNCTION double beta(const double& r, const double& theta) const
    {
        return m_evaluator(Coord<R, Theta>(r, theta), m_coeff_beta);
    }

    /// Required for the concept, only used in custom mesh generation (refinement_radius); not needed here.
    static double getAlphaJump()
    {
        assert(false);
        return 0.0;
    }
};

} // namespace GMGPolarTools

/**
 * @brief A Poisson-like solver using the GMGPolar multigrid library.
 *
 * Solves -∇·(α∇φ) + βφ = ρ on a polar domain with homogeneous Dirichlet BCs
 * at the outer boundary, using an across-the-origin discretisation at r = 0.
 *
 * @tparam ToPhysicalMapping    Mapping from (r,θ) to (x,y).
 * @tparam GridR                Discrete radial grid.
 * @tparam GridTheta            Discrete poloidal grid.
 * @tparam BSplinesR            Radial B-spline space.
 * @tparam BSplinesTheta        Poloidal B-spline space.
 * @tparam SplineBuilder        2D spline builder for (GridR × GridTheta).
 * @tparam SplineEvaluator      2D spline evaluator for (BSplinesR × BSplinesTheta).
 */
template <
        class ToPhysicalMapping,
        class GridR,
        class GridTheta,
        class BSplinesR,
        class BSplinesTheta,
        class SplineBuilder,
        class SplineEvaluator>
class GMGPolarPoissonLikeSolver
    : public IPolarPoissonLikeSolver<IdxRange<GridR, GridTheta>, IdxRange<GridR, GridTheta>>
{
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxStepRTheta = IdxStep<GridR, GridTheta>;

    using SplineRThetaMem = DFieldMem<IdxRange<BSplinesR, BSplinesTheta>>;

    using DomainGeometry = GMGPolarTools::MappingToDomainGeometry<ToPhysicalMapping>;
    using DensityCoeffs = GMGPolarTools::
            PolarPoissonLikeCoefficients<SplineEvaluator, BSplinesR, BSplinesTheta>;

private:
    DomainGeometry const m_domain_geom;
    SplineBuilder const& m_builder;
    SplineEvaluator const& m_evaluator;
    ExtrapolationType const& m_extrapolation_rule;
    SplineRThetaMem m_coeff_alpha;
    SplineRThetaMem m_coeff_beta;
    DensityCoeffs const m_density_coeffs;
    int m_max_iterations;
    double m_absTol;
    double m_relTol;


public:
    /**
     * @brief Construct a GMGPolarPoissonLikeSolver.
     *
     * @param[in] to_physical The mapping from the logical to the physical domain.
     * @param[in] builder A builder to construct the coefficients of the interpolation.
     * @param[in] evaluator The evaluator for the interpolation.
     * @param[in] extrapolation_rule A parameter to pass extrapolation rule to GMGPolar, default ExtrapolationType::NONE.
     * @param[in] max_iterations The maximum number of iterations that the solver should carry out.
     * @param[in] absTol The absolute tolerance for the convergence of the solver.
     * @param[in] relTol The relative tolerance for the convergence of the solver.
     */
    GMGPolarPoissonLikeSolver(
            ToPhysicalMapping to_physical,
            SplineBuilder const& builder,
            SplineEvaluator const& evaluator,
            ExtrapolationType const& extrapolation_rule = ExtrapolationType::NONE,
            std::optional<int> max_iterations = std::nullopt,
            std::optional<double> absTol = std::nullopt,
            std::optional<double> relTol = std::nullopt)
        : m_domain_geom(to_physical)
        , m_builder(builder)
        , m_evaluator(evaluator)
        , m_extrapolation_rule(extrapolation_rule)
        , m_coeff_alpha(get_spline_idx_range(m_builder))
        , m_coeff_beta(get_spline_idx_range(m_builder))
        , m_density_coeffs(
                  m_evaluator,
                  get_const_field(m_coeff_alpha),
                  get_const_field(m_coeff_beta))
        , m_max_iterations(max_iterations.value_or(100))
        , m_absTol(absTol.value_or(1e-10))
        , m_relTol(relTol.value_or(1e-7))
    {
    }

    /**
     * @brief Rebuild the internal spline representations of α and β from grid values.
     * @param[in] alpha Values of α at the grid interpolation points.
     * @param[in] beta  Values of β at the grid interpolation points.
     */
    void update_coefficients(DConstField<IdxRangeRTheta> alpha, DConstField<IdxRangeRTheta> beta)
            override
    {
        m_builder(get_field(m_coeff_alpha), get_const_field(alpha));
        m_builder(get_field(m_coeff_beta), get_const_field(beta));
    }

    /**
     * @brief Solve the Poisson-like equation.
     *
     * @param[out] phi The solution @f$\phi@f$ on the grid.
     * @param[in]  rho The right-hand side @f$\rho@f$ on the grid.
     */
    void operator()(DField<IdxRangeRTheta> phi, DConstField<IdxRangeRTheta> rho) const override
    {
        IdxRangeRTheta idx_range = get_idx_range(phi);
        IdxRangeR idx_range_r(idx_range);
        IdxRangeTheta idx_range_theta(idx_range);
        IdxRangeTheta idx_range_theta_with_poloidal_point(
                idx_range_theta.front(),
                idx_range_theta.extents() + 1);

        DFieldMem<IdxRangeR> r_coords(idx_range_r);
        DFieldMem<IdxRangeTheta> theta_coords(idx_range_theta_with_poloidal_point);
        ddcHelper::dump_coordinates(Kokkos::DefaultExecutionSpace(), get_field(r_coords));
        ddcHelper::dump_coordinates(Kokkos::DefaultExecutionSpace(), get_field(theta_coords));

        // --- Create a cartesian grid representation of the polar grid for GMGPolar --- //
        gmgpolar::PolarGrid const polar_grid(
                r_coords.allocation_kokkos_view(),
                theta_coords.allocation_kokkos_view());

        // --- Create GMGPolar solver for the selected geometry and coefficients --- //
        gmgpolar::GMGPolar<DomainGeometry, DensityCoeffs>
                solver(polar_grid, m_domain_geom, m_density_coeffs);

        // ------------------//
        // Solver parameters //
        // ------------------//

        // --- General solver output and visualisation settings --- //
        solver.verbose(0); // Enable/disable verbose output
        solver.paraview(false); // Enable/disable ParaView output

        // --- Numerical method setup --- //
        // Are boundary conditions provided on the interior. False = Use Across-the-origin discretisation
        solver.DirBC_Interior(false);
        // Stencil distribution strategy: Take, Give
        solver.stencilDistributionMethod(StencilDistributionMethod::TAKE);
        // Cache density profile coefficients: alpha, beta
        solver.cacheDensityProfileCoefficients(true);
        // Cache domain geometry data: arr, att, art, detDF
        solver.cacheDomainGeometry(true);

        // --- Multigrid settings --- //
        solver.extrapolation(m_extrapolation_rule); // Select extrapolation
        solver.maxLevels(-1); // Max multigrid levels (-1 = use deepest possible)
        solver.preSmoothingSteps(1); // Smoothing before coarse-grid correction
        solver.postSmoothingSteps(1); // Smoothing after coarse-grid correction
        solver.multigridCycle(MultigridCycleType::V_CYCLE); // Multigrid cycle type
        solver.FMG(true); // Full Multigrid mode on/off
        solver.FMG_iterations(2); // FMG iteration count
        solver.FMG_cycle(MultigridCycleType::F_CYCLE); // FMG cycle type

        // --- Preconditioned Conjugate Gradient settings --- //
        solver.PCG(false); // Preconditioned Conjugate Gradient mode on/off
        solver.PCG_FMG(true); // Use FMG as preconditioner for PCG
        solver.PCG_FMG_iterations(1); // FMG iterations for PCG preconditioner
        solver.PCG_FMG_cycle(MultigridCycleType::V_CYCLE); // FMG cycle type for PCG preconditioner
        solver.PCG_MG_iterations(1); // Multigrid iterations for PCG preconditioner
        solver.PCG_MG_cycle(MultigridCycleType::V_CYCLE); // Multigrid cycle type for PCG iterations

        // --- Iterative solver controls --- //
        solver.maxIterations(m_max_iterations); // Max number of iterations
        solver.residualNormType(ResidualNormType::WEIGHTED_EUCLIDEAN); // Residual norm type
        solver.absoluteTolerance(m_absTol); // Absolute residual tolerance
        solver.relativeTolerance(m_relTol); // Relative residual tolerance

        // --- Finalise solver setup --- //
        solver.setup();

        // Source term: maps GMGPolar (i_r, i_theta) indices to rho grid values
        GMGPolarTools::HomogeneousDirichletBoundaryConditions const bcs;
        solver.solve(bcs, rho.allocation_kokkos_view());

        // Copy solution back to phi
        //Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solution = solver.solution();
        Kokkos::View<double*> solution = solver.solution();

        ddc::parallel_for_each(
                idx_range,
                KOKKOS_LAMBDA(IdxRTheta idx) {
                    IdxStepRTheta offset(idx - idx_range.front());
                    int i_r = ddc::select<GridR>(offset);
                    int i_theta = ddc::select<GridTheta>(offset);
                    phi(idx) = solution[polar_grid.index(i_r, i_theta)];
                });
    }
};
