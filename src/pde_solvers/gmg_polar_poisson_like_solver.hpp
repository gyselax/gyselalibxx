// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <vector>

#include <GMGPolar/gmgpolar.h>

#include "ddc_alias_inline_functions.hpp"
#include "i_interpolation_builder.hpp"
#include "ipolar_poisson_like_solver.hpp"

namespace GMGPolarTools {

/**
 * @brief Wraps a gyselalibxx coordinate mapping to satisfy the GMGPolar DomainGeometry concept.
 * @tparam ToPhysicalMapping A mapping from (r, theta) curvilinear coordinates to (x, y) Cartesian.
 */
template <class ToPhysicalMapping>
class MappingToDomainGeometry
{
    using R = typename CoordWithOPoint<typename ToPhysicalMapping::CoordArg>::curvilinear_tag_r;
    using Theta = typename CoordWithOPoint<typename ToPhysicalMapping::CoordArg>::curvilinear_tag_theta;

    using X = ddc::
            type_seq_element_t<0, ddc::to_type_seq_t<typename ToPhysicalMapping::CoordResult>>;
    using Y = ddc::
            type_seq_element_t<1, ddc::to_type_seq_t<typename ToPhysicalMapping::CoordResult>>;

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
 * @brief Wraps gyselalibxx interpolation-represented coefficients to satisfy the GMGPolar
 *        DensityProfileCoefficients concept.
 * @tparam EvaluatorType A 2D evaluator for the representation described by IdxRangeCoeff.
 */
template <class EvaluatorType, class IdxRangeCoeff, class CoordRTheta>
class PolarPoissonLikeCoefficients
{
    using DConstCoeffRTheta = DConstField<IdxRangeCoeff>;

private:
    EvaluatorType m_evaluator;
    DConstCoeffRTheta m_coeff_alpha;
    DConstCoeffRTheta m_coeff_beta;

public:
    /// Build the class instance
    PolarPoissonLikeCoefficients(
            EvaluatorType evaluator,
            DConstCoeffRTheta coeff_alpha,
            DConstCoeffRTheta coeff_beta)
        : m_evaluator(evaluator)
        , m_coeff_alpha(coeff_alpha)
        , m_coeff_beta(coeff_beta)
    {
    }

    /// The coefficient alpha in the Poisson-like equation
    KOKKOS_INLINE_FUNCTION double alpha(const double& r, const double& theta) const
    {
        return m_evaluator(CoordRTheta(r, theta), m_coeff_alpha);
    }
    /// The coefficient beta in the Poisson-like equation
    KOKKOS_INLINE_FUNCTION double beta(const double& r, const double& theta) const
    {
        return m_evaluator(CoordRTheta(r, theta), m_coeff_beta);
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
 * @tparam InterpolatorType     2D interpolator for (GridR × GridTheta).
 */
template <class ToPhysicalMapping, class GridR, class GridTheta, class InterpolatorType>
class GMGPolarPoissonLikeSolver
    : public IPolarPoissonLikeSolver<IdxRange<GridR, GridTheta>, IdxRange<GridR, GridTheta>>
{
    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;

    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxStepRTheta = IdxStep<GridR, GridTheta>;

    using BuilderType = typename InterpolatorType::BuilderType;
    using EvaluatorType = typename InterpolatorType::EvaluatorType;

    using IdxRangeCoeff = typename InterpolationBuilderTraits<BuilderType>::coeff_idx_range_type;
    using CoeffRThetaMem = DFieldMem<IdxRangeCoeff>;

    using DomainGeometry = GMGPolarTools::MappingToDomainGeometry<ToPhysicalMapping>;
    using DensityCoeffs = GMGPolarTools::
            PolarPoissonLikeCoefficients<EvaluatorType, IdxRangeCoeff, Coord<R, Theta>>;

private:
    DomainGeometry const m_domain_geom;
    BuilderType const& m_builder;
    EvaluatorType const& m_evaluator;
    ExtrapolationType const m_extrapolation_rule;
    CoeffRThetaMem m_coeff_alpha;
    CoeffRThetaMem m_coeff_beta;
    DensityCoeffs const m_density_coeffs;
    int m_max_iterations;
    double m_absTol;
    double m_relTol;

    DFieldMem<IdxRangeR> m_r_coords;
    DFieldMem<IdxRangeTheta> m_theta_coords;
    gmgpolar::PolarGrid m_polar_grid;
    std::unique_ptr<gmgpolar::GMGPolar<DomainGeometry, DensityCoeffs>> m_solver;


public:
    /**
     * @brief Construct a GMGPolarPoissonLikeSolver.
     *
     * @param[in] to_physical The mapping from the logical to the physical domain.
     * @param[in] interpolator An interpolator to construct and evaluate the coefficients of the interpolation.
     * @param[in] extrapolation_rule A parameter to pass extrapolation rule to GMGPolar, default ExtrapolationType::NONE.
     * @param[in] max_iterations The maximum number of iterations that the solver should carry out.
     * @param[in] absTol The absolute tolerance for the convergence of the solver.
     * @param[in] relTol The relative tolerance for the convergence of the solver.
     */
    GMGPolarPoissonLikeSolver(
            ToPhysicalMapping to_physical,
            InterpolatorType const& interpolator,
            ExtrapolationType const extrapolation_rule = ExtrapolationType::NONE,
            std::optional<int> max_iterations = std::nullopt,
            std::optional<double> absTol = std::nullopt,
            std::optional<double> relTol = std::nullopt)
        : m_domain_geom(to_physical)
        , m_builder(interpolator.get_builder())
        , m_evaluator(interpolator.get_evaluator())
        , m_extrapolation_rule(extrapolation_rule)
        , m_coeff_alpha(get_spline_idx_range(m_builder))
        , m_coeff_beta(get_spline_idx_range(m_builder))
        , m_density_coeffs(
                  m_evaluator,
                  get_const_field(m_coeff_alpha),
                  get_const_field(m_coeff_beta))
        , m_max_iterations(max_iterations.value_or(100))
        , m_absTol(absTol.value_or(1e-10))
        , m_relTol(relTol.value_or(1e-6))
    {
        IdxRangeRTheta idx_range(m_builder.interpolation_domain());
        IdxRangeR idx_range_r(idx_range);
        IdxRangeTheta idx_range_theta(idx_range);
        IdxRangeTheta idx_range_theta_with_poloidal_point(
                idx_range_theta.front(),
                idx_range_theta.extents() + 1);
        m_r_coords = DFieldMem<IdxRangeR>(idx_range_r);
        m_theta_coords = DFieldMem<IdxRangeTheta>(idx_range_theta_with_poloidal_point);
        ddcHelper::dump_coordinates(Kokkos::DefaultExecutionSpace(), get_field(m_r_coords));
        ddcHelper::dump_coordinates(Kokkos::DefaultExecutionSpace(), get_field(m_theta_coords));

        // --- Create a cartesian grid representation of the polar grid for GMGPolar --- //
        m_polar_grid = gmgpolar::PolarGrid(
                m_r_coords.allocation_kokkos_view(),
                m_theta_coords.allocation_kokkos_view());
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

        // --- Create GMGPolar solver for the selected geometry and coefficients --- //
        m_solver = std::make_unique<gmgpolar::GMGPolar<
                DomainGeometry,
                DensityCoeffs>>(m_polar_grid, m_domain_geom, m_density_coeffs);

        // ------------------//
        // Solver parameters //
        // ------------------//

        // --- General solver output and visualisation settings --- //
        m_solver->verbose(0); // Enable/disable verbose output
        m_solver->paraview(false); // Enable/disable ParaView output

        // --- Numerical method setup --- //
        // Are boundary conditions provided on the interior. False = Use Across-the-origin discretisation
        m_solver->DirBC_Interior(false);
        // Stencil distribution strategy: Take, Give
        m_solver->stencilDistributionMethod(StencilDistributionMethod::TAKE);
        // Cache density profile coefficients: alpha, beta
        m_solver->cacheDensityProfileCoefficients(true);
        // Cache domain geometry data: arr, att, art, detDF
        m_solver->cacheDomainGeometry(true);

        // --- Multigrid settings --- //
        m_solver->extrapolation(m_extrapolation_rule); // Select extrapolation
        m_solver->maxLevels(-1); // Max multigrid levels (-1 = use deepest possible)
        m_solver->preSmoothingSteps(1); // Smoothing before coarse-grid correction
        m_solver->postSmoothingSteps(1); // Smoothing after coarse-grid correction
        m_solver->multigridCycle(MultigridCycleType::V_CYCLE); // Multigrid cycle type
        m_solver->FMG(true); // Full Multigrid mode on/off
        m_solver->FMG_iterations(2); // FMG iteration count
        m_solver->FMG_cycle(MultigridCycleType::F_CYCLE); // FMG cycle type

        // --- Preconditioned Conjugate Gradient settings --- //
        m_solver->PCG(false); // Preconditioned Conjugate Gradient mode on/off
        m_solver->PCG_FMG(true); // Use FMG as preconditioner for PCG
        m_solver->PCG_FMG_iterations(1); // FMG iterations for PCG preconditioner
        m_solver->PCG_FMG_cycle(
                MultigridCycleType::V_CYCLE); // FMG cycle type for PCG preconditioner
        m_solver->PCG_MG_iterations(2); // Multigrid iterations for PCG preconditioner
        m_solver->PCG_MG_cycle(
                MultigridCycleType::V_CYCLE); // Multigrid cycle type for PCG iterations

        // --- Iterative solver controls --- //
        m_solver->maxIterations(m_max_iterations); // Max number of iterations
        m_solver->residualNormType(ResidualNormType::WEIGHTED_EUCLIDEAN); // Residual norm type
        m_solver->absoluteTolerance(m_absTol); // Absolute residual tolerance
        m_solver->relativeTolerance(m_relTol); // Relative residual tolerance

        // --- Finalise solver setup --- //
        m_solver->setup();
    }

    /**
     * @brief Solve the Poisson-like equation.
     *
     * @param[out] phi The solution @f$\phi@f$ on the grid.
     * @param[in]  rho The right-hand side @f$\rho@f$ on the grid.
     */
    void operator()(DField<IdxRangeRTheta> phi, DConstField<IdxRangeRTheta> rho) const override
    {
        // Source term: maps GMGPolar (i_r, i_theta) indices to rho grid values
        GMGPolarTools::HomogeneousDirichletBoundaryConditions const bcs;
        m_solver->solve(bcs, rho.allocation_kokkos_view());

        // Copy solution back to phi
        Kokkos::View<double*> solution = m_solver->solution();

        gmgpolar::PolarGrid const& polar_grid = m_polar_grid;
        IdxRangeRTheta idx_range(get_idx_range(phi));

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
