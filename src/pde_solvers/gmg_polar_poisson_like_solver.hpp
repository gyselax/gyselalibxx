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
    double Fx(const double& r, const double& theta) const
    {
        return Coord<X>(m_to_physical(Coord<R, Theta>(r, theta)));
    }
    /// Y(r, theta)
    double Fy(const double& r, const double& theta) const
    {
        return Coord<Y>(m_to_physical(Coord<R, Theta>(r, theta)));
    }
    /// d/dr X(r, theta)
    double dFx_dr(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<X, R_cov>(Coord<R, Theta>(r, theta));
    }
    /// d/dr Y(r, theta)
    double dFy_dr(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<Y, R_cov>(Coord<R, Theta>(r, theta));
    }
    /// d/(d theta) X(r, theta)
    double dFx_dt(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<X, Theta_cov>(Coord<R, Theta>(r, theta));
    }
    /// d/(d theta) Y(r, theta)
    double dFy_dt(const double& r, const double& theta) const
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
    double u_D(const double& r, const double& theta) const
    {
        return 0.0;
    }
    /// The value of the solution on the inner boundary (at r=rmin). Required for the concept, not needed here.
    double u_D_Interior(const double& r, const double& theta) const
    {
        // Only used if DirBC_Interior = true
        assert(false);
        return 0.0;
    }
};

/**
 * @brief Wraps gyselalibxx spline-represented coefficients to satisfy the GMGPolar
 *        DensityProfileCoefficients concept.
 * @tparam SplineEvaluator_host A 2D host-space spline evaluator for (BSplinesR, BSplinesTheta).
 * @tparam BSplinesR The radial B-spline type.
 * @tparam BSplinesTheta The poloidal B-spline type.
 */
template <class SplineEvaluator_host, class BSplinesR, class BSplinesTheta>
class PolarPoissonLikeCoefficients
{
    static_assert(
            std::is_same_v<typename SplineEvaluator_host::memory_space, Kokkos::HostSpace>,
            "SplineEvaluator_host must operate on Kokkos::HostSpace");

    using R = typename BSplinesR::continuous_dimension_type;
    using Theta = typename BSplinesTheta::continuous_dimension_type;

    using DConstSplineRTheta_host
            = DConstField<IdxRange<BSplinesR, BSplinesTheta>, Kokkos::HostSpace>;

private:
    SplineEvaluator_host m_evaluator;
    DConstSplineRTheta_host m_coeff_alpha;
    DConstSplineRTheta_host m_coeff_beta;

public:
    /// Build th class instance
    PolarPoissonLikeCoefficients(
            SplineEvaluator_host evaluator,
            DConstSplineRTheta_host coeff_alpha,
            DConstSplineRTheta_host coeff_beta)
        : m_evaluator(evaluator)
        , m_coeff_alpha(coeff_alpha)
        , m_coeff_beta(coeff_beta)
    {
    }

    /// The coefficient alpha in the Poisson-like equation
    double alpha(const double& r, const double& theta) const
    {
        return m_evaluator(Coord<R, Theta>(r, theta), m_coeff_alpha);
    }
    /// The coefficient beta in the Poisson-like equation
    double beta(const double& r, const double& theta) const
    {
        return m_evaluator(Coord<R, Theta>(r, theta), m_coeff_beta);
    }

    /// Required for the concept, only used in custom mesh generation (refinement_radius); not needed here.
    double getAlphaJump() const
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
 * @tparam SplineBuilder_host   2D host-space spline builder for (GridR × GridTheta).
 * @tparam SplineEvaluator_host 2D host-space spline evaluator for (BSplinesR × BSplinesTheta).
 */
template <
        class ToPhysicalMapping,
        class GridR,
        class GridTheta,
        class BSplinesR,
        class BSplinesTheta,
        class SplineBuilder_host,
        class SplineEvaluator_host>
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

    using SplineRThetaMem_host = DFieldMem<IdxRange<BSplinesR, BSplinesTheta>, Kokkos::HostSpace>;

    using DomainGeometry = GMGPolarTools::MappingToDomainGeometry<ToPhysicalMapping>;
    using DensityCoeffs = GMGPolarTools::
            PolarPoissonLikeCoefficients<SplineEvaluator_host, BSplinesR, BSplinesTheta>;

private:
    DomainGeometry const m_domain_geom;
    SplineBuilder_host const& m_builder;
    SplineEvaluator_host const& m_evaluator;
    SplineRThetaMem_host m_coeff_alpha;
    SplineRThetaMem_host m_coeff_beta;
    DensityCoeffs const m_density_coeffs;


public:
    GMGPolarPoissonLikeSolver(
            ToPhysicalMapping to_physical,
            SplineBuilder_host const& builder,
            SplineEvaluator_host const& evaluator)
        : m_domain_geom(to_physical)
        , m_builder(builder)
        , m_evaluator(evaluator)
        , m_coeff_alpha(get_spline_idx_range(m_builder))
        , m_coeff_beta(get_spline_idx_range(m_builder))
        , m_density_coeffs(
                  m_evaluator,
                  get_const_field(m_coeff_alpha),
                  get_const_field(m_coeff_beta))
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
        auto alpha_host = ddc::create_mirror_view_and_copy(alpha);
        auto beta_host = ddc::create_mirror_view_and_copy(beta);
        m_builder(get_field(m_coeff_alpha), get_const_field(alpha_host));
        m_builder(get_field(m_coeff_beta), get_const_field(beta_host));
    }

    /**
     * @brief Solve the Poisson-like equation.
     *
     * @param[out] phi The solution @f$\phi@f$ on the grid.
     * @param[in]  rho The right-hand side @f$\rho@f$ on the grid.
     */
    void operator()(DField<IdxRangeRTheta> phi, DConstField<IdxRangeRTheta> rho) const override
    {
        // Copy rho to host
        auto rho_host = ddc::create_mirror_view_and_copy(rho);

        IdxRangeRTheta idx_range = get_idx_range(phi);
        IdxRangeR idx_range_r(idx_range);
        IdxRangeTheta idx_range_theta(idx_range);
        IdxRangeTheta idx_range_theta_with_poloidal_point(
                idx_range_theta.front(),
                idx_range_theta.extents() + 1);

        host_t<DFieldMem<IdxRangeR>> r_coords(idx_range_r);
        host_t<DFieldMem<IdxRangeTheta>> theta_coords(idx_range_theta_with_poloidal_point);
        ddcHelper::dump_coordinates(Kokkos::DefaultHostExecutionSpace(), get_field(r_coords));
        ddcHelper::dump_coordinates(Kokkos::DefaultHostExecutionSpace(), get_field(theta_coords));

        gmgpolar::PolarGrid const polar_grid(
                r_coords.allocation_kokkos_view(),
                theta_coords.allocation_kokkos_view());

        gmgpolar::GMGPolar<DomainGeometry, DensityCoeffs>
                solver(polar_grid, m_domain_geom, m_density_coeffs);

        // ----------------------------------------------------------------
        // Solver parameters
        solver.DirBC_Interior(false); // Use across-origin discretisation
        solver.FMG(true);
        solver.FMG_iterations(3);
        solver.FMG_cycle(MultigridCycleType::F_CYCLE);
        solver.extrapolation(ExtrapolationType::IMPLICIT_EXTRAPOLATION);
        solver.maxLevels(7);
        solver.preSmoothingSteps(1);
        solver.postSmoothingSteps(1);
        solver.multigridCycle(MultigridCycleType::F_CYCLE);
        solver.maxIterations(150);
        solver.residualNormType(ResidualNormType::EUCLIDEAN);
        solver.absoluteTolerance(1e-50);
        solver.relativeTolerance(1e-6);
        // ----------------------------------------------------------------

        solver.setup();

        // Source term: maps GMGPolar (i_r, i_theta) indices to rho grid values
        GMGPolarTools::HomogeneousDirichletBoundaryConditions const bcs;
        solver.solve(bcs, rho_host.allocation_kokkos_view());

        // Copy solution back to phi
        auto phi_host = ddc::create_mirror_view(Kokkos::DefaultHostExecutionSpace(), phi);
        //Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solution = solver.solution();
        Kokkos::View<double*> solution = solver.solution();

        ddc::host_for_each(idx_range, [&](IdxRTheta idx) {
            IdxStepRTheta offset(idx - idx_range.front());
            int i_r = ddc::select<GridR>(offset);
            int i_theta = ddc::select<GridTheta>(offset);
            phi_host(idx) = solution[polar_grid.index(i_r, i_theta)];
        });

        ddc::parallel_deepcopy(phi, get_const_field(phi_host));
    }
};
