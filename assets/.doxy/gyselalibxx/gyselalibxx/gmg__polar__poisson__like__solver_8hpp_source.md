

# File gmg\_polar\_poisson\_like\_solver.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**gmg\_polar\_poisson\_like\_solver.hpp**](gmg__polar__poisson__like__solver_8hpp.md)

[Go to the documentation of this file](gmg__polar__poisson__like__solver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <vector>

#include <GMGPolar/gmgpolar.h>

#include "ddc_alias_inline_functions.hpp"
#include "i_interpolation_builder.hpp"
#include "ipolar_poisson_like_solver.hpp"

namespace GMGPolarTools {

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
    explicit MappingToDomainGeometry(ToPhysicalMapping to_physical) : m_to_physical(to_physical) {}

    KOKKOS_INLINE_FUNCTION double Fx(const double& r, const double& theta) const
    {
        return Coord<X>(m_to_physical(Coord<R, Theta>(r, theta)));
    }
    KOKKOS_INLINE_FUNCTION double Fy(const double& r, const double& theta) const
    {
        return Coord<Y>(m_to_physical(Coord<R, Theta>(r, theta)));
    }
    KOKKOS_INLINE_FUNCTION double dFx_dr(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<X, R_cov>(Coord<R, Theta>(r, theta));
    }
    KOKKOS_INLINE_FUNCTION double dFy_dr(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<Y, R_cov>(Coord<R, Theta>(r, theta));
    }
    KOKKOS_INLINE_FUNCTION double dFx_dtheta(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<X, Theta_cov>(Coord<R, Theta>(r, theta));
    }
    KOKKOS_INLINE_FUNCTION double dFy_dtheta(const double& r, const double& theta) const
    {
        return m_to_physical.template jacobian_component<Y, Theta_cov>(Coord<R, Theta>(r, theta));
    }
};

class HomogeneousDirichletBoundaryConditions
{
public:
    static KOKKOS_INLINE_FUNCTION double u_D(const double& r, const double& theta)
    {
        return 0.0;
    }
    static KOKKOS_INLINE_FUNCTION double u_D_Interior(const double& r, const double& theta)
    {
        // Only used if DirBC_Interior = true
        assert(false);
        return 0.0;
    }
};

template <class EvaluatorType, class IdxRangeCoeff, class CoordRTheta>
class PolarPoissonLikeCoefficients
{
    using DConstCoeffRTheta = DConstField<IdxRangeCoeff>;

private:
    EvaluatorType m_evaluator;
    DConstCoeffRTheta m_coeff_alpha;
    DConstCoeffRTheta m_coeff_beta;

public:
    PolarPoissonLikeCoefficients(
            EvaluatorType evaluator,
            DConstCoeffRTheta coeff_alpha,
            DConstCoeffRTheta coeff_beta)
        : m_evaluator(evaluator)
        , m_coeff_alpha(coeff_alpha)
        , m_coeff_beta(coeff_beta)
    {
    }

    KOKKOS_INLINE_FUNCTION double alpha(const double& r, const double& theta) const
    {
        return m_evaluator(CoordRTheta(r, theta), m_coeff_alpha);
    }
    KOKKOS_INLINE_FUNCTION double beta(const double& r, const double& theta) const
    {
        return m_evaluator(CoordRTheta(r, theta), m_coeff_beta);
    }

    static double getAlphaJump()
    {
        assert(false);
        return 0.0;
    }
};

} // namespace GMGPolarTools

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
```


