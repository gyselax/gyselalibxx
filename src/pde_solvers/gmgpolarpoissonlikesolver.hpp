#include <GMGPolar/gmgpolar.h>

namespace GMGPolarTools {

template <class ToPhysicalMapping>
class MappingToDomainGeometry : public DomainGeometry
{
    using R = typename ToPhysicalMapping::curvilinear_tag_r;
    using Theta = typename ToPhysicalMapping::curvilinear_tag_r;

    using X = typename ToPhysicalMapping::cartesian_tag_x;
    using Y = typename ToPhysicalMapping::cartesian_tag_y;

    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

private:
    ToPhysicalMapping m_to_physical;

public:
    explicit MappingToDomainGeometry(ToPhysicalMapping to_physical) : m_to_physical(to_physical) {}

    // In earlier versions denoted by 'x' and 'y'.
    double Fx(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const final
    {
        return Coord<X>(m_to_physical(Coord<R, Theta>(r, theta)));
    }
    double Fy(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const final
    {
        return Coord<Y>(m_to_physical(Coord<R, Theta>(r, theta)));
    }

    // In earlier versions denoted by 'Jrr', 'Jtr', 'Jrt' and 'Jtt'.
    double dFx_dr(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const final
    {
        return m_to_physical.template jacobian_component<X, R_cov>(Coord<R, Theta>(r, theta));
    }
    double dFy_dr(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const final
    {
        return m_to_physical.template jacobian_component<Y, R_cov>(Coord<R, Theta>(r, theta));
    }
    double dFx_dt(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const final
    {
        return m_to_physical.template jacobian_component<X, Theta_cov>(Coord<R, Theta>(r, theta));
    }
    double dFy_dt(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const final
    {
        return m_to_physical.template jacobian_component<Y, Theta_cov>(Coord<R, Theta>(r, theta));
    }
};

template <class R, class Theta, class RHSFunction>
class SourceTermWrapper : public SourceTerm
{
private:
    RHSFunction const& m_rhs;

public:
    SourceTermWrapper(RHSFunction const& rhs) : m_rhs(rhs) {}

    inline double rhs_f(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const final
    {
        return m_rhs(Coord<R, Theta>(r, theta));
    }
};

class HomogeneousDirichletBoundaryConditions : public BoundaryConditions
{
public:
    double u_D(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const final
    {
        return 0.0;
    }
    // Only used if DirBC_Interior = true
    double u_D_Interior(
            const double& r,
            const double& theta,
            const double& sin_theta,
            const double& cos_theta) const
    {
        assert(false);
        return 0.0;
    }
};

template <class RadialSplineEvaluator_host, class BSplinesR>
class PolarPoissonLikeCoefficients : DensityProfileCoefficients
{
    static_assert(std::is_same_v<RadialSplineEvaluator_host::memory_space, Kokkos::HostSpace>);
    using R = typename BSplinesR::continuous_dimension_type;
    using CoordR = Coord<R>;

    using DConstSplineR_host = DConstField<IdxRange<BSplinesR>, Kokkos::HostSpace>;

private:
    RadialSplineEvaluator_host m_evaluator_r;
    DConstSplineR_host m_coeff_alpha;
    DConstSplineR_host m_coeff_beta;

public:
    PolarPoissonLikeCoefficients(
            RadialSplineEvaluator_host evaluator_r,
            DConstSplineR_host coeff_alpha,
            DConstSplineR_host coeff_beta)
        : m_evaluator_r(evaluator_r)
        , m_coeff_alpha(coeff_alpha)
        , m_coeff_beta(coeff_beta)
    {
    }

    double alpha(const double& r) const final
    {
        return m_evaluator_r(CoordR(r), m_coeff_alpha);
    }
    double beta(const double& r) const final
    {
        return m_evaluator_r(CoordR(r), m_coeff_beta);
    }

    // Only used in custom mesh generation -> refinement_radius
    double getAlphaJump() const
    {
        assert(false);
    }
};

} // namespace GMGPolarTools

template <
        class ToPhysicalMapping,
        class GridR,
        class GridTheta,
        class BSplinesR,
        class RadialSplineEvaluator_host>
class GMGPolarPoissonLikeSolver
{
    using R = typename ToPhysicalMapping::curvilinear_tag_r;
    using Theta = typename ToPhysicalMapping::curvilinear_tag_r;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxStepRTheta = IdxStep<GridR, GridTheta>;

    using DFieldRTheta = DField<IdxRangeRTheta>;

    using DConstSplineR_host
            = DConstField<IdxRange<BSplinesR>, Kokkos::DefaultExecutionSpace::memory_space>;

private:
    ToPhysicalMapping m_to_physical;
    IdxRangeRTheta m_idx_range;

    RadialSplineEvaluator_host m_evaluator_r;
    DConstSplineR_host m_coeff_alpha;
    DConstSplineR_host m_coeff_beta;

    std::string m_radial_grid_file;
    std::string m_poloidal_grid_file;

public:
    GMGPolarPoissonLikeSolver(
            DConstSplineR_host coeff_alpha,
            DConstSplineR_host coeff_beta,
            ToPhysicalMapping to_physical,
            IdxRangeRTheta idx_range,
            RadialSplineEvaluator_host evaluator_r)
        : m_to_physical(to_physical)
        , m_idx_range(idx_range)
        , m_evaluator_r(evaluator_r)
        , m_coeff_alpha(coeff_alpha)
        , m_coeff_beta(coeff_beta)
        , m_radial_grid_file("gmgpolar_input_r.txt")
        , m_poloidal_grid_file("gmgpolar_input_r.txt")
    {
        std::ofstream r_grid_file(m_radial_grid_file);
        ddc::for_each(ddc::select<GridR>(idx_range), [&](IdxR ir) {
            double r_coord = ddc::coordinate(ir);
            r_grid_file << r_coord << std::endl;
        });
        r_grid_file.close();
        std::ofstream theta_grid_file(m_poloidal_grid_file);

        ddc::for_each(ddc::select<GridTheta>(idx_range), [&](IdxTheta it) {
            double theta_coord = ddc::coordinate(it);
            theta_grid_file << theta_coord << std::endl;
        });
        theta_grid_file.close();
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
     *      class. It must be an object with an operator() which evaluates a
     *      CoordRTheta and can be called from GPU.
     * @param[inout] phi
     *      The values of the solution @f$\phi@f$ on the given coords_eval.
     */
    template <class CPURHSFunction>
    void operator()(CPURHSFunction const& rhs, DFieldRTheta phi) const
    {
        static_assert(
                std::is_invocable_r_v<double, RHSFunction, CoordRTheta>,
                "RHSFunction must have an operator() which takes a coordinate and returns a "
                "double");
        assert(get_idx_range(phi) == m_idx_range);
        GMGPolar
                solver(std::make_unique<MappingToDomainGeometry<ToPhysicalMapping>>(m_to_physical),
                       std::make_unique<PolarPoissonLikeCoefficients<
                               RadialSplineEvaluator_host,
                               BSplinesR>>(m_evaluator_r, m_coeff_alpha, m_coeff_beta),
                       std::make_unique<HomogeneousDirichletBoundaryConditions>(),
                       std::make_unique<SourceTermWrapper<R, Theta, RHSFunction>>(rhs));

        solver.R0(ddc::coordinate(get_idx_range<GridR>(phi).front()));
        solver.Rmax(ddc::coordinate(get_idx_range<GridR>(phi).back()));

        // To use a grid provided by the user, currently we must read from file
        solver.load_grid_file(true);
        solver.file_grid_radii(m_radial_grid_file);
        solver.file_grid_angles(m_poloidal_grid_file);

        solver.DirBC_Interior(false); // Use across-origin calculation

        // ----------------------------------------------------------------
        // Parameters to be identified and values confirmed
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
        solver.relativeTolerance(1e-50);
        // ----------------------------------------------------------------

        // Perform setup and solve
        solver.setup();
        solver.solve();

        auto phi_host = ddc::create_mirror_view(phi);

        Vector<double>& solution = solver.solution();
        const PolarGrid& grid = solver.grid();
        ddc::parallel_for_each(
                Kokkos::DefaultHostExecutionSpace(),
                m_idx_range,
                [&](IdxRTheta idx) {
                    IdxStepRTheta offset(idx - m_idx_range.front());
                    int i_r = ddc::select<GridR>(offset);
                    int i_theta = ddc::select<GridTheta>(offset);
                    phi_host(idx) = solution[grid.index(i_r, i_theta)];
                });

        ddc::parallel_deepcopy(phi, get_const_field(phi_host));
    }
};
