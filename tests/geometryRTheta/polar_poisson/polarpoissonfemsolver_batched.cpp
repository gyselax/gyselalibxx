// SPDX-License-Identifier: MIT
//
// Tests the batched PolarSplineFEMPoissonLikeSolver by comparing the result of a batched
// solve with the results of unbatched solves carried out for each batch index.
//
// Each batch uses different values of alpha, beta and rho.

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "discrete_poloidal_cs_spline_mapping.hpp"
#include "discrete_poloidal_cs_spline_mapping_builder.hpp"
#include "geometry_r_theta.hpp"
#include "mesh_builder.hpp"
#include "polar_spline_fem_poisson_like_solver.hpp"
#include "spline_definitions_r_theta.hpp"

namespace {

struct A
{
};
struct B
{
};

struct GridA : UniformGridBase<A>
{
};

struct GridB : UniformGridBase<B>
{
};

using Mapping = CircularToCartesian<R, Theta, X, Y>;
using DiscreteMappingBuilder
        = DiscretePoloidalCSSplineMappingBuilder<X, Y, SplineInterpolatorRTheta>;
using DiscreteMapping = typename DiscreteMappingBuilder::MappingType;

using IdxRangeBatch = IdxRange<GridA, GridB>;
using IdxBatch = Idx<GridA, GridB>;
using IdxRangeBatched = IdxRange<GridA, GridB, GridR, GridTheta>;
using IdxBatched = Idx<GridA, GridB, GridR, GridTheta>;

using DFieldMemBatched = DFieldMem<IdxRangeBatched>;
using DFieldBatched = DField<IdxRangeBatched>;

using PoissonSolver = PolarSplineFEMPoissonLikeSolver<
        GridR,
        GridTheta,
        PolarBSplinesRTheta,
        SplineInterpolatorRTheta,
        DiscreteMapping>;

using BatchedPoissonSolver = PolarSplineFEMPoissonLikeSolver<
        GridR,
        GridTheta,
        PolarBSplinesRTheta,
        SplineInterpolatorRTheta,
        DiscreteMapping,
        IdxRangeBatched>;

constexpr double res_tol = 1e-13;

KOKKOS_FUNCTION double alpha_func(double r, double a)
{
    return (1.0 + 0.5 * a) * Kokkos::exp(-Kokkos::tanh((r - 0.7) / 0.05));
}

KOKKOS_FUNCTION double beta_func(double r, double a, double b)
{
    return (1.0 + b) / alpha_func(r, a);
}

KOKKOS_FUNCTION double rho_func(double r, double theta, double a, double b)
{
    return (1.0 + a) * (1.0 - r * r) * (1.0 + 0.3 * (b + 1.0) * r * Kokkos::cos(theta));
}

/// The RHS of the equation for a given batch index.
class BatchedRHS
{
public:
    KOKKOS_FUNCTION double operator()(IdxBatch idx_batch, CoordRTheta const& coord) const
    {
        double const a = ddc::coordinate(Idx<GridA>(idx_batch));
        double const b = ddc::coordinate(Idx<GridB>(idx_batch));
        return rho_func(ddc::get<R>(coord), ddc::get<Theta>(coord), a, b);
    }
};

/// The RHS of the equation for a fixed batch.
class UnbatchedRHS
{
    double m_a;
    double m_b;

public:
    UnbatchedRHS(double a, double b) : m_a(a), m_b(b) {}

    KOKKOS_FUNCTION double operator()(CoordRTheta const& coord) const
    {
        return rho_func(ddc::get<R>(coord), ddc::get<Theta>(coord), m_a, m_b);
    }
};

} // namespace

void test_PolarPoissonFEMBatched__MatchesUnbatched()
{
    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_ncells(16);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(32);

    std::vector<CoordR> r_break_points = build_uniform_break_points(r_min, r_max, r_ncells);
    std::vector<CoordTheta> theta_break_points
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);

    ddc::init_discrete_space<GridA>(Coord<A>(0.0), 1.0);
    ddc::init_discrete_space<GridB>(Coord<B>(0.0), 1.0);
    ddc::init_discrete_space<BSplinesR>(r_break_points);
    ddc::init_discrete_space<BSplinesTheta>(theta_break_points);
    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR const idx_range_r(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta const idx_range_theta(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta const idx_range(idx_range_r, idx_range_theta);

    IdxRangeBatch const idx_range_batch(IdxBatch(0, 0), IdxStep<GridA, GridB>(2, 2));
    IdxRangeBatched const idx_range_batched(idx_range_batch, idx_range);

    SplineInterpolatorRTheta const interpolator(idx_range);

    const Mapping mapping;
    DiscreteMappingBuilder const
            discrete_mapping_builder(Kokkos::DefaultExecutionSpace(), mapping, interpolator);
    DiscreteMapping const discrete_mapping = discrete_mapping_builder();

    ddc::init_discrete_space<PolarBSplinesRTheta>(discrete_mapping);

    // Initialise the batched coefficients and RHS
    DFieldMemBatched alpha_alloc(idx_range_batched);
    DFieldMemBatched beta_alloc(idx_range_batched);
    DFieldMemBatched rho_alloc(idx_range_batched);
    DFieldBatched alpha = get_field(alpha_alloc);
    DFieldBatched beta = get_field(beta_alloc);
    DFieldBatched rho = get_field(rho_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_batched,
            KOKKOS_LAMBDA(IdxBatched const idx) {
                double const a = ddc::coordinate(Idx<GridA>(idx));
                double const b = ddc::coordinate(Idx<GridB>(idx));
                double const r = ddc::coordinate(Idx<GridR>(idx));
                double const theta = ddc::coordinate(Idx<GridTheta>(idx));
                alpha(idx) = alpha_func(r, a);
                beta(idx) = beta_func(r, a, b);
                rho(idx) = rho_func(r, theta, a, b);
            });

    // Solve the batched equations
    BatchedPoissonSolver
            batched_solver(discrete_mapping, interpolator, idx_range_batched, 1000, res_tol);
    batched_solver.update_coefficients(get_const_field(alpha), get_const_field(beta));

    DFieldMemBatched phi_from_field_alloc(idx_range_batched);
    DFieldMemBatched phi_from_func_alloc(idx_range_batched);
    batched_solver(get_field(phi_from_field_alloc), get_const_field(rho));
    batched_solver(get_field(phi_from_func_alloc), BatchedRHS());

    auto phi_from_field_host = ddc::create_mirror_view_and_copy(get_field(phi_from_field_alloc));
    auto phi_from_func_host = ddc::create_mirror_view_and_copy(get_field(phi_from_func_alloc));

    // Compare with unbatched solves
    PoissonSolver solver(discrete_mapping, interpolator, idx_range, 1000, res_tol);
    DFieldMemRTheta alpha_ref(idx_range);
    DFieldMemRTheta beta_ref(idx_range);
    DFieldMemRTheta rho_ref(idx_range);
    DFieldMemRTheta phi_from_field_ref(idx_range);
    DFieldMemRTheta phi_from_func_ref(idx_range);

    ddc::host_for_each(idx_range_batch, [&](IdxBatch idx_batch) {
        double const a = ddc::coordinate(Idx<GridA>(idx_batch));
        double const b = ddc::coordinate(Idx<GridB>(idx_batch));

        ddc::parallel_deepcopy(get_field(alpha_ref), alpha[idx_batch]);
        ddc::parallel_deepcopy(get_field(beta_ref), beta[idx_batch]);
        ddc::parallel_deepcopy(get_field(rho_ref), rho[idx_batch]);

        solver.update_coefficients(get_const_field(alpha_ref), get_const_field(beta_ref));
        solver(get_field(phi_from_field_ref), get_const_field(rho_ref));
        solver(get_field(phi_from_func_ref), UnbatchedRHS(a, b));

        auto phi_from_field_ref_host
                = ddc::create_mirror_view_and_copy(get_field(phi_from_field_ref));
        auto phi_from_func_ref_host
                = ddc::create_mirror_view_and_copy(get_field(phi_from_func_ref));

        double max_phi = 0.0;
        ddc::host_for_each(idx_range, [&](IdxRTheta idx) {
            max_phi = std::max(max_phi, std::abs(phi_from_func_ref_host(idx)));
        });
        // Check that the batches are not trivial
        EXPECT_GT(max_phi, 1e-3);

        ddc::host_for_each(idx_range, [&](IdxRTheta idx) {
            EXPECT_NEAR(
                    phi_from_field_host(idx_batch, idx),
                    phi_from_field_ref_host(idx),
                    1e-10 * max_phi);
            EXPECT_NEAR(
                    phi_from_func_host(idx_batch, idx),
                    phi_from_func_ref_host(idx),
                    1e-10 * max_phi);
        });
    });
}

TEST(PolarPoissonFEMBatched, MatchesUnbatched)
{
    test_PolarPoissonFEMBatched__MatchesUnbatched();
}
