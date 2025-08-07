

# File vortex\_merger\_equilibrium.hpp

[**File List**](files.md) **>** [**geometryRTheta**](dir_e9f169004bcfe9f3cb1f8a27ce024e59.md) **>** [**initialisation**](dir_1b70d60e6147eeeade38d183e3e9d318.md) **>** [**vortex\_merger\_equilibrium.hpp**](vortex__merger__equilibrium_8hpp.md)

[Go to the documentation of this file](vortex__merger__equilibrium_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "l_norm_tools.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"


template <class Mapping>
class VortexMergerEquilibria
{
private:
    Mapping const& m_mapping;
    IdxRangeRTheta const& m_grid;
    SplineRThetaBuilder const& m_builder;
    SplineRThetaEvaluatorNullBound const& m_evaluator;
    PolarSplineFEMPoissonLikeSolver<
            GridR,
            GridTheta,
            PolarBSplinesRTheta,
            SplineRThetaEvaluatorNullBound> const& m_poisson_solver;

public:
    VortexMergerEquilibria(
            Mapping const& mapping,
            IdxRangeRTheta const& grid,
            SplineRThetaBuilder const& builder,
            SplineRThetaEvaluatorNullBound const& evaluator,
            PolarSplineFEMPoissonLikeSolver<
                    GridR,
                    GridTheta,
                    PolarBSplinesRTheta,
                    SplineRThetaEvaluatorNullBound> const& poisson_solver)
        : m_mapping(mapping)
        , m_grid(grid)
        , m_builder(builder)
        , m_evaluator(evaluator)
        , m_poisson_solver(poisson_solver) {};

    void find_equilibrium(
            host_t<DFieldRTheta> sigma_host,
            host_t<DFieldRTheta> phi_eq_host,
            host_t<DFieldRTheta> rho_eq_host,
            std::function<double(double const)> const& function,
            double const phi_max, // ToDo: ADD CASE WHERE RHO_MAX IS GIVEN.
            double const tau,
            int count_max = 25) const
    {
        DFieldMemRTheta phi_star_alloc(m_grid);
        host_t<DFieldMemRTheta> ci_alloc_host(m_grid);

        IdxRangeBSRTheta idx_range_bsplinesRTheta = get_spline_idx_range(m_builder);
        Spline2DMem rho_coef_alloc(idx_range_bsplinesRTheta);

        auto rho_eq_alloc = ddc::create_mirror_view(Kokkos::DefaultExecutionSpace(), rho_eq_host);
        DFieldRTheta rho_eq(rho_eq_alloc);

        double difference_sigma(0.);
        int count = 0;

        do {
            count += 1;
            // STEP 1: compute rho^i
            ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
                rho_eq_host(irtheta) = sigma_host(irtheta) * function(phi_eq_host(irtheta));
            });


            // STEP 2: compute phi_star^i with PDE solver
            ddc::parallel_deepcopy(rho_eq, get_const_field(rho_eq_host));
            m_builder(get_field(rho_coef_alloc), get_const_field(rho_eq));
            PoissonLikeRHSFunction poisson_rhs(get_const_field(rho_coef_alloc), m_evaluator);
            m_poisson_solver(poisson_rhs, get_field(phi_star_alloc));
            auto phi_star_host = ddc::create_mirror_view_and_copy(get_field(phi_star_alloc));

            // STEP 3: compute c^i
            // If phi_max is given:
            double norm_Linf_phi_star(0.);
            ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
                double const abs_phi_star = fabs(phi_star_host(irtheta));
                norm_Linf_phi_star
                        = norm_Linf_phi_star > abs_phi_star ? norm_Linf_phi_star : abs_phi_star;
            });

            ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
                ci_alloc_host(irtheta) = phi_max / norm_Linf_phi_star;
            });


            // STEP 4: update sigma and phi
            difference_sigma = 0.;
            ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
                double const abs_diff_sigma
                        = fabs(sigma_host(irtheta) - ci_alloc_host(irtheta) * sigma_host(irtheta));
                difference_sigma
                        = difference_sigma > abs_diff_sigma ? difference_sigma : abs_diff_sigma;

                sigma_host(irtheta) = ci_alloc_host(irtheta) * sigma_host(irtheta);
                phi_eq_host(irtheta) = ci_alloc_host(irtheta) * phi_star_host(irtheta);
            });

        } while ((difference_sigma > tau) and (count < count_max));


        // STEP 1: compute rho^i
        ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
            rho_eq_host(irtheta) = sigma_host(irtheta) * function(phi_eq_host(irtheta));
        });

        // Unify at the centre point:
        IdxRangeR r_idx_range = get_idx_range<GridR>(rho_eq_host);
        IdxRangeTheta theta_idx_range = get_idx_range<GridTheta>(rho_eq_host);
        if (std::fabs(ddc::coordinate(r_idx_range.front())) < 1e-15) {
            ddc::for_each(theta_idx_range, [&](const IdxTheta itheta) {
                rho_eq_host(r_idx_range.front(), itheta)
                        = rho_eq_host(r_idx_range.front(), theta_idx_range.front());
            });
        }
    };


    void set_equilibrium(
            host_t<DFieldRTheta> rho_eq_host,
            std::function<double(double const)> function,
            double const phi_max,
            double const tau)
    {
        IdxRangeRTheta grid = get_idx_range<GridR, GridTheta>(rho_eq_host);

        // Equilibrium:
        host_t<DFieldMemRTheta> sigma_0_alloc_host(grid);
        host_t<DFieldMemRTheta> phi_eq_alloc_host(grid);
        const double sig = 0.3;
        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            const CoordRTheta coord_rtheta(ddc::coordinate(irtheta));
            const CoordXY coord_xy(m_mapping(coord_rtheta));
            const double x = ddc::get<X>(coord_xy);
            const double y = ddc::get<Y>(coord_xy);
            sigma_0_alloc_host(irtheta) = sig;
            phi_eq_alloc_host(irtheta) = std::exp(-(x * x + y * y) / (2 * sig * sig));
        });
        find_equilibrium(
                get_field(sigma_0_alloc_host),
                get_field(phi_eq_alloc_host),
                get_field(rho_eq_host),
                function,
                phi_max,
                tau);
    };
};
```


