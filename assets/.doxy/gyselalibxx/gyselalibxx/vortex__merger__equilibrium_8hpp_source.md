

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
    SplineRThetaBuilder_host const& m_builder;
    SplineRThetaEvaluatorNullBound_host const& m_evaluator;
    PolarSplineFEMPoissonLikeSolver<
            GridR,
            GridTheta,
            PolarBSplinesRTheta,
            SplineRThetaEvaluatorNullBound> const& m_poisson_solver;

public:
    VortexMergerEquilibria(
            Mapping const& mapping,
            IdxRangeRTheta const& grid,
            SplineRThetaBuilder_host const& builder,
            SplineRThetaEvaluatorNullBound_host const& evaluator,
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
            host_t<DFieldRTheta> sigma,
            host_t<DFieldRTheta> phi_eq,
            host_t<DFieldRTheta> rho_eq,
            std::function<double(double const)> const& function,
            double const phi_max, // ToDo: ADD CASE WHERE RHO_MAX IS GIVEN.
            double const tau,
            int count_max = 25) const
    {
        DFieldMemRTheta phi_star(m_grid);
        host_t<DFieldMemRTheta> ci(m_grid);

        IdxRangeBSRTheta idx_range_bsplinesRTheta = get_spline_idx_range(m_builder);
        host_t<Spline2DMem> rho_coef(idx_range_bsplinesRTheta);

        double difference_sigma(0.);
        int count = 0;

        do {
            count += 1;
            // STEP 1: compute rho^i
            ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
                rho_eq(irtheta) = sigma(irtheta) * function(phi_eq(irtheta));
            });


            // STEP 2: compute phi_star^i with PDE solver
            m_builder(get_field(rho_coef), get_const_field(rho_eq));
            PoissonLikeRHSFunction poisson_rhs(get_const_field(rho_coef), m_evaluator);
            m_poisson_solver(poisson_rhs, get_field(phi_star));
            auto phi_star_host = ddc::create_mirror_view_and_copy(get_field(phi_star));

            // STEP 3: compute c^i
            // If phi_max is given:
            double norm_Linf_phi_star(0.);
            ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
                double const abs_phi_star = fabs(phi_star_host(irtheta));
                norm_Linf_phi_star
                        = norm_Linf_phi_star > abs_phi_star ? norm_Linf_phi_star : abs_phi_star;
            });

            ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
                ci(irtheta) = phi_max / norm_Linf_phi_star;
            });


            // STEP 4: update sigma and phi
            difference_sigma = 0.;
            ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
                double const abs_diff_sigma = fabs(sigma(irtheta) - ci(irtheta) * sigma(irtheta));
                difference_sigma
                        = difference_sigma > abs_diff_sigma ? difference_sigma : abs_diff_sigma;

                sigma(irtheta) = ci(irtheta) * sigma(irtheta);
                phi_eq(irtheta) = ci(irtheta) * phi_star_host(irtheta);
            });

        } while ((difference_sigma > tau) and (count < count_max));


        // STEP 1: compute rho^i
        ddc::for_each(m_grid, [&](IdxRTheta const irtheta) {
            rho_eq(irtheta) = sigma(irtheta) * function(phi_eq(irtheta));
        });

        // Unify at the centre point:
        IdxRangeR r_idx_range = get_idx_range<GridR>(rho_eq);
        IdxRangeTheta theta_idx_range = get_idx_range<GridTheta>(rho_eq);
        if (std::fabs(ddc::coordinate(r_idx_range.front())) < 1e-15) {
            ddc::for_each(theta_idx_range, [&](const IdxTheta itheta) {
                rho_eq(r_idx_range.front(), itheta)
                        = rho_eq(r_idx_range.front(), theta_idx_range.front());
            });
        }
    };


    void set_equilibrium(
            host_t<DFieldRTheta> rho_eq,
            std::function<double(double const)> function,
            double const phi_max,
            double const tau)
    {
        IdxRangeRTheta grid = get_idx_range<GridR, GridTheta>(rho_eq);

        // Equilibrium:
        host_t<DFieldMemRTheta> sigma_0(grid);
        host_t<DFieldMemRTheta> phi_eq(grid);
        const double sig = 0.3;
        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            const CoordRTheta coord_rtheta(ddc::coordinate(irtheta));
            const CoordXY coord_xy(m_mapping(coord_rtheta));
            const double x = ddc::get<X>(coord_xy);
            const double y = ddc::get<Y>(coord_xy);
            sigma_0(irtheta) = sig;
            phi_eq(irtheta) = std::exp(-(x * x + y * y) / (2 * sig * sig));
        });
        find_equilibrium(
                get_field(sigma_0),
                get_field(phi_eq),
                get_field(rho_eq),
                function,
                phi_max,
                tau);
    };
};
```


