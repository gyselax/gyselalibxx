// SPDX-License-Identifier: MIT

#include <ddc/discretization>
#include <ddc/for_each>

#include "singlemodeperturbinitialization.hpp"

SingleModePerturbInitialization::SingleModePerturbInitialization(
        SpeciesInformation const& species_info,
        ViewSp<int> const init_perturb_mode,
        DViewSp const init_perturb_amplitude)
    : m_species_info(species_info)
    , m_init_perturb_mode(init_perturb_mode)
    , m_init_perturb_amplitude(init_perturb_amplitude)
{
}

DSpanSpXVx SingleModePerturbInitialization::operator()(DSpanSpXVx const allfdistribu) const
{
    IDomainX const gridx = allfdistribu.domain<IDimX>();
    IDomainVx const gridvx = allfdistribu.domain<IDimVx>();
    IDomainSp const gridsp = allfdistribu.domain<IDimSp>();

    // Initialization of the perturbation
    DFieldX perturbation(gridx);
    for_each(gridsp, [&](IndexSp const isp) {
        perturbation_initialization(
                perturbation,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));

        // Initialization of the distribution function --> fill values
        for_each(gridsp, [&](IndexSp const isp) {
            for_each(gridx, [&](IndexX const ix) {
                for_each(gridvx, [&](IndexVx const iv) {
                    double fdistribu_val
                            = m_species_info.maxw_values()(isp, iv) * (1. + perturbation(ix));
                    if (fdistribu_val < 1.e-60) {
                        fdistribu_val = 1.e-60;
                    }
                    allfdistribu(isp, ix, iv) = fdistribu_val;
                });
            });
        });
    });
    return allfdistribu;
}

void SingleModePerturbInitialization::perturbation_initialization(
        DSpanX const perturbation,
        int const mode,
        double const perturb_amplitude) const
{
    IDomainX const gridx = perturbation.domain();
#if defined(ENABLE_PERIODIC_RDIMX)
#if ENABLE_PERIODIC_RDIMX
    double const Lx = fabs(rlength(gridx) + step<IDimX>());
#else
    double const Lx = fabs(rlength(gridx));
#endif
#endif
    double const kx = mode * 2. * M_PI / Lx;
    for (IndexX const ix : gridx) {
        CoordX const x = to_real(ix);
        perturbation(ix) = perturb_amplitude * cos(kx * x);
    }
}
