#include "singlemodeperturbinitialization.hpp"

SingleModePerturbInitialization::SingleModePerturbInitialization(
        SpeciesInformation const& species_info,
        ViewSp<int> init_perturb_mode,
        DViewSp init_perturb_amplitude)
    : m_species_info(species_info)
    , m_init_perturb_mode(init_perturb_mode)
    , m_init_perturb_amplitude(init_perturb_amplitude)
{
}

DSpanSpXVx SingleModePerturbInitialization::operator()(DSpanSpXVx allfdistribu) const
{
    MDomainX gridx = allfdistribu.domain<MeshX>();
    MDomainVx gridvx = allfdistribu.domain<MeshVx>();
    MDomainSp gridsp = allfdistribu.domain<MeshSp>();

    // Initialization of the perturbation
    DBlockX perturbation(gridx);
    for (MCoordSp isp : gridsp) {
        perturbation_initialization(
                perturbation,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));

        // Initialization of the distribution function --> fill values
        for (MCoordSp isp : gridsp) {
            for (MCoordX ix : gridx) {
                for (MCoordVx iv : gridvx) {
                    double fdistribu_val
                            = m_species_info.maxw_values()(isp, iv) * (1. + perturbation(ix));
                    if (fdistribu_val < 1.e-60) {
                        fdistribu_val = 1.e-60;
                    }
                    allfdistribu(isp, ix, iv) = fdistribu_val;
                }
            }
        }
    }
    return allfdistribu;
}

void SingleModePerturbInitialization::perturbation_initialization(
        DSpanX perturbation,
        const int mode,
        const double perturb_amplitude) const
{
    static_assert(Dim::X::PERIODIC, "this computation for Lx is only valid for X periodic");
    auto gridx = perturbation.domain();
    const double Lx = fabs(gridx.mesh<MeshX>().step() + gridx.rmax() - gridx.rmin());
    const double kx = mode * 2. * M_PI / Lx;
    for (MCoordX ix : gridx) {
        const RCoordX x = gridx.to_real(ix);
        perturbation(ix) = perturb_amplitude * cos(kx * x);
    }
}
