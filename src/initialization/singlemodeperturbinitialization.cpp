// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "singlemodeperturbinitialization.hpp"

SingleModePerturbInitialization::SingleModePerturbInitialization(
        DViewSpVx fequilibrium,
        ViewSp<int> const init_perturb_mode,
        DViewSp const init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
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
                    double fdistribu_val = m_fequilibrium(isp, iv) * (1. + perturbation(ix));
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

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
template <class IDim>
constexpr std::enable_if_t<!IDim::continuous_dimension_type::PERIODIC, double>
total_interval_length(DiscreteDomain<IDim> const& dom)
{
    return fabs(rlength(dom));
}

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
template <class RDim>
constexpr std::enable_if_t<RDim::PERIODIC, double> total_interval_length(
        DiscreteDomain<UniformPointSampling<RDim>> const& dom)
{
    return fabs(rlength(dom) + step<UniformPointSampling<RDim>>());
}

void SingleModePerturbInitialization::perturbation_initialization(
        DSpanX const perturbation,
        int const mode,
        double const perturb_amplitude) const
{
    IDomainX const gridx = perturbation.domain();
    double const Lx = total_interval_length(gridx);

    double const kx = mode * 2. * M_PI / Lx;
    for (IndexX const ix : gridx) {
        CoordX const x = coordinate(ix);
        perturbation(ix) = perturb_amplitude * cos(kx * x);
    }
}
