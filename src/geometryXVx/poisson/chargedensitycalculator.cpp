// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "chargedensitycalculator.hpp"
#include "quadrature.hpp"
#include "simpson_quadrature.hpp"



ChargeDensityCalculator::ChargeDensityCalculator(Quadrature<IDimVx> const& quad) : m_quad(quad) {}

device_t<DSpanX> ChargeDensityCalculator::operator()(
        device_t<DSpanX> const rho_device,
        device_t<DViewSpXVx> const allfdistribu_device) const
{
    Kokkos::Profiling::pushRegion("ChargeDensityCalculator");
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu_device);
    auto rho_alloc = ddc::create_mirror_view_and_copy(rho_device);
    ddc::ChunkSpan const allfdistribu = allfdistribu_alloc.span_cview();
    ddc::ChunkSpan const rho = rho_alloc.span_view();

    DFieldVx f_vx_slice(allfdistribu.domain<IDimVx>());
    IndexSp const last_kin_species = allfdistribu.domain<IDimSp>().back();
    IndexSp const last_species = ddc::discrete_space<IDimSp>().charges().domain().back();
    double chargedens_adiabspecies = 0.;
    if (last_kin_species != last_species) {
        chargedens_adiabspecies = double(charge(last_species));
    }

    ddc::for_each(rho.domain(), [&](IndexX const ix) {
        rho(ix) = chargedens_adiabspecies;
        ddc::for_each(ddc::get_domain<IDimSp>(allfdistribu), [&](IndexSp const isp) {
            ddc::deepcopy(f_vx_slice, allfdistribu[isp][ix]);

            rho(ix) += charge(isp) * m_quad(f_vx_slice.span_view());
        });
    });
    ddc::deepcopy(rho_device, rho);
    Kokkos::Profiling::popRegion();

    return rho_device;
}
