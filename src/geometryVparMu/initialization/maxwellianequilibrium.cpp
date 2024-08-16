// SPDX-License-Identifier: MIT

#include <cmath>

#include <ddc/ddc.hpp>

#include "maxwellianequilibrium.hpp"

MaxwellianEquilibrium::MaxwellianEquilibrium(
        host_t<DFieldMemSp> mass,
        host_t<DFieldMemSp> density_eq,
        host_t<DFieldMemSp> temperature_eq,
        host_t<DFieldMemSp> mean_velocity_eq,
        double magnetic_field = 1.0)
    : m_mass(std::move(mass))
    , m_density_eq(std::move(density_eq))
    , m_temperature_eq(std::move(temperature_eq))
    , m_mean_velocity_eq(std::move(mean_velocity_eq))
    , m_magnetic_field(magnetic_field)
{
}

DFieldSpVparMu MaxwellianEquilibrium::operator()(DFieldSpVparMu const allfequilibrium) const
{
    IdxRangeSp const idxrange_sp = get_idx_range<Species>(allfequilibrium);
    IdxRangeVparMu const idxrange_vparmu = get_idx_range<GridVpar, GridMu>(allfequilibrium);

    // Initialization of the maxwellian
    DFieldMemVparMu maxwellian_alloc(idxrange_vparmu);
    DFieldVparMu maxwellian = get_field(maxwellian_alloc);
    double const magnetic_field(1.0);
    ddc::for_each(idxrange_sp, [&](IdxSp const isp) {
        compute_maxwellian(
                maxwellian,
                m_mass(isp),
                m_density_eq(isp),
                m_temperature_eq(isp),
                m_mean_velocity_eq(isp),
                magnetic_field);

        ddc::parallel_deepcopy(allfequilibrium[isp], maxwellian);
    });
    return allfequilibrium;
}

MaxwellianEquilibrium MaxwellianEquilibrium::init_from_input(
        IdxRangeSp idx_range_kinsp,
        PC_tree_t const& yaml_input_file)
{
    host_t<DFieldMemSp> mass(idx_range_kinsp);
    host_t<DFieldMemSp> density_eq(idx_range_kinsp);
    host_t<DFieldMemSp> temperature_eq(idx_range_kinsp);
    host_t<DFieldMemSp> mean_velocity_eq(idx_range_kinsp);
    double const magnetic_field = 1.0;

    for (IdxSp const isp : idx_range_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(yaml_input_file, ".SpeciesInfo[%d]", isp.uid());

        mass(isp) = PCpp_double(conf_isp, ".mass");
        density_eq(isp) = PCpp_double(conf_isp, ".density_eq");
        temperature_eq(isp) = PCpp_double(conf_isp, ".temperature_eq");
        mean_velocity_eq(isp) = PCpp_double(conf_isp, ".mean_velocity_eq");
    }

    return MaxwellianEquilibrium(
            std::move(mass),
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq),
            magnetic_field);
}

void MaxwellianEquilibrium::compute_maxwellian(
        DFieldVparMu const fMaxwellian,
        double const mass,
        double const density,
        double const temperature,
        double const mean_velocity,
        double const magnetic_field)
{
    double const mass_on_2piT = mass / (2. * M_PI * temperature);
    double const coeff_maxw = Kokkos::sqrt(mass_on_2piT) * mass_on_2piT;
    IdxRangeVparMu const idxrange_vparmu = get_idx_range<GridVpar, GridMu>(fMaxwellian);

    // Compute fM(vpar,mu) = (2*PI*T)**1.5*n*exp(-energy) with
    //  - n the density, T the temperature and Upar the mean velocity
    //  - B the magnetic field and
    //  - energy = (0.5*(vpar-Upar)**2+mu*B)/T
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_vparmu,
            KOKKOS_LAMBDA(IdxVparMu const ivparmu) {
                double const vpar = ddc::coordinate(ddc::select<GridVpar>(ivparmu));
                double const mu = ddc::coordinate(ddc::select<GridMu>(ivparmu));
                double const energy = (0.5 * (vpar - mean_velocity) * (vpar - mean_velocity)
                                       + mu * magnetic_field)
                                      / temperature;
                fMaxwellian(ivparmu) = coeff_maxw * density * Kokkos::exp(-energy);
            });
}
