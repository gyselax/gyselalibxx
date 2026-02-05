// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "momentscalculator.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "species_info.hpp"

MomentsCalculator::MomentsCalculator(DConstFieldVxVyVz coeffs, DConstFieldVz coeffs1D) 
: m_quadrature(coeffs) 
, m_quadrature_1D(coeffs1D) 
{}

void MomentsCalculator::operator()(DFieldX rho, DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("Hybrid_model_ChargeDensityCalculator1d3v");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
    host_t<DConstFieldSp> const kinetic_charges_host
            = charges_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);
    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho,
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += -kinetic_charges(isp) * allfdistribu(isp, idx);
                }
                return sum;
            });

    // I comment the following code, as in the hybrid model, we only compute the ion charge density.
    /*
    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }
    */
    

    Kokkos::Profiling::popRegion();
}


void MomentsCalculator::operator()(DFieldX mean_current_x, DFieldX mean_current_y, DFieldX mean_current_z, DFieldX rho, DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("Hybrid_model_MeanCurrentCalculator1d3v");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
    host_t<DConstFieldSp> const kinetic_charges_host
            = charges_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            mean_current_x,
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVx ivx(idx);
                double const vx = ddc::coordinate(ivx);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += -kinetic_charges(isp) * vx * allfdistribu(isp, idx) / rho(ix);
                }
                return sum;
            });

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            mean_current_y,
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVy ivy(idx);
                double const vy = ddc::coordinate(ivy);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += -kinetic_charges(isp) * vy * allfdistribu(isp, idx) / rho(ix);
                }
                return sum;
            });

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            mean_current_z,
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVz ivz(idx);
                double const vz = ddc::coordinate(ivz);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += -kinetic_charges(isp) * vz * allfdistribu(isp, idx) / rho(ix);
                }
                return sum;
            });

    // I comment the following code, as in the hybrid model, we only compute the ion contribution.
    /*
    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }
    */

    Kokkos::Profiling::popRegion();
}



void MomentsCalculator::operator()(DFieldX momentum_x, DFieldX momentum_y, DFieldX momentum_z, DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("Hybrid_model_MeanCurrentCalculator1d3v");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const masses_host = ddc::host_discrete_space<Species>().masses();
    host_t<DConstFieldSp> const kinetic_masses_host
            = masses_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_masses_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_masses_host);
    DConstFieldSp kinetic_masses = get_const_field(kinetic_masses_alloc);

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            momentum_x,
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVx ivx(idx);
                double const vx = ddc::coordinate(ivx);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += kinetic_masses(isp) * vx * allfdistribu(isp, idx);
                }
                return sum;
            });

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            momentum_y,
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVy ivy(idx);
                double const vy = ddc::coordinate(ivy);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += kinetic_masses(isp) * vy * allfdistribu(isp, idx);
                }
                return sum;
            });

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            momentum_z,
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVz ivz(idx);
                double const vz = ddc::coordinate(ivz);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += kinetic_masses(isp) * vz * allfdistribu(isp, idx);
                }
                return sum;
            });

    // I comment the following code, as in the hybrid model, we only compute the ion contribution.
    /*
    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }
    */

    Kokkos::Profiling::popRegion();
}



void MomentsCalculator::operator()(DFieldX kinetic_density, DConstFieldSpVxVyVzX allfdistribu, char ki) const
{
    Kokkos::Profiling::pushRegion("Hybrid_model_KineticDensityCalculator1d3v");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const masses_host = ddc::host_discrete_space<Species>().masses();
    host_t<DConstFieldSp> const kinetic_masses_host
            = masses_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_masses_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_masses_host);
    DConstFieldSp kinetic_masses = get_const_field(kinetic_masses_alloc);

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            kinetic_density,
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVx ivx(idx);
                IdxVy ivy(idx);
                IdxVz ivz(idx);
                double const vx = ddc::coordinate(ivx);
                double const vy = ddc::coordinate(ivy);
                double const vz = ddc::coordinate(ivz);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += kinetic_masses(isp) * 0.5 * (vx*vx + vy*vy + vz*vz) * allfdistribu(isp, idx);
                }
                return sum;
            });

    // I comment the following code, as in the hybrid model, we only compute the ion charge density.
    /*
    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }
    */
    

    Kokkos::Profiling::popRegion();
}



void MomentsCalculator::operator()(DFieldSpX rho, DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("Hybrid_model_ChargeDensityCalculator");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
    host_t<DConstFieldSp> const kinetic_charges_host
            = charges_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);

    /*
    for (IdxSp isp : kin_species_idx_range) {
        m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho[isp],
            KOKKOS_LAMBDA(IdxXYVxVy idx) {
                double sum = kinetic_charges(isp) * allfdistribu(isp, idx);
                return sum;
            });   
    }
    */

    ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
        m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho[isp],
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                double sum = -kinetic_charges(isp) * allfdistribu(isp, idx);
                return sum;
            }); 
    });

    

    // I comment the following code, as in the hybrid model, we only compute the ion charge density.
    /*
    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }
    */

    Kokkos::Profiling::popRegion();
}


void MomentsCalculator::operator()(DFieldSpX mean_current_x, DFieldSpX mean_current_y, DFieldSpX mean_current_z, DFieldSpX rho, DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("Hybrid_model_MeanCurrentCalculator1d3v");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
    host_t<DConstFieldSp> const kinetic_charges_host
            = charges_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);

    ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
        m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            mean_current_x[isp],
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVx ivx(idx);
                double const vx = ddc::coordinate(ivx);
                double sum = -kinetic_charges(isp) * vx * allfdistribu(isp, idx) / rho(isp,ix);
                return sum;
            });   
    });

    ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
        m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            mean_current_y[isp],
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVy ivy(idx);
                double const vy = ddc::coordinate(ivy);
                double sum = -kinetic_charges(isp) * vy * allfdistribu(isp, idx) / rho(isp,ix);
                return sum;
            });   
    });

    ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
        m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            mean_current_z[isp],
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVz ivz(idx);
                double const vz = ddc::coordinate(ivz);
                double sum = -kinetic_charges(isp) * vz * allfdistribu(isp, idx) / rho(isp,ix);
                return sum;
            });   
    });
    /*
    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            mean_velocity_x,
            KOKKOS_LAMBDA(IdxXYVxVy idx) {
                IdxXY ixy(idx);
                IdxVx ivx(idx);
                double const vx = ddc::coordinate(ivx);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += kinetic_charges(isp) * vx * allfdistribu(isp, idx) / rho(ixy);
                }
                return sum;
            });

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            mean_velocity_y,
            KOKKOS_LAMBDA(IdxXYVxVy idx) {
                IdxXY ixy(idx);
                IdxVy ivy(idx);
                double const vy = ddc::coordinate(ivy);
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += kinetic_charges(isp) * vy * allfdistribu(isp, idx) / rho(ixy);
                }
                return sum;
            });
    */

    // I comment the following code, as in the hybrid model, we only compute the ion contribution.
    /*
    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }
    */

    Kokkos::Profiling::popRegion();
}



void MomentsCalculator::operator()(DFieldSpXVxVy rho, DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("Hybrid_model_Density_after_1D_int_Calculator");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
    host_t<DConstFieldSp> const kinetic_charges_host
            = charges_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);

    /*
    for (IdxSp isp : kin_species_idx_range) {
        m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho[isp],
            KOKKOS_LAMBDA(IdxXYVxVy idx) {
                double sum = kinetic_charges(isp) * allfdistribu(isp, idx);
                return sum;
            });   
    }
    */

    ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
        m_quadrature_1D(
            Kokkos::DefaultExecutionSpace(),
            rho[isp],
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                double sum = allfdistribu(isp, idx);
                return sum;
            }); 
    });

    

    // I comment the following code, as in the hybrid model, we only compute the ion charge density.
    /*
    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }
    */

    Kokkos::Profiling::popRegion();
}



void MomentsCalculator::operator()(DFieldSpX parallel_temperature, DFieldSpX perpendicular_temperature, 
                                    DConstFieldSpX rho, DConstFieldX Bx, DConstFieldX By, DConstFieldX Bz, 
                                    DConstFieldSpX ux, DConstFieldSpX uy, DConstFieldSpX uz, 
                                    DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("Hybrid_model_TemperatureCalculator1d3v");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
    host_t<DConstFieldSp> const kinetic_charges_host
            = charges_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);

    ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
        m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            parallel_temperature[isp],
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVx ivx(idx);
                IdxVy ivy(idx);
                IdxVz ivz(idx);
                double const vx = ddc::coordinate(ivx);
                double const vy = ddc::coordinate(ivy);
                double const vz = ddc::coordinate(ivz);
                double B_norm = Kokkos::sqrt(Bx(ix) * Bx(ix) + By(ix) * By(ix) + Bz(ix) * Bz(ix) );
                double bx = Bx(ix) / B_norm;
                double by = By(ix) / B_norm;
                double bz = Bz(ix) / B_norm;
                double p_xx = (vx - ux(isp,ix))*(vx-ux(isp,ix));
                double p_xy = (vx - ux(isp,ix))*(vy-uy(isp,ix));
                double p_xz = (vx - ux(isp,ix))*(vz-uz(isp,ix));
                double p_yy = (vy - uy(isp,ix))*(vy-uy(isp,ix));
                double p_yz = (vy - uy(isp,ix))*(vz-uz(isp,ix));
                double p_zz = (vz - uz(isp,ix))*(vz-uz(isp,ix));

                double sum =  p_xx*bx*bx;
                       sum += p_yy*by*by;
                       sum += p_zz*bz*bz;
                       sum += 2.0 * p_xy*bx*by;
                       sum += 2.0 * p_xz*bx*bz;
                       sum += 2.0 * p_yz*by*bz;
                       sum *=  allfdistribu(isp, idx) / rho(isp,ix);
                return sum;
            });   
    });

    ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
        m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            perpendicular_temperature[isp],
            KOKKOS_LAMBDA(IdxXVxVyVz idx) {
                IdxX ix(idx);
                IdxVx ivx(idx);
                IdxVy ivy(idx);
                IdxVz ivz(idx);
                double const vx = ddc::coordinate(ivx);
                double const vy = ddc::coordinate(ivy);
                double const vz = ddc::coordinate(ivz);
                double B_norm = Kokkos::sqrt(Bx(ix) * Bx(ix) + By(ix) * By(ix) + Bz(ix) * Bz(ix) );
                double bx = Bx(ix) / B_norm;
                double by = By(ix) / B_norm;
                double bz = Bz(ix) / B_norm;
                double p_xx = (vx - ux(isp,ix))*(vx-ux(isp,ix));
                double p_xy = (vx - ux(isp,ix))*(vy-uy(isp,ix));
                double p_xz = (vx - ux(isp,ix))*(vz-uz(isp,ix));
                double p_yy = (vy - uy(isp,ix))*(vy-uy(isp,ix));
                double p_yz = (vy - uy(isp,ix))*(vz-uz(isp,ix));
                double p_zz = (vz - uz(isp,ix))*(vz-uz(isp,ix));

                double sum =  p_xx*bx*bx;
                       sum += p_yy*by*by;
                       sum += p_zz*bz*bz;
                       sum += 2.0 * p_xy*bx*by;
                       sum += 2.0 * p_xz*bx*bz;
                       sum += 2.0 * p_yz*by*bz;
                       sum  = 0.5 * (p_xx + p_yy + p_zz) - 0.5 * sum;
                       sum *=  allfdistribu(isp, idx) / rho(isp,ix);
                return sum;
            });   
    });

    Kokkos::Profiling::popRegion();
}




