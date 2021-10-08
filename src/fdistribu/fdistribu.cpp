#include <cmath>
#include <iostream>

#include "fdistribu.h"

using std::cos, std::fabs, std::sqrt, std::exp;

/*
 Initialization of the distribution function
*/
/*DistributionFunction::DistributionFunction(MDomainXVx const &domXVx):
    DBlockXVx values(domXVx)
{
    fdistribu_initialization(values);
}*/
DistributionFunction::DistributionFunction(
        int const species_charge,
        double const species_mass,
        double const species_n_eq,
        double const species_T_eq,
        double const species_u_eq,
        MDomainXVx const& domXVx)
    : charge(species_charge)
    , mass(species_mass)
    , density_eq(species_n_eq)
    , temperature_eq(species_T_eq)
    , mean_velocity_eq(species_u_eq)
    , Maxw_values(select<MeshVx>(domXVx))
    , values(domXVx)
{
}

/*
 Computing the non-centered Maxwellian function as
   fM(v) = n/(sqrt(2*PI*T))*exp(-(v-u)**2/(2*T))
  with n the density and T the temperature and
  where u is the mean velocity
*/
void Maxwellian_initialization(
        const double density,
        const double temperature,
        const double mean_velocity,
        DSpanVx fMaxwellian)
{
    const double inv_sqrt_2piT = 1. / sqrt(2. * M_PI * temperature);
    auto gridvx = fMaxwellian.domain<MeshVx>();
    for (MCoordVx iv : gridvx) {
        const RCoordVx v = gridvx.to_real(iv);
        fMaxwellian(iv) = density * inv_sqrt_2piT
                          * exp(-(v - mean_velocity) * (v - mean_velocity) / (2. * temperature));
    }
}

/*
 Initialization of the perturbation
*/
void perturbation_initialization(
        const int mode,
        const double perturb_amplitude,
        DSpanX perturbation)
{
    auto gridx = perturbation.domain();
    const double Lx = fabs(gridx.mesh<MeshX>().step() + gridx.rmax() - gridx.rmin());
    const double kx = mode * 2. * M_PI / Lx;
    for (MCoordX ix : gridx) {
        const RCoordX x = gridx.to_real(ix);
        perturbation(ix) = perturb_amplitude * cos(kx * x);
    }
}

/*
 Distribution function initialization :
  f(x,v) = fM*(1+perturb*cos(kx*x))
      with kx = 2*pi*mode/Lx
   where fM is a Maxwellian :
 The boundary conditions are  :
    f(x=x0,v) = f(x=xmax,v) = fM(v)
*/
void DistributionFunction::init()
{
    MDomainX gridx = values.domain<MeshX>();
    MDomainVx gridvx = values.domain<MeshVx>();

    // Initialization of the perturbation
    const int mode = 1;
    const double perturb_amplitude = 0.001;
    DBlockX perturbation(gridx);
    perturbation_initialization(mode, perturb_amplitude, perturbation);

    // Initialization of the Maxwellian --> fill Maxw_values
    Maxwellian_initialization(density_eq, temperature_eq, mean_velocity_eq, Maxw_values);

    // Initialization of the distribution function --> fill values
    for (MCoordX ix : gridx) {
        for (MCoordVx iv : gridvx) {
            double fdistribu_val = Maxw_values(iv) * (1. + perturbation(ix));
            if (fdistribu_val < 1.e-60) {
                fdistribu_val = 1.e-60;
            }
            values(ix, iv) = fdistribu_val;
        }
    }
}
