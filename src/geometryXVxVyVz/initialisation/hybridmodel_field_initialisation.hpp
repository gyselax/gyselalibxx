// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>
#include "geometry.hpp"
#include "ifieldinitialisation.hpp"
#include "paraconfpp.hpp"

/// Fields for the hybrid model. This initialises magnetic field and pressure field.
class Hybridmodel_field_initialisation : public IFieldInitialisation
{
    int m_magnetic_init_perturb_mode;

    double m_magnetic_init_perturb_amplitude;

    int m_pressure_init_perturb_mode;

    double m_pressure_init_perturb_amplitude;

public:
    /**
     * @brief The constructor for the Hybridmodel_field_initialisation class.
     * @param[in] magnetic_init_perturb_mode The perturbation mode of the magnetic field
     * @param[in] magnetic_init_perturb_amplitude The perturbation amplitude of the magnetic field
     * @param[in] pressure_init_perturb_mode The perturbation mode of the pressure field
     * @param[in] pressure_init_perturb_amplitude The perturbation amplitude of the pressure field
     */
    Hybridmodel_field_initialisation(
            int magnetic_init_perturb_mode,
            double magnetic_init_perturb_amplitude,
            int pressure_init_perturb_mode,
            double pressure_init_perturb_amplitude);

    ~Hybridmodel_field_initialisation() override = default;

    /**
     * @brief Initialises magnetic_field and pressure.
     * @param[out] allfequilibrium Magnetic_field and pressure.
     * @return Magnetic_field and pressure.
     */
    DFieldX operator()(DFieldX magnetic_field_x, DFieldX magnetic_field_y, DFieldX magnetic_field_z, DFieldX pressure) const override;

};
