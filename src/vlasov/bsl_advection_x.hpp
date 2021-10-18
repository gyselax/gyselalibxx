#pragma once

#include <geometry.hpp>

#include "iadvectionx.hpp"

class IPreallocatableInterpolatorX;
class SpeciesInformation;

class BslAdvectionX : public IAdvectionX
{
private:
    IPreallocatableInterpolatorX const& m_interpolator;

    SpeciesInformation const& m_species_info;

public:
    BslAdvectionX(
            SpeciesInformation const& species_info,
            IPreallocatableInterpolatorX const& interpolator);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
