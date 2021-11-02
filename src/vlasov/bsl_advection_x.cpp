#include <cassert>

#include <ddc/discretization>

#include <species_info.hpp>

#include "bsl_advection_x.hpp"
#include "i_interpolator_x.hpp"

using namespace std;

BslAdvectionX::BslAdvectionX(
        SpeciesInformation const& species_info,
        IPreallocatableInterpolatorX const& interpolator)
    : m_interpolator(interpolator)
    , m_species_info(species_info)
{
}

DSpanSpXVx BslAdvectionX::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    IDomainX const& x_dom = get_domain<IDimX>(allfdistribu);
    IDomainVx const& v_dom = get_domain<IDimVx>(allfdistribu);
    IDomainSp const& sp_dom = get_domain<IDimSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DFieldX feet_coords(x_dom);
    DFieldX contiguous_slice(x_dom);
    InterpolatorXProxy const interpolator = m_interpolator.preallocate();

    for (IndexSp const isp : sp_dom) {
        double const sqrt_me_on_mspecies = std::sqrt(
                m_species_info.mass()(m_species_info.ielec()) / m_species_info.mass()(isp));

        for (IndexVx const iv : v_dom) {
            // compute the displacement
            double const dx = sqrt_me_on_mspecies * dt * to_real(iv);

            // compute the coordinates of the feet
            for (IndexX const ix : x_dom) {
                feet_coords(ix) = to_real(ix) - dx;
            }

            // copy the slice in contiguous memory
            deepcopy(contiguous_slice, allfdistribu[isp][iv]);

            // interpolate the function at the feet using the provided interpolator
            interpolator(contiguous_slice, feet_coords);

            // copy back
            deepcopy(allfdistribu[isp][iv], contiguous_slice);
        }
    }
    return allfdistribu;
}
