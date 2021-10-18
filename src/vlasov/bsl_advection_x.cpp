#include <cassert>

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

DSpanSpXVx BslAdvectionX::operator()(DSpanSpXVx allfdistribu, double dt) const
{
    const IDomainX& x_dom = get_domain<IDimX>(allfdistribu);
    const IDomainVx& v_dom = get_domain<IDimVx>(allfdistribu);
    const IDomainSp& sp_dom = get_domain<IDimSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DFieldX feet_coords(x_dom);
    DFieldX contiguous_slice(x_dom);
    InterpolatorXProxy interpolator = m_interpolator.preallocate();

    for (IndexSp isp : sp_dom) {
        const double sqrt_me_on_mspecies = std::sqrt(
                m_species_info.mass()(m_species_info.ielec()) / m_species_info.mass()(isp));

        for (IndexVx iv : v_dom) {
            // compute the displacement
            const double dx = sqrt_me_on_mspecies * dt * v_dom.to_real(iv);

            // compute the coordinates of the feet
            for (IndexX ix : x_dom) {
                feet_coords(ix) = x_dom.to_real(ix) - dx;
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
