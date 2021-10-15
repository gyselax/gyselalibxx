#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/BlockSpan>
#include <ddc/ProductMDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>

#include <sll/null_boundary_value.hpp>
#include <sll/spline_evaluator.hpp>

#include <species_info.hpp>

#include "splineadvectionx.hpp"

class BoundaryValue;

using namespace std;
using namespace std::experimental;

SplineAdvectionX::SplineAdvectionX(
        SpeciesInformation const& species_info,
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_species_info(species_info)
{
}

DSpanSpXVx SplineAdvectionX::operator()(DSpanSpXVx allfdistribu, double dt) const
{
    assert(get_domain<MeshX>(allfdistribu) == m_spline_x_builder.interpolation_domain());

    const MDomainX& x_dom = get_domain<MeshX>(allfdistribu);
    const MDomainVx& v_dom = get_domain<MeshVx>(allfdistribu);
    const MDomainSp& sp_dom = get_domain<MeshSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DBlockX feet_coords(x_dom);
    DBlockX contiguous_slice(x_dom);

    // Construct a domain over the bounded basis and allocate memory on this support
    Block<double, BSDomainX> spline_coef(m_spline_x_builder.spline_domain());

    for (MCoordSp isp : sp_dom) {
        const double sqrt_me_on_mspecies = std::sqrt(
                m_species_info.mass()(m_species_info.ielec()) / m_species_info.mass()(isp));

        for (MCoordVx iv : v_dom) {
            // compute the displacement
            const double dx = sqrt_me_on_mspecies * dt * v_dom.to_real(iv);

            // compute the coordinates of the feet
            for (MCoordX ix : x_dom) {
                feet_coords(ix) = x_dom.to_real(ix) - dx;
            }

            // copy the slice in contiguous memory
            deepcopy(contiguous_slice, allfdistribu[isp][iv]);

            // build a spline representation of the data
            m_spline_x_builder(spline_coef, contiguous_slice);

            // evaluate the function at the feet using the spline
            m_spline_x_evaluator(contiguous_slice.view(), feet_coords.cview(), spline_coef.cview());

            // copy back
            deepcopy(allfdistribu[isp][iv], contiguous_slice);
        }
    }
    return allfdistribu;
}
