#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/BlockSpan>
#include <ddc/ProductMDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>

#include <sll/null_boundary_value.h>
#include <sll/spline_evaluator.h>

#include <species_info.hpp>

#include "splineadvectionx.h"

class BoundaryValue;

using namespace std;
using namespace std::experimental;

SplineAdvectionX::SplineAdvectionX(
        SpeciesInformation const& species_info,
        const BSplinesX& bspl,
        const SplineXBuilder& spl_interp)
    : SplineAdvectionX(
            species_info,
            bspl,
            spl_interp,
            NullBoundaryValue::value,
            NullBoundaryValue::value)
{
}

SplineAdvectionX::SplineAdvectionX(
        SpeciesInformation const& species_info,
        const BSplinesX& bspl,
        const SplineXBuilder& spl_interp,
        const BoundaryValue& bc_left,
        const BoundaryValue& bc_right)
    : m_x_spline_basis(bspl)
    , m_spline_x_builder(spl_interp)
    , m_spline_x_evaluator(bspl, bc_left, bc_right)
    , m_species_info(species_info)
{
    assert(bspl.is_periodic());
}

DSpanSpXVx SplineAdvectionX::operator()(DSpanSpXVx fdistribu, double dt) const
{
    assert(get_domain<MeshX>(fdistribu) == m_spline_x_builder.interpolation_domain());

    const MDomainX& x_dom = get_domain<MeshX>(fdistribu);
    const MDomainVx& v_dom = get_domain<MeshVx>(fdistribu);
    const MDomainSp& sp_dom = get_domain<MeshSp>(fdistribu);

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
            deepcopy(contiguous_slice, fdistribu[isp][iv]);

            // build a spline representation of the data
            m_spline_x_builder(spline_coef, contiguous_slice);

            // evaluate the function at the feet using the spline
            m_spline_x_evaluator(contiguous_slice.view(), feet_coords.cview(), spline_coef.cview());

            // copy back
            deepcopy(fdistribu[isp][iv], contiguous_slice);
        }
    }
    return fdistribu;
}
