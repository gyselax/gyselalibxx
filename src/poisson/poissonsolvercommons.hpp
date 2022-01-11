#include <iostream>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>
#include <species_info.hpp>

//===========================================================================
// Compute rho
//===========================================================================
void compute_rho(
        DSpanX rho,
        SpeciesInformation const& species_info,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator,
        Span1D<double> derivs_vxmin,
        Span1D<double> derivs_vxmax,
        DViewSpXVx allfdistribu)
{
    DFieldVx f_vx_slice(allfdistribu.domain<IDimVx>());
    Chunk<double, BSDomainVx> vx_spline_coef(spline_vx_builder.spline_domain());

    double chargedens_adiabspecies = 0.;
    IndexSp last_kin_species = allfdistribu.domain<IDimSp>().back();
    IndexSp last_species = species_info.charge().domain().back();
    if (last_kin_species != last_species) {
        chargedens_adiabspecies = double(species_info.charge()(last_species));
    }
    for (IndexX ix : rho.domain()) {
        rho(ix) = chargedens_adiabspecies;
        for (IndexSp const isp : get_domain<IDimSp>(allfdistribu)) {
            deepcopy(f_vx_slice, allfdistribu[isp][ix]);
            spline_vx_builder(
                    vx_spline_coef.span_view(),
                    f_vx_slice.span_cview(),
                    &derivs_vxmin,
                    &derivs_vxmax);
            rho(ix) += species_info.charge()(isp)
                       * spline_vx_evaluator.integrate(vx_spline_coef.span_cview());
        }
    }
}


//===========================================================================
// Compute efield = -dPhi/dx where Phi is the electrostatic potential
//  input : Phi values
//===========================================================================
void compute_electric_field_fromvalues(
        DSpanX electric_field,
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        DViewX const electrostatic_potential)
{
    IDomainX const& x_dom = electrostatic_potential.domain();
    Chunk<double, BSDomainX> elecpot_spline_coef(spline_x_builder.spline_domain());
    spline_x_builder(elecpot_spline_coef, electrostatic_potential);

    for (IndexX const ix : x_dom) {
        electric_field(ix) = -spline_x_evaluator.deriv(to_real(ix), elecpot_spline_coef);
    }
}
