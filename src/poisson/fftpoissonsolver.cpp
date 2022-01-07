#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/Chunk>
#include <ddc/ChunkSpan>
#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/NonUniformDiscretization>
#include <ddc/UniformDiscretization>

#include <sll/bsplines_uniform.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>
#include <species_info.hpp>

#include "fftpoissonsolver.hpp"
#include "poissonsolvercommons.hpp"

FftPoissonSolver::FftPoissonSolver(
        SpeciesInformation const& species_info,
        IFourierTransform<RDimX> const& fft,
        IInverseFourierTransform<RDimX> const& ifft,
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_species_info(species_info)
    , m_fft(fft)
    , m_ifft(ifft)
    , m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_spline_vx_builder(spline_vx_builder)
    , m_spline_vx_evaluator(spline_vx_evaluator)
    , m_derivs_vxmin_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_data.data(), m_derivs_vxmin_data.size())
    , m_derivs_vxmax_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_data.data(), m_derivs_vxmax_data.size())
{
}

// 1- Inner solvers shall be passed in the constructor
// 2- Should it take an array of distribution functions ?
void FftPoissonSolver::operator()(
        DSpanX const electrostatic_potential,
        DSpanX const electric_field,
        DViewSpXVx const allfdistribu) const
{
    assert(electrostatic_potential.domain() == get_domain<IDimX>(allfdistribu));
    IDomainX const x_dom = electrostatic_potential.domain();

    // Compute the RHS of the Poisson equation.
    Chunk<double, IDomainX> rho(x_dom);
    DFieldVx contiguous_slice_vx(allfdistribu.domain<IDimVx>());
    Chunk<double, BSDomainVx> vx_spline_coef(m_spline_vx_builder.spline_domain());
    compute_rho(
            rho,
            m_species_info,
            m_spline_vx_builder,
            m_spline_vx_evaluator,
            m_derivs_vxmin,
            m_derivs_vxmax,
            allfdistribu);

    // Build a mesh in the fourier space, for N points
    IDimFx const mesh_fx = m_fft.compute_fourier_domain(x_dom);
    IDomainFx const dom_fx(DiscreteVector<IDimFx>(mesh_fx.size()));

    // Compute FFT(rho)
    Chunk<std::complex<double>, IDomainFx> complex_Phi_fx(dom_fx);
    m_fft(complex_Phi_fx, rho.span_view());

    // Solve Poisson's equation -d2Phi/dx2 = rho
    //   in Fourier space as -kx*kx*FFT(Phi)=FFT(rho))
    complex_Phi_fx(dom_fx.front()) = 0.;
    for (auto it_freq = dom_fx.cbegin() + 1; it_freq != dom_fx.cend(); ++it_freq) {
        double const kx = 2. * M_PI * mesh_fx.to_real(*it_freq);
        complex_Phi_fx(*it_freq) /= kx * kx;
    }

    // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
    m_ifft(electrostatic_potential, complex_Phi_fx);

    // Compute efield = -dPhi/dx where Phi is the electrostatic potential
    compute_electric_field_fromvalues(
            electric_field,
            m_spline_x_builder,
            m_spline_x_evaluator,
            electrostatic_potential);
}
