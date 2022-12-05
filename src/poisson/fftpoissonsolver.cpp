// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>

#include <sll/bsplines_uniform.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "fftpoissonsolver.hpp"

FftPoissonSolver::FftPoissonSolver(
        IFourierTransform<RDimX> const& fft,
        IInverseFourierTransform<RDimX> const& ifft,
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_fft(fft)
    , m_ifft(ifft)
    , m_compute_rho(spline_vx_builder, spline_vx_evaluator)
    , m_electric_field(spline_x_builder, spline_x_evaluator)
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
    m_compute_rho(rho, allfdistribu);

    // Build a mesh in the fourier space, for N points
    IDomainFx const dom_fx(IndexFx(0), IVectFx(discrete_space<IDimFx>().size()));

    // Compute FFT(rho)
    Chunk<std::complex<double>, IDomainFx> complex_Phi_fx(dom_fx);
    m_fft(complex_Phi_fx.span_view(), rho.span_view());

    // Solve Poisson's equation -d2Phi/dx2 = rho
    //   in Fourier space as -kx*kx*FFT(Phi)=FFT(rho))
    complex_Phi_fx(dom_fx.front()) = 0.;
    for_each(dom_fx.remove_first(IVectFx(1)), [&](IndexFx const ifreq) {
        double const kx = 2. * M_PI * coordinate(ifreq);
        complex_Phi_fx(ifreq) /= kx * kx;
    });

    // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
    m_ifft(electrostatic_potential, complex_Phi_fx);

    // Compute efield = -dPhi/dx where Phi is the electrostatic potential
    m_electric_field(electric_field, electrostatic_potential);
}
