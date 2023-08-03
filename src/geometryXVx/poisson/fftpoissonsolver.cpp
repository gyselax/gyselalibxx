// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include <sll/bsplines_uniform.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "fftpoissonsolver.hpp"

FftPoissonSolver::FftPoissonSolver(
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_compute_rho(spline_vx_builder, spline_vx_evaluator)
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
    assert(electrostatic_potential.domain() == ddc::get_domain<IDimX>(allfdistribu));
    IDomainX const x_dom = electrostatic_potential.domain();

    // Compute the RHS of the Poisson equation.
    ddc::Chunk<double, IDomainX> rho(x_dom);
    DFieldVx contiguous_slice_vx(allfdistribu.domain<IDimVx>());
    m_compute_rho(rho, allfdistribu);

    // Build a mesh in the fourier space, for N points
    IDomainFx const k_mesh = ddc::FourierMesh(x_dom, false);

    ddc::Chunk<Kokkos::complex<double>, IDomainFx> intermediate_chunk
            = ddc::Chunk(k_mesh, ddc::HostAllocator<Kokkos::complex<double>>());

    // Compute FFT(rho)
    ddc::FFT_Normalization norm = ddc::FFT_Normalization::BACKWARD;
    ddc::
            fft(Kokkos::DefaultHostExecutionSpace(),
                intermediate_chunk.span_view(),
                rho.span_view(),
                ddc::kwArgs_fft {norm});

    // Solve Poisson's equation -d2Phi/dx2 = rho
    //   in Fourier space as -kx*kx*FFT(Phi)=FFT(rho))
    intermediate_chunk(k_mesh.front()) = 0.;
    ddc::for_each(k_mesh.remove_first(IVectFx(1)), [&](IndexFx const ikx) {
        intermediate_chunk(ikx) = intermediate_chunk(ikx) / (coordinate(ikx) * coordinate(ikx));
    });

    // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
    ddc::
            ifft(Kokkos::DefaultHostExecutionSpace(),
                 electrostatic_potential.span_view(),
                 intermediate_chunk.span_view(),
                 ddc::kwArgs_fft {norm});

    // Compute efield = -dPhi/dx where Phi is the electrostatic potential
    m_electric_field(electric_field, electrostatic_potential);
}
