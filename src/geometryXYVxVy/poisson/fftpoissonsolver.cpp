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
        SplineXYBuilder const& spline_xy_builder,
        SplineXYEvaluator const& spline_xy_evaluator,
        SplineVxVyBuilder const& spline_vxvy_builder,
        SplineVxVyEvaluator const& spline_vxvy_evaluator)
    : m_compute_rho(spline_vxvy_builder, spline_vxvy_evaluator)
    , m_electric_field(spline_xy_builder, spline_xy_evaluator)
{
}

// 1- Inner solvers shall be passed in the constructor
// 2- Should it take an array of distribution functions ?
void FftPoissonSolver::operator()(
        DSpanXY const electrostatic_potential,
        DSpanXY const electric_field_x,
        DSpanXY const electric_field_y,
        DViewSpXYVxVy const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("PoissonSolver");
    assert((electrostatic_potential.domain() == ddc::get_domain<IDimX, IDimY>(allfdistribu)));
    IDomainXY const xy_dom = electrostatic_potential.domain();

    // Compute the RHS of the Poisson equation.
    ddc::Chunk<double, IDomainXY> rho(xy_dom);
    DFieldVxVy contiguous_slice_vxvy(allfdistribu.domain<IDimVx, IDimVy>());
    m_compute_rho(rho, allfdistribu);

    // Build a mesh in the fourier space, for N points
    IDomainFxFy const k_mesh = ddc::FourierMesh(xy_dom, false);

    ddc::Chunk<Kokkos::complex<double>, IDomainFxFy> intermediate_chunk
            = ddc::Chunk(k_mesh, ddc::HostAllocator<Kokkos::complex<double>>());

    // Compute FFT(rho)
    ddc::FFT_Normalization norm = ddc::FFT_Normalization::BACKWARD;
    ddc::
            fft(Kokkos::DefaultHostExecutionSpace(),
                intermediate_chunk.span_view(),
                rho.span_view(),
                (ddc::kwArgs_fft) {.normalization = norm});

    // Solve Poisson's equation -d2Phi/dx2 = rho
    //   in Fourier space as -kx*kx*FFT(Phi)=FFT(rho))
    intermediate_chunk(k_mesh.front()) = 0.;
    ddc::for_each(k_mesh, [&](IndexFxFy const ikxky) {
        IndexFx const ikx = ddc::select<IDimFx>(ikxky);
        IndexFy const iky = ddc::select<IDimFy>(ikxky);

        if (ikxky != k_mesh.front()) {
            intermediate_chunk(ikxky) = intermediate_chunk(ikxky)
                                        / ((double)coordinate(ikx) * (double)coordinate(ikx)
                                           + (double)coordinate(iky) * (double)coordinate(iky));
        }
    });

    // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
    ddc::
            ifft(Kokkos::DefaultHostExecutionSpace(),
                 electrostatic_potential.span_view(),
                 intermediate_chunk.span_view(),
                 (ddc::kwArgs_fft) {.normalization = norm});

    // Compute efield = -dPhi/dx where Phi is the electrostatic potential
    m_electric_field(electric_field_x, electric_field_y, electrostatic_potential);
    Kokkos::Profiling::popRegion();
}
