// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include <geometry.hpp>

#include "fftqnsolver.hpp"

FftQNSolver::FftQNSolver(IChargeDensityCalculator const& compute_rho) : m_compute_rho(compute_rho)
{
}

// 1- Inner solvers shall be passed in the constructor
// 2- Should it take an array of distribution functions ?
void FftQNSolver::operator()(
        DSpanX const electrostatic_potential,
        DSpanX const electric_field,
        DViewSpXVx const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("QNSolver");
    assert(electrostatic_potential.domain() == ddc::get_domain<IDimX>(allfdistribu));
    IDomainX const x_dom = electrostatic_potential.domain();
    // Compute the RHS of the Quasi-Neutrality equation.
    DFieldX rho(x_dom);

    m_compute_rho(rho, allfdistribu);

    // Build a mesh in the fourier space, for N points
    IDomainFx const k_mesh = ddc::FourierMesh<IDimFx>(x_dom, false);
    device_t<ddc::Chunk<Kokkos::complex<double>, IDomainFx>> intermediate_chunk_alloc(k_mesh);
    ddc::ChunkSpan intermediate_chunk = intermediate_chunk_alloc.span_view();
    // Compute FFT(rho)
    ddc::FFT_Normalization norm = ddc::FFT_Normalization::BACKWARD;
    ddc::
            fft(Kokkos::DefaultExecutionSpace(),
                intermediate_chunk,
                rho.span_view(),
                ddc::kwArgs_fft {norm});

    // Solve Poisson's equation -d2Phi/dx2 = rho
    //   in Fourier space as kx*kx*FFT(Phi)=FFT(rho)

    // First, 0 mode of Phi is set to 0 to avoid divergency. Note: to allow writing in the GPU memory, starting a device kernel is necessary which is performed by iterating on a 0D domain (single element).
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            ddc::DiscreteDomain<>(),
            KOKKOS_LAMBDA(ddc::DiscreteElement<>) {
                intermediate_chunk(k_mesh.front()) = Kokkos::complex<double>(0.);
            });

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            k_mesh.remove_first(IVectFx(1)),
            KOKKOS_LAMBDA(IndexFx const ikx) {
                intermediate_chunk(ikx)
                        = intermediate_chunk(ikx) / (coordinate(ikx) * coordinate(ikx));
            });

    // Find the electric field in Fourier space
    // FFT(efield) = -kx*i*FFT(Phi)
    Kokkos::Profiling::pushRegion("ElectricField");
    Kokkos::complex<double> imaginary_unit(0.0, 1.0);
    device_t<ddc::Chunk<Kokkos::complex<double>, IDomainFx>> fourier_efield_alloc(k_mesh);
    ddc::ChunkSpan fourier_efield = fourier_efield_alloc.span_view();

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            k_mesh,
            KOKKOS_LAMBDA(IndexFx const ikx) {
                fourier_efield(ikx) = -imaginary_unit * coordinate(ikx) * intermediate_chunk(ikx);
            });

    // Perform the inverse 1D FFT of the solution to deduce the electric field
    ddc::
            ifft(Kokkos::DefaultExecutionSpace(),
                 electric_field,
                 fourier_efield,
                 ddc::kwArgs_fft {norm});
    Kokkos::Profiling::popRegion();

    // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
    ddc::
            ifft(Kokkos::DefaultExecutionSpace(),
                 electrostatic_potential.span_view(),
                 intermediate_chunk,
                 ddc::kwArgs_fft {norm});
    Kokkos::Profiling::popRegion();
}
