// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>

#include <geometry.hpp>

#include "fftpoissonsolver.hpp"

FftPoissonSolver::FftPoissonSolver(IChargeDensityCalculator const& compute_rho)
    : m_compute_rho(compute_rho)
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
    DFieldXY rho(xy_dom);
    DFieldVxVy contiguous_slice_vxvy(allfdistribu.domain<IDimVx, IDimVy>());
    m_compute_rho(rho, allfdistribu);

    // Build a mesh in the fourier space, for N points
    IDomainFxFy const k_mesh = ddc::FourierMesh<IDimFx, IDimFy>(xy_dom, false);

    device_t<ddc::Chunk<Kokkos::complex<double>, IDomainFxFy>> intermediate_chunk_alloc(k_mesh);
    ddc::ChunkSpan intermediate_chunk = intermediate_chunk_alloc.span_view();

    // Compute FFT(rho)
    ddc::FFT_Normalization norm = ddc::FFT_Normalization::BACKWARD;
    ddc::
            fft(Kokkos::DefaultExecutionSpace(),
                intermediate_chunk.span_view(),
                rho.span_view(),
                ddc::kwArgs_fft {norm});

    // Solve Poisson's equation -d2Phi/dx2 = rho
    //   in Fourier space as -kx*kx*FFT(Phi)=FFT(rho))
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            k_mesh,
            KOKKOS_LAMBDA(IndexFxFy const ikxky) {
                IndexFx const ikx = ddc::select<IDimFx>(ikxky);
                IndexFy const iky = ddc::select<IDimFy>(ikxky);

                if (ikxky != k_mesh.front()) {
                    intermediate_chunk(ikxky)
                            = intermediate_chunk(ikxky)
                              / ((double)coordinate(ikx) * (double)coordinate(ikx)
                                 + (double)coordinate(iky) * (double)coordinate(iky));
                } else {
                    intermediate_chunk(ikxky) = 0.;
                }
            });

    // Find the electric field in Fourier space
    // FFT(efield) = [-kx*i*FFT(Phi), -ky*i*FFT(Phi)]
    Kokkos::Profiling::pushRegion("ElectricField");
    Kokkos::complex<double> imaginary_unit(0.0, 1.0);
    ddc::Chunk fourier_efield_alloc
            = ddc::Chunk(k_mesh, ddc::DeviceAllocator<Kokkos::complex<double>>());
    ddc::ChunkSpan fourier_efield = fourier_efield_alloc.span_view();

    // Calculate x component
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            k_mesh,
            KOKKOS_LAMBDA(IndexFxFy const ikxky) {
                IndexFx const ikx = ddc::select<IDimFx>(ikxky);
                fourier_efield(ikxky)
                        = -imaginary_unit * coordinate(ikx) * intermediate_chunk(ikxky);
            });

    // Perform the inverse 1D FFT of the solution to deduce the electric field
    ddc::
            ifft(Kokkos::DefaultExecutionSpace(),
                 electric_field_x.span_view(),
                 fourier_efield.span_view(),
                 ddc::kwArgs_fft {norm});

    // Calculate y component
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            k_mesh,
            KOKKOS_LAMBDA(IndexFxFy const ikxky) {
                IndexFy const iky = ddc::select<IDimFy>(ikxky);
                fourier_efield(ikxky)
                        = -imaginary_unit * coordinate(iky) * intermediate_chunk(ikxky);
            });

    // Perform the inverse 1D FFT of the solution to deduce the electric field
    ddc::
            ifft(Kokkos::DefaultExecutionSpace(),
                 electric_field_y.span_view(),
                 fourier_efield.span_view(),
                 ddc::kwArgs_fft {norm});
    Kokkos::Profiling::popRegion();

    // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
    ddc::
            ifft(Kokkos::DefaultExecutionSpace(),
                 electrostatic_potential.span_view(),
                 intermediate_chunk.span_view(),
                 ddc::kwArgs_fft {norm});

    Kokkos::Profiling::popRegion();
}
