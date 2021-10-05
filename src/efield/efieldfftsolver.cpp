#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/Block>
#include <ddc/BlockSpan>
#include <ddc/MCoord>
#include <ddc/NonUniformMesh>
#include <ddc/ProductMDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>
#include <ddc/UniformMesh>

#include <sll/bsplines_uniform.h>
#include <sll/null_boundary_value.h>
#include <sll/spline_builder.h>
#include <sll/spline_evaluator.h>

#include <geometry.h>

#include "efieldfftsolver.h"

static double compute_dens(DViewVx const& fdistribu)
{
    auto&& dom_vx = get_domain<MeshVx>(fdistribu);
    double const dv = dom_vx.mesh<MeshVx>().step();
    double rho = 0.;
    for (MCoordVx ivx : dom_vx) {
        rho += fdistribu(ivx) * dv;
    }
    return rho;
}

static void compute_rho(DSpanX const& rho, DViewXVx const& fdistribu)
{
    DBlockVx contiguous_slice(fdistribu.domain<MeshVx>());
    for (MCoordX ix : rho.domain()) {
        deepcopy(contiguous_slice, fdistribu[ix]);
        double const dens_elec = compute_dens(contiguous_slice);
        double const dens_ion = 1.;
        rho(ix) = (dens_ion - dens_elec);
    }
}

EfieldFftSolver::EfieldFftSolver(
        IFourierTransform<Dim::X> const& fft,
        IInverseFourierTransform<Dim::X> const& ifft,
        BSplinesX const& bsplines_x,
        SplineXBuilder const& spline_x_builder)
    : m_fft(fft)
    , m_ifft(ifft)
    , m_spline_x_basis(bsplines_x)
    , m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(bsplines_x, NullBoundaryValue::value, NullBoundaryValue::value)
{
}

// 1- Inner solvers sall be passed in the constructor
// 2- Should it take an array of distribution functions ?
DSpanX EfieldFftSolver::operator()(DSpanX ex, DViewXVx fdistribu) const
{
    assert(ex.domain() == get_domain<MeshX>(fdistribu));
    UniformMDomainX dom_x = ex.domain();

    // Compute the RHS of the Poisson equation.
    Block<double, MDomainX> rho(dom_x);
    compute_rho(rho, fdistribu);

    // Copy RHS into a complex `Block`.
    Block<std::complex<double>, MDomainX> complex_rho(dom_x);
    deepcopy(complex_rho, rho);

    // Build a mesh in the fourier space, for N points
    MeshFx const mesh_fx = m_fft.compute_fourier_domain(dom_x);
    MDomainFx const dom_fx(mesh_fx, MLength<MeshFx>(mesh_fx.size()));
    Block<std::complex<double>, MDomainFx> complex_phi_fx(dom_fx);
    m_fft(complex_phi_fx, complex_rho);

    // Solve Poisson's equation -d^2 Phi/dx^2 = rho in the Fourier space.
    complex_phi_fx(0) = 0.;
    for (std::size_t ifreq = 1; ifreq < mesh_fx.size(); ++ifreq) {
        double const kx = 2. * M_PI * mesh_fx.to_real(ifreq);
        complex_phi_fx(ifreq) /= kx * kx;
    }

    // Perform the inverse 1D FFT of the solution.
    Block<std::complex<double>, MDomainX> complex_phi_x(dom_x);
    m_ifft(complex_phi_x, complex_phi_fx);

    // Duplicate the result in the (real) phi_x 1D array.
    Block<double, MDomainX> phi_x(dom_x);
    // deepcopy(phi_x, complex_phi_x);
    for (MCoordX ix : dom_x) {
        phi_x(ix) = std::real(complex_phi_x(ix));
    }

    // Construct a domain over the bounded basis and allocate memory on this support
    BSDomainX const dom_bsx(m_spline_x_basis, MLength<BSplinesX>(m_spline_x_basis.size()));
    Block<double, BSDomainX> phi_spline_coef(dom_bsx);
    m_spline_x_builder(phi_spline_coef, phi_x);

    // Compute the electric field E = -d Phi / dx using a spline representation
    for (MCoordX ix : dom_x) {
        ex(ix) = -m_spline_x_evaluator.deriv(dom_x.to_real(ix), phi_spline_coef.cview());
    }

    return ex;
}
