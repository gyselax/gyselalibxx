#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include "block.h"
#include "block_span.h"
#include "block_spline.h"
#include "bsplines_uniform.h"
#include "deepcopy.h"
#include "efieldfftsolver.h"
#include "geometry.h"
#include "mcoord.h"
#include "mdomain.h"
#include "non_uniform_mesh.h"
#include "null_boundary_value.h"
#include "product_mdomain.h"
#include "product_mesh.h"
#include "rcoord.h"
#include "spline_builder.h"
#include "spline_evaluator.h"
#include "taggedvector.h"
#include "uniform_mesh.h"

static double compute_dens(DViewVx const& fdistribu)
{
    auto&& dom_vx = get_domain<MeshVx>(fdistribu);
    double const dv = get<MeshVx>(dom_vx.mesh()).step();
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

static void compute_efield(DSpanX const& ex, DViewX const& phi_x)
{
    auto&& dom_x = phi_x.domain();
    BSplinesX bsplines(dom_x.get<MeshX>());
    SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC> builder(bsplines);
    Block<BSplinesX, double> coef(bsplines);
    builder(coef, phi_x, nullptr, nullptr);

    SplineEvaluator evaluator(coef, NullBoundaryValue::value, NullBoundaryValue::value);
    for (MCoordX ix : dom_x) {
        ex(ix) = -evaluator.deriv(dom_x.to_real(ix));
    }
}

EfieldFftSolver::EfieldFftSolver(
        IFourierTransform<Dim::X> const& fft,
        IInverseFourierTransform<Dim::X> const& ifft)
    : m_fft(fft)
    , m_ifft(ifft)
{
}

// 1- Inner solvers sall be passed in the constructor
// 2- Should it take an array of distribution functions ?
DSpanX EfieldFftSolver::operator()(DSpanX ex, DViewXVx fdistribu) const
{
    assert(ex.domain() == get_domain<MeshX>(fdistribu));
    UniformMDomainX dom_x = ex.domain();

    // Compute the RHS of the Poisson equation.
    Block<MDomainX, double> rho(dom_x);
    compute_rho(rho, fdistribu);

    // Copy RHS into a complex `Block`.
    Block<MDomainX, std::complex<double>> complex_rho(dom_x);
    deepcopy(complex_rho, rho);

    // Build a mesh in the fourier space, for N points
    MeshFx const dom_fx = m_fft.compute_fourier_domain(dom_x);
    auto&& dom_fx_prod_version = ProductMesh(dom_fx);
    Block<MDomainFx, std::complex<double>> complex_phi_fx(
            ProductMDomain(dom_fx_prod_version, MCoord<MeshFx>(dom_fx.size())));
    m_fft(complex_phi_fx, complex_rho);

    // Solve Poisson's equation -d^2 Phi/dx^2 = rho in the Fourier space.
    complex_phi_fx(0) = 0.;
    for (std::size_t ifreq = 1; ifreq < dom_fx.size(); ++ifreq) {
        double const kx = 2. * M_PI * dom_fx.to_real(ifreq);
        complex_phi_fx(ifreq) /= kx * kx;
    }

    // Perform the inverse 1D FFT of the solution.
    Block<MDomainX, std::complex<double>> complex_phi_x(dom_x);
    m_ifft(complex_phi_x, complex_phi_fx);

    // Duplicate the result in the (real) phi_x 1D array.
    Block<MDomainX, double> phi_x(dom_x);
    // deepcopy(phi_x, complex_phi_x);
    for (MCoordX ix : dom_x) {
        phi_x(ix) = std::real(complex_phi_x(ix));
    }

    compute_efield(ex, phi_x);

    return ex;
}
