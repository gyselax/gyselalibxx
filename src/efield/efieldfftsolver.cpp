#include <iostream>

#include "experimental/block_spline.h"
#include "experimental/blockview_spline.h"
#include "experimental/spline_builder.hpp"
#include "experimental/spline_evaluator.hpp"

#include "block.h"
#include "blockview.h"
#include "efieldfftsolver.h"
#include "fftw.h"
#include "geometry.h"
#include "null_boundary_value.h"

constexpr bool two_species = false;
constexpr bool froz_e = false;

template <bool frozen>
static double compute_dens(DBlockViewVx const& fdistribu)
{
    if constexpr (frozen) {
        auto&& dom_vx = get_domain<Dim::Vx>(fdistribu);
        double const dv = dom_vx.mesh().step();
        double rho = 0.;
        for (MCoordVx ivx : dom_vx) {
            rho += fdistribu(ivx) * dv;
        }
        return rho;
    } else {
        return 1.;
    }
}

static void compute_rho(DBlockSpanX const& rho, DBlockViewXVx const& fdistribu)
{
    for (MCoordX ix : rho.domain()) {
        double const dens_elec = compute_dens<froz_e>(fdistribu[ix]);
        double const dens_ion = 1.;
        rho(ix) = (dens_ion - dens_elec);
    }
}

static void compute_phi(DBlockSpanX const& ex, DBlockViewX const& phi)
{
    using namespace experimental;
    // TODO: solve the conflict of meshes, the problem comes from the SplineBuilder
    BSplinesX bsplines(phi.domain());
    SplineBuilder builder(bsplines, BoundCond::PERIODIC, BoundCond::PERIODIC);
    Block<BSplinesX, double> coef(bsplines);
    // builder(coef, phi_x, nullptr, nullptr);

    // SplineEvaluator evaluator(coef, NullBoundaryValue::value, NullBoundaryValue::value);
    // for (MCoordX ix : dom_x) {
    //     ex(ix) = -evaluator.deriv(dom_x.to_real(ix));
    // }
}

// 1- Inner solvers sall be passed in the constructor
// 2- Should it take an array of distribution functions ?
DBlockSpanX EfieldFftSolver::operator()(DBlockSpanX ex, DBlockViewXVx fdistribu) const
{
    assert(ex.domain() == get_domain<Dim::X>(fdistribu));
    UniformMDomainX dom_x = ex.domain();

    // Compute the RHS of the Poisson equation.
    Block<MDomainX, double> rho(dom_x);
    compute_rho(rho, fdistribu);

    // Copy RHS into a complex `Block`.
    Block<MDomainX, std::complex<double>> complex_rho(dom_x);
    deepcopy(complex_rho, rho);

    // Build a mesh in the fourier space, for N points
    // Minimal freq: 0
    // Step freq: 1/Lx
    // Maximal freq: (N-1)/Lx = 1/dx
    MDomainFx
            dom_fx(RCoordFx(0.),
                   RCoordFx(1.0 / dom_x.mesh().step()),
                   MCoordFx(0),
                   MCoordFx(dom_x.size()));
    // Perform the 1D FFT in of the RHS and store it in complex_phi_fx.
    Block<MDomainFx, std::complex<double>> complex_phi_fx(dom_fx);
    m_fft(complex_phi_fx, complex_rho);

    // Solve Poisson's equation -d^2 Phi/dx^2 = rho in the Fourier space.
    std::vector<double> frequencies = m_fft.ifftshift(dom_fx);
    complex_phi_fx(0) = 0.;
    for (std::size_t ifreq = 1; ifreq < frequencies.size(); ++ifreq) {
        double const kx = 2. * M_PI * frequencies[ifreq];
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

    compute_phi(ex, phi_x);

    return ex;
}
