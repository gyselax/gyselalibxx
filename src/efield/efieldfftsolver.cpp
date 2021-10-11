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

using namespace std::complex_literals;

EfieldFftSolver::EfieldFftSolver(
        IFourierTransform<Dim::X> const& fft,
        IInverseFourierTransform<Dim::X> const& ifft,
        BSplinesVx const& bsplines_vx,
        SplineVxBuilder const& spline_vx_builder)
    : m_fft(fft)
    , m_ifft(ifft)
    , m_spline_vx_basis(bsplines_vx)
    , m_spline_vx_builder(spline_vx_builder)
    , m_spline_vx_evaluator(bsplines_vx, NullBoundaryValue::value, NullBoundaryValue::value)
    , m_derivs_vxmin_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_data.data(), m_derivs_vxmin_data.size())
    , m_derivs_vxmax_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_data.data(), m_derivs_vxmax_data.size())
{
}

// 1- Inner solvers sall be passed in the constructor
// 2- Should it take an array of distribution functions ?
DSpanX EfieldFftSolver::operator()(DSpanX efield, DViewXVx fdistribu) const
{
    assert(efield.domain() == get_domain<MeshX>(fdistribu));
    UniformMDomainX dom_x = efield.domain();

    // Compute the RHS of the Poisson equation.
    Block<double, MDomainX> rho(dom_x);
    DBlockVx contiguous_slice_vx(fdistribu.domain<MeshVx>());
    Block<double, BSDomainVx> vx_spline_coef(m_spline_vx_builder.spline_domain());
    for (MCoordX ix : rho.domain()) {
        deepcopy(contiguous_slice_vx, fdistribu[ix]);
        m_spline_vx_builder(
                vx_spline_coef.view(),
                contiguous_slice_vx.cview(),
                &m_derivs_vxmin,
                &m_derivs_vxmax);
        double const dens_elec = m_spline_vx_evaluator.integrate(vx_spline_coef.cview());
        double const dens_ion = 1.;
        rho(ix) = (dens_ion - dens_elec);
    }

    // Build a mesh in the fourier space, for N points
    MeshFx const mesh_fx = m_fft.compute_fourier_domain(dom_x);
    MDomainFx const dom_fx(mesh_fx, MLength<MeshFx>(mesh_fx.size()));

    // Compute FFT(rho)
    Block<std::complex<double>, MDomainFx> complex_Ex_fx(dom_fx);
    m_fft(complex_Ex_fx, rho);

    // Solve Poisson's equation -d2Phi/dx2 = rho by solving
    //   dEx/dx = rho in the Fourier space (i.e i*kx*FFT(E)=FFT(rho)) with Ex = -dPhi/dx.
    complex_Ex_fx(dom_fx.front()) = 0.;
    for (auto it_freq = dom_fx.cbegin() + 1; it_freq != dom_fx.cend(); ++it_freq) {
        double const kx = 2. * M_PI * mesh_fx.to_real(*it_freq);
        complex_Ex_fx(*it_freq) /= 1.0i * kx;
    }

    // Perform the inverse 1D FFT of the solution.
    m_ifft(efield, complex_Ex_fx);

    return efield;
}
