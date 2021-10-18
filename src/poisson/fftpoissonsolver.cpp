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

using namespace std::complex_literals;

FftPoissonSolver::FftPoissonSolver(
        SpeciesInformation const& species_info,
        IFourierTransform<RDimX> const& fft,
        IInverseFourierTransform<RDimX> const& ifft,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_fft(fft)
    , m_ifft(ifft)
    , m_spline_vx_builder(spline_vx_builder)
    , m_spline_vx_evaluator(spline_vx_evaluator)
    , m_derivs_vxmin_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_data.data(), m_derivs_vxmin_data.size())
    , m_derivs_vxmax_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_data.data(), m_derivs_vxmax_data.size())
    , m_species_info(species_info)
{
}

// 1- Inner solvers sall be passed in the constructor
// 2- Should it take an array of distribution functions ?
DSpanX FftPoissonSolver::operator()(DSpanX electric_potential, DViewSpXVx allfdistribu) const
{
    assert(electric_potential.domain() == get_domain<IDimX>(allfdistribu));
    IDomainX dom_x = electric_potential.domain();

    // Compute the RHS of the Poisson equation.
    Chunk<double, IDomainX> rho(dom_x);
    DFieldVx contiguous_slice_vx(allfdistribu.domain<IDimVx>());
    Chunk<double, BSDomainVx> vx_spline_coef(m_spline_vx_builder.spline_domain());
    for (IndexX ix : rho.domain()) {
        rho(ix) = m_species_info.charge()(m_species_info.ielec());
        for (IndexSp isp : get_domain<IDimSp>(allfdistribu)) {
            deepcopy(contiguous_slice_vx, allfdistribu[isp][ix]);
            m_spline_vx_builder(
                    vx_spline_coef.span_view(),
                    contiguous_slice_vx.span_cview(),
                    &m_derivs_vxmin,
                    &m_derivs_vxmax);
            rho(ix) += m_species_info.charge()(isp)
                       * m_spline_vx_evaluator.integrate(vx_spline_coef.span_cview());
        }
    }

    // Build a mesh in the fourier space, for N points
    IDimFx const mesh_fx = m_fft.compute_fourier_domain(dom_x);
    IDomainFx const dom_fx(mesh_fx, DiscreteVector<IDimFx>(mesh_fx.size()));

    // Compute FFT(rho)
    Chunk<std::complex<double>, IDomainFx> complex_Phi_fx(dom_fx);
    m_fft(complex_Phi_fx, rho);

    // Solve Poisson's equation -d2Phi/dx2 = rho
    //   in Fourier space as -kx*kx*FFT(Phi)=FFT(rho))
    complex_Phi_fx(dom_fx.front()) = 0.;
    for (auto it_freq = dom_fx.cbegin() + 1; it_freq != dom_fx.cend(); ++it_freq) {
        double const kx = 2. * M_PI * mesh_fx.to_real(*it_freq);
        complex_Phi_fx(*it_freq) /= kx * kx;
    }

    // Perform the inverse 1D FFT of the solution.
    m_ifft(electric_potential, complex_Phi_fx);

    return electric_potential;
}
