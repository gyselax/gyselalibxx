#pragma once

#include <sll/spline_evaluator.hpp>

#include <fft.hpp>
#include <geometry.hpp>
#include <ifft.hpp>

#include "ipoissonsolver.hpp"

class SpeciesInformation;

class FftPoissonSolver : public IPoissonSolver
{
    IFourierTransform<Dim::X> const& m_fft;

    IInverseFourierTransform<Dim::X> const& m_ifft;

    SplineVxBuilder const& m_spline_vx_builder;

    SplineEvaluator<BSplinesVx> m_spline_vx_evaluator;

    std::vector<double> m_derivs_vxmin_data;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_data;

    Span1D<double> m_derivs_vxmax;

    SpeciesInformation const& m_species_info;

public:
    FftPoissonSolver(
            SpeciesInformation const& species_info,
            IFourierTransform<Dim::X> const& fft,
            IInverseFourierTransform<Dim::X> const& ifft,
            BSplinesVx const& bsplines_vx,
            SplineVxBuilder const& spline_vx_builder);

    DSpanX operator()(DSpanX electric_potential, DViewSpXVx allfdistribu) const override;
};
