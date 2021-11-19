#pragma once

#include <sll/spline_evaluator.hpp>

#include <fft.hpp>
#include <geometry.hpp>
#include <ifft.hpp>

#include "ipoissonsolver.hpp"

class SpeciesInformation;

class FftPoissonSolver : public IPoissonSolver
{
    IFourierTransform<RDimX> const& m_fft;

    IInverseFourierTransform<RDimX> const& m_ifft;

    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> const& m_spline_x_evaluator;

    SplineVxBuilder const& m_spline_vx_builder;

    SplineEvaluator<BSplinesVx> const& m_spline_vx_evaluator;

    std::vector<double> m_derivs_vxmin_data;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_data;

    Span1D<double> m_derivs_vxmax;

    SpeciesInformation const& m_species_info;

public:
    FftPoissonSolver(
            SpeciesInformation const& species_info,
            IFourierTransform<RDimX> const& fft,
            IInverseFourierTransform<RDimX> const& ifft,
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    ~FftPoissonSolver() override = default;

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;
};
