#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <irighthandside.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>


class CollisionsInter : public IRightHandSide
{
private:
    double m_nustar0;
    DFieldSpX m_nustar_profile;

public:
    CollisionsInter(IDomainSpXVx const& mesh, double nustar0);

    ~CollisionsInter() = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;

    double get_nustar0() const;

private:
    void compute_rhs(DSpanSpVx rhs, DViewSp nustar_profile, DViewSpVx allfdistribu) const;
};
