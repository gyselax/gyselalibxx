

# File diocotron\_initialisation\_equilibrium.hpp

[**File List**](files.md) **>** [**geometryRTheta**](dir_e9f169004bcfe9f3cb1f8a27ce024e59.md) **>** [**initialisation**](dir_1b70d60e6147eeeade38d183e3e9d318.md) **>** [**diocotron\_initialisation\_equilibrium.hpp**](diocotron__initialisation__equilibrium_8hpp.md)

[Go to the documentation of this file](diocotron__initialisation__equilibrium_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

class DiocotronDensitySolution
{
private:
    double const m_W1;
    double const m_R1;
    double const m_R2;
    double const m_W2;

    double const m_Q;
    int const m_l;
    double const m_eps;

    double m_omega_Re;
    double m_omega_Im;

    const double m_r_bar;
    const double m_d;
    const double m_p;


public:
    DiocotronDensitySolution(
            double const W1,
            double const R1,
            double const R2,
            double const W2,
            double const Q,
            int const l,
            double const eps)
        : m_W1(W1)
        , m_R1(R1)
        , m_R2(R2)
        , m_W2(W2)
        , m_Q(Q)
        , m_l(l)
        , m_eps(eps)
        , m_omega_Re(0.)
        , m_omega_Im(0.)
        , m_r_bar((R1 + R2) / 2.)
        , m_d((R2 - R1) / 2.)
        , m_p(50.)
    {
        assert(0 <= m_W1);
        assert(m_W1 <= m_R1);
        assert(m_R1 < m_R2);
        assert(m_R2 <= m_W2);
        assert(m_l != 0);

        // Computation of omega: -----------------------------------------------------------------
        dispersion_relation(m_omega_Re, m_omega_Im);
    };

    double get_frequency() const;

    double get_slope() const;


    double initialisation(Coord<R, Theta> const& coord) const;

    double equilibrium(Coord<R, Theta> const& coord) const;

private:
    void dispersion_relation(double& omega_Re, double& omega_Im) const;
};
```


