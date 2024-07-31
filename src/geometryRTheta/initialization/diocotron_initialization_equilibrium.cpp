#include <cmath>

#include "diocotron_initialization_equilibrium.hpp"
#include "geometry.hpp"



void DiocotronDensitySolution::dispersion_relation(double& omega_Re, double& omega_Im) const
{
    // Rotation velocity induced by the radial electric field:
    double const omega_q = -2 * m_Q / m_R1 / m_R1; // c = 1; Bo = 1.
    // Diocotron frequency:
    double const omega_D = 1. / 2.;
    double const omega_frac = omega_q / omega_D;

    double const div_W1_R1 = m_R1 > 0 ? m_W1 / m_R1 : 1.;
    double const div_W1_R2 = m_W1 / m_R2;
    double const div_W1_W2 = m_W1 / m_W2;

    double const div_R1_R2 = m_R1 / m_R2;
    double const div_R1_W2 = m_R1 / m_W2;

    double const div_R2_W2 = m_R2 / m_W2;


    double const div_W1_R1_2l = ipow(div_W1_R1, 2 * m_l);
    double const div_W1_R2_2l = ipow(div_W1_R2, 2 * m_l);
    double const div_W1_W2_2l = ipow(div_W1_W2, 2 * m_l);

    double const div_R1_R2_2l = ipow(div_R1_R2, 2 * m_l);
    double const div_R1_W2_2l = ipow(div_R1_W2, 2 * m_l);

    double const div_R2_W2_2l = ipow(div_R2_W2, 2 * m_l);

    double const bl
            = 1. / (1 - div_W1_W2_2l) * m_l
                      * (1 - div_R1_R2 * div_R1_R2 + omega_frac * (1 + div_R1_R2 * div_R1_R2))
                      * (1 - div_W1_W2_2l)
              + (1 - div_R1_R2_2l) * (div_R2_W2_2l - div_W1_R1_2l);

    double const cl = 1. / (1 - div_W1_W2_2l) * m_l * m_l * omega_frac
                              * (1 - div_R1_R2 * div_R1_R2 + omega_frac * div_R1_R2 * div_R1_R2)
                              * (1 - div_W1_W2_2l)
                      - m_l * omega_frac * (1 - div_W1_R2_2l) * (1 - div_R2_W2_2l)
                      + m_l * (1 - div_R1_R2 * div_R1_R2 + omega_frac * div_R1_R2 * div_R1_R2)
                                * (1 - div_R1_W2_2l) * (1 - div_W1_R1_2l)
                      - (1 - div_R2_W2_2l) * (1 - div_W1_R1_2l) * (1 - div_R1_R2_2l);



    double const Delta = bl * bl - 4 * cl;
    if (Delta >= 0) {
        std::cout << "! WARNING: stable oscillations with (W1, R1, R2, W2, l) = (" << m_W1 << ", "
                  << m_R1 << ", " << m_R2 << ", " << m_W2 << ", " << m_l << ")." << std::endl;
        omega_Re = omega_D * (bl + std::sqrt(Delta)) / 2;
        omega_Im = omega_D * 0.;
    } else {
        omega_Re = omega_D * bl / 2;
        omega_Im = omega_D * std::sqrt(-Delta) / 2; // should be positive for the instability.
    }
}


double DiocotronDensitySolution::get_frequency() const
{
    return m_omega_Re;
}

double DiocotronDensitySolution::get_slope() const
{
    return m_omega_Im;
}


double DiocotronDensitySolution::initialisation(CoordRTheta const& coord) const
{
    const double r = ddc::get<R>(coord);
    const double theta = ddc::get<Theta>(coord);

    if ((m_R1 <= r) and (r <= m_R2)) {
        return (1. + m_eps * std::cos(m_l * theta))
               * std::exp(-(std::pow((r - m_r_bar) / m_d, m_p)));
    } else {
        return 0.;
    }
}

double DiocotronDensitySolution::equilibrium(CoordRTheta const& coord) const
{
    const double r = ddc::get<R>(coord);

    if ((m_R1 <= r) and (r <= m_R2)) {
        return 1. * std::exp(-(std::pow((r - m_r_bar) / m_d, m_p)));
    } else {
        return 0.;
    }
}
