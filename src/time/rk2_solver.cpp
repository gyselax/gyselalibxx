// SPDX-License-Identifier: MIT

#include <geometry.hpp>

#include "rk2_solver.hpp"


RK2_solver::RK2_solver(
        std::function<
                void(DSpanVx df, DViewSpXVx allfdistribu, double const time, IndexSpX const ispx)>
                rhs,
        double const deltat)
    : m_rhs(rhs)
    , m_deltat(deltat)
{
}

DSpanSpXVx RK2_solver::operator()(DSpanSpXVx allfdistribu, int const steps) const
{
    DFieldSpXVx allfdistribu_half(allfdistribu.domain());
    DFieldSpXVx df(allfdistribu.domain());
    // RK2 first half step
    for_each(
            policies::parallel_host,
            get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) { m_rhs(df[ispx], allfdistribu, 0, ispx); });
    for_each(policies::parallel_host, allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        allfdistribu_half(ispxvx) = allfdistribu(ispxvx) + df(ispxvx) * m_deltat / 2.0;
    });

    // RK2 final step
    for_each(
            policies::parallel_host,
            get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) { m_rhs(df[ispx], allfdistribu_half, m_deltat / 2.0, ispx); });
    for_each(policies::parallel_host, allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        allfdistribu(ispxvx) = allfdistribu(ispxvx) + df(ispxvx) * m_deltat;
    });

    return allfdistribu;
}
