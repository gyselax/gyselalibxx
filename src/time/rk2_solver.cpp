// SPDX-License-Identifier: MIT

#include <geometry.hpp>

#include "rk2_solver.hpp"


RK2_solver::RK2_solver(std::function<
                       void(DSpanVx const df,
                            DViewSpXVx const allfdistribu,
                            double const time,
                            IndexSpX const& ispx)> rhs)
    : m_rhs(std::move(rhs))
{
}

DSpanSpXVx RK2_solver::operator()(
        DSpanSpXVx const allfdistribu,
        double const deltat,
        int const steps) const
{
    DFieldSpXVx allfdistribu_half(allfdistribu.domain());
    DFieldSpXVx df(allfdistribu.domain());
    for (int iter = 0; iter < steps; ++iter) {
        // RK2 first half step
        for_each(
                policies::parallel_host,
                get_domain<IDimSp, IDimX>(allfdistribu),
                [&](IndexSpX const ispx) { m_rhs(df[ispx], allfdistribu, 0, ispx); });

        for_each(policies::parallel_host, allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
            allfdistribu_half(ispxvx) = allfdistribu(ispxvx) + df(ispxvx) * deltat / 2.0;
        });

        // RK2 final step
        for_each(
                policies::parallel_host,
                get_domain<IDimSp, IDimX>(allfdistribu),
                [&](IndexSpX const ispx) {
                    m_rhs(df[ispx], allfdistribu_half, deltat / 2.0, ispx);
                });

        for_each(policies::parallel_host, allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
            allfdistribu(ispxvx) = allfdistribu(ispxvx) + df(ispxvx) * deltat;
        });
    }

    return allfdistribu;
}
