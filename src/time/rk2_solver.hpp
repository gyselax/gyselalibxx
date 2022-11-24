// SPDX-License-Identifier: MIT

#pragma once

#include <functional>

#include <ddc/ddc.hpp>

#include <itimesolver.hpp>

class RK2_solver : public ITimeSolver
{
private:
    std::function<void(DSpanVx, DViewSpXVx, const double, const IndexSpX)> m_rhs;
    double const m_deltat;

public:
    RK2_solver(
            std::function<
                    void(DSpanVx df,
                         DViewSpXVx allfdistribu,
                         double const time,
                         IndexSpX const ispx)> rhs,
            double const deltat);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, int const steps) const override;
};
