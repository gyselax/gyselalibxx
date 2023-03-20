// SPDX-License-Identifier: MIT

#pragma once

#include <functional>

#include <ddc/ddc.hpp>

#include <itimesolver.hpp>

class RK2_solver : public ITimeSolver
{
private:
    std::function<void(DSpanVx, DViewSpXVx, double, IndexSpX const&)> m_rhs;

public:
    RK2_solver(std::function<
               void(DSpanVx df, DViewSpXVx allfdistribu, double time, IndexSpX const& ispx)> rhs);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double deltat, int steps = 1) const override;
};
