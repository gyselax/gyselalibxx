// SPDX-License-Identifier: MIT
#include <tuple>
#include <type_traits>
#include <utility>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <crank_nicolson.hpp>

using namespace ddc;


struct RDimX
{
    bool PERIODIC = false;
};

struct IDimX : UniformPointSampling<RDimX>
{
};

TEST(CrankNicolsonFixture, CrankNicolsonOrder)
{
    using CoordX = Coordinate<RDimX>;
    using IndexX = DiscreteElement<IDimX>;
    using IVectX = DiscreteVector<IDimX>;
    using IDomainX = DiscreteDomain<IDimX>;
    using DChunkX = Chunk<double, IDomainX>;

    CoordX x_min(0.0);
    CoordX x_max(1.0);
    IVectX x_size(10);

    IndexX start(0);

    int constexpr Ntests = 5;

    double dt(0.01);
    int Nt(1);

    std::array<double, Ntests> error;
    std::array<double, Ntests - 1> order;

    init_discrete_space<IDimX>(IDimX::init(x_min, x_max, x_size));
    IDomainX dom(start, x_size);

    CrankNicolson<DChunkX> crank_nicolson(dom);

    DChunkX vals(dom);
    DChunkX result(dom);

    double exp_val = exp(5.0 * dt * Nt);
    ddc::for_each(dom, [&](IndexX ix) {
        double const C = (double(ix.uid()) - 0.6);
        result(ix) = C * exp_val + 0.6;
    });

    for (int j(0); j < Ntests; ++j) {
        ddc::for_each(dom, [&](IndexX ix) { vals(ix) = double(ix.uid()); });

        for (int i(0); i < Nt; ++i) {
            crank_nicolson
                    .update(vals,
                            dt,
                            [&](ChunkSpan<double, IDomainX> dy,
                                ChunkSpan<double const, IDomainX> y) {
                                ddc::for_each(dom, [&](IndexX ix) { dy(ix) = 5.0 * y(ix) - 3.0; });
                            });
        }

        double linf_err = 0.0;
        ddc::for_each(dom, [&](IndexX ix) {
            double const err = abs(result(ix) - vals(ix));
            linf_err = err > linf_err ? err : linf_err;
        });
        error[j] = linf_err;

        dt *= 0.5;
        Nt *= 2;
    }
    for (int j(0); j < Ntests - 1; ++j) {
        order[j] = log(error[j] / error[j + 1]) / log(2.0);
        EXPECT_NEAR(order[j], 2., 1e-1);
    }
}
