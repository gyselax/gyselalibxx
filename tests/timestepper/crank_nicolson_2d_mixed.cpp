// SPDX-License-Identifier: MIT
#include <tuple>
#include <type_traits>
#include <utility>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <crank_nicolson.hpp>
#include <directional_tag.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

#include "utils_tools.hpp"

using namespace ddc;


struct RDimX
{
    bool PERIODIC = false;
};

struct RDimY
{
    bool PERIODIC = false;
};

TEST(CrankNicolson2DFixtureMixedTypes, CrankNicolson2DOrderMixedTypes)
{
    using CoordX = Coordinate<RDimX>;
    using CoordY = Coordinate<RDimY>;
    using CoordXY = Coordinate<RDimX, RDimY>;
    using IDimX = UniformPointSampling<RDimX>;
    using IDimY = UniformPointSampling<RDimY>;
    using IndexX = DiscreteElement<IDimX>;
    using IndexY = DiscreteElement<IDimY>;
    using IVectX = DiscreteVector<IDimX>;
    using IVectY = DiscreteVector<IDimY>;
    using IDomainX = DiscreteDomain<IDimX>;
    using IDomainY = DiscreteDomain<IDimY>;
    using IndexXY = DiscreteElement<IDimX, IDimY>;
    using IDomainXY = DiscreteDomain<IDimX, IDimY>;
    using AdvectionField = VectorField<double, IDomainXY, NDTag<RDimX, RDimY>>;
    using CChunkXY = Chunk<CoordXY, IDomainXY>;
    using Method = CrankNicolson<CChunkXY, AdvectionField>;
    using AdvectionFieldSpan = VectorFieldSpan<double, IDomainXY, NDTag<RDimX, RDimY>>;
    using AdvectionFieldView = VectorFieldView<double, IDomainXY, NDTag<RDimX, RDimY>>;

    CoordX x_min(-1.0);
    CoordX x_max(1.0);
    IVectX x_size(5);

    CoordY y_min(-1.0);
    CoordY y_max(1.0);
    IVectY y_size(5);

    IndexX start_x(0);
    IndexY start_y(0);

    int constexpr Ntests = 2;

    double const xc = 0.25;
    double const yc = 0.0;
    double const omega = 2 * M_PI;

    double dt(0.1);
    int Nt(1);

    std::array<double, Ntests> error;
    std::array<double, Ntests - 1> order;

    init_discrete_space(IDimX::init(x_min, x_max, x_size));
    IDomainX dom_x(start_x, x_size);
    init_discrete_space(IDimY::init(y_min, y_max, y_size));
    IDomainY dom_y(start_y, y_size);

    IDomainXY dom(dom_x, dom_y);

    Method crank_nicolson(dom);

    CChunkXY vals(dom);
    CChunkXY result(dom);

    double cos_val = std::cos(omega * dt * Nt);
    double sin_val = std::sin(omega * dt * Nt);
    ddc::for_each(dom, [&](IndexXY ixy) {
        double const dist_x = (coordinate(select<IDimX>(ixy)) - xc);
        double const dist_y = (coordinate(select<IDimY>(ixy)) - yc);

        ddc::get<RDimX>(result(ixy)) = xc + dist_x * cos_val - dist_y * sin_val;
        ddc::get<RDimY>(result(ixy)) = yc + dist_x * sin_val + dist_y * cos_val;
    });

    for (int j(0); j < Ntests; ++j) {
        ddc::for_each(dom, [&](IndexXY ixy) { vals(ixy) = coordinate(ixy); });

        for (int i(0); i < Nt; ++i) {
            crank_nicolson.update(
                    Kokkos::DefaultHostExecutionSpace(),
                    vals,
                    dt,
                    [yc, xc, &dom, omega](AdvectionFieldSpan dy, ChunkView<CoordXY, IDomainXY> y) {
                        ddc::for_each(dom, [&](IndexXY ixy) {
                            ddcHelper::get<RDimX>(dy)(ixy) = omega * (yc - ddc::get<RDimY>(y(ixy)));
                            ddcHelper::get<RDimY>(dy)(ixy) = omega * (ddc::get<RDimX>(y(ixy)) - xc);
                        });
                    },
                    [&dom](ChunkSpan<CoordXY, IDomainXY> y, AdvectionFieldView dy, double dt) {
                        ddc::for_each(dom, [&](IndexXY ixy) { y(ixy) += dt * dy(ixy); });
                    });
        }

        double linf_err = 0.0;
        ddc::for_each(dom, [&](IndexXY ixy) {
            double const err_x = ddc::get<RDimX>(result(ixy) - vals(ixy));
            double const err_y = ddc::get<RDimY>(result(ixy) - vals(ixy));
            double const err = std::sqrt(err_x * err_x + err_y * err_y);
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
