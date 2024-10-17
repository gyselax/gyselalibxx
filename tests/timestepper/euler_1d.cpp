// SPDX-License-Identifier: MIT
#include <tuple>
#include <type_traits>
#include <utility>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "euler.hpp"


struct X
{
    bool PERIODIC = false;
};

struct GridX : UniformGridBase<X>
{
};

TEST(EulerFixture, EulerOrder)
{
    using CoordX = Coord<X>;
    using IdxX = Idx<GridX>;
    using IdxStepX = IdxStep<GridX>;
    using IdxRangeX = IdxRange<GridX>;
    using DFieldMemX = host_t<DFieldMem<IdxRangeX>>;

    CoordX x_min(0.0);
    CoordX x_max(1.0);
    IdxStepX x_size(5);

    IdxX start(0);

    int constexpr Ntests = 5;

    double dt(0.01);
    int Nt(1);

    std::array<double, Ntests> error;
    std::array<double, Ntests - 1> order;

    ddc::init_discrete_space<GridX>(GridX::init(x_min, x_max, x_size));
    IdxRangeX idx_range(start, x_size);

    Euler<DFieldMemX, DFieldMemX, Kokkos::DefaultHostExecutionSpace> euler(idx_range);

    DFieldMemX vals(idx_range);
    DFieldMemX result(idx_range);

    double exp_val = exp(5.0 * dt * Nt);
    ddc::for_each(idx_range, [&](IdxX ix) {
        double const C = (double(ix - IdxX(0)) - 0.6);
        result(ix) = C * exp_val + 0.6;
    });

    for (int j(0); j < Ntests; ++j) {
        ddc::for_each(idx_range, [&](IdxX ix) { vals(ix) = double(ix - IdxX(0)); });

        for (int i(0); i < Nt; ++i) {
            euler
                    .update(vals,
                            dt,
                            [&](host_t<DField<IdxRangeX>> dy, host_t<DConstField<IdxRangeX>> y) {
                                ddc::for_each(idx_range, [&](IdxX ix) {
                                    dy(ix) = 5.0 * y(ix) - 3.0;
                                });
                            });
        }

        double linf_err = 0.0;
        ddc::for_each(idx_range, [&](IdxX ix) {
            double const err = abs(result(ix) - vals(ix));
            linf_err = err > linf_err ? err : linf_err;
        });
        error[j] = linf_err;

        dt *= 0.5;
        Nt *= 2;
    }
    for (int j(0); j < Ntests - 1; ++j) {
        order[j] = log(error[j] / error[j + 1]) / log(2.0);
        EXPECT_NEAR(order[j], 1., 1e-1);
    }
}

void EulerOrderGPUTest()
{
    using CoordX = Coord<X>;
    using IdxX = Idx<GridX>;
    using IdxStepX = IdxStep<GridX>;
    using IdxRangeX = IdxRange<GridX>;
    using DFieldMemX = DFieldMem<IdxRangeX>;
    using DFieldX = DField<IdxRangeX>;

    CoordX x_min(0.0);
    CoordX x_max(1.0);
    IdxStepX x_size(5);

    IdxX start(0);

    int constexpr Ntests = 5;

    double dt(0.01);
    int Nt(1);

    std::array<double, Ntests> error;
    std::array<double, Ntests - 1> order;

    ddc::init_discrete_space<GridX>(GridX::init(x_min, x_max, x_size));
    IdxRangeX idx_range(start, x_size);

    Euler<DFieldMemX> euler(idx_range);

    DFieldMemX vals_alloc(idx_range);
    host_t<DFieldMemX> result(idx_range);
    DFieldX vals = get_field(vals_alloc);

    double exp_val = exp(5.0 * dt * Nt);
    ddc::for_each(idx_range, [&](IdxX ix) {
        double const C = (double(ix - IdxX(0)) - 0.6);
        result(ix) = C * exp_val + 0.6;
    });

    for (int j(0); j < Ntests; ++j) {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxX ix) { vals(ix) = double(ix - IdxX(0)); });

        for (int i(0); i < Nt; ++i) {
            euler.update(vals, dt, [&](DField<IdxRangeX> dy, DConstField<IdxRangeX> y) {
                ddc::parallel_for_each(
                        Kokkos::DefaultExecutionSpace(),
                        idx_range,
                        KOKKOS_LAMBDA(IdxX ix) { dy(ix) = 5.0 * y(ix) - 3.0; });
            });
        }

        auto vals_host = ddc::create_mirror_view_and_copy(vals);

        double linf_err = 0.0;
        ddc::for_each(idx_range, [&](IdxX ix) {
            double const err = abs(result(ix) - vals_host(ix));
            linf_err = err > linf_err ? err : linf_err;
        });
        error[j] = linf_err;

        dt *= 0.5;
        Nt *= 2;
    }
    for (int j(0); j < Ntests - 1; ++j) {
        order[j] = log(error[j] / error[j + 1]) / log(2.0);
        EXPECT_NEAR(order[j], 1., 1e-1);
    }
}

TEST(EulerFixture, EulerOrderGPU)
{
    EulerOrderGPUTest();
}
