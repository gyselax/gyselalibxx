// SPDX-License-Identifier: MIT
#include <tuple>
#include <type_traits>
#include <utility>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "rk2.hpp"
#include "rk3.hpp"
#include "rk4.hpp"


template <class T>
class RungeKuttaFixture;

template <std::size_t ORDER>
class RungeKuttaFixture<std::tuple<std::integral_constant<std::size_t, ORDER>>>
    : public testing::Test
{
public:
    static int constexpr order = ORDER;

    struct X
    {
        static bool constexpr PERIODIC = false;
    };
    using CoordX = Coord<X>;
    struct GridX : UniformGridBase<X>
    {
    };
    using IdxX = Idx<GridX>;
    using IdxStepX = IdxStep<GridX>;
    using IdxRangeX = IdxRange<GridX>;
    using DFieldMemX = host_t<DFieldMem<IdxRangeX>>;
    using RungeKutta = std::conditional_t<
            ORDER == 2,
            RK2<DFieldMemX>,
            std::conditional_t<ORDER == 3, RK3<DFieldMemX>, RK4<DFieldMemX>>>;
};

using runge_kutta_types = testing::Types<
        std::tuple<std::integral_constant<std::size_t, 2>>,
        std::tuple<std::integral_constant<std::size_t, 3>>,
        std::tuple<std::integral_constant<std::size_t, 4>>>;

TYPED_TEST_SUITE(RungeKuttaFixture, runge_kutta_types);

TYPED_TEST(RungeKuttaFixture, RungeKuttaOrder)
{
    using CoordX = typename TestFixture::CoordX;
    using GridX = typename TestFixture::GridX;
    using IdxX = typename TestFixture::IdxX;
    using IdxStepX = typename TestFixture::IdxStepX;
    using IdxRangeX = typename TestFixture::IdxRangeX;
    using DFieldMemX = typename TestFixture::DFieldMemX;
    using RungeKutta = typename TestFixture::RungeKutta;

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
    IdxRangeX dom(start, x_size);

    RungeKutta runge_kutta(dom);

    DFieldMemX vals(dom);
    DFieldMemX result(dom);

    double exp_val = exp(5.0 * dt * Nt);
    ddc::for_each(dom, [&](IdxX ix) {
        double const C = (double(ix.uid()) - 0.6);
        result(ix) = C * exp_val + 0.6;
    });

    for (int j(0); j < Ntests; ++j) {
        ddc::for_each(dom, [&](IdxX ix) { vals(ix) = double(ix.uid()); });

        for (int i(0); i < Nt; ++i) {
            runge_kutta
                    .update(vals,
                            dt,
                            [&](host_t<DField<IdxRangeX>> dy, host_t<DConstField<IdxRangeX>> y) {
                                ddc::for_each(dom, [&](IdxX ix) { dy(ix) = 5.0 * y(ix) - 3.0; });
                            });
        }

        double linf_err = 0.0;
        ddc::for_each(dom, [&](IdxX ix) {
            double const err = abs(result(ix) - vals(ix));
            linf_err = err > linf_err ? err : linf_err;
        });
        error[j] = linf_err;

        dt *= 0.5;
        Nt *= 2;
    }
    for (int j(0); j < Ntests - 1; ++j) {
        order[j] = log(error[j] / error[j + 1]) / log(2.0);
        EXPECT_NEAR(order[j], double(TestFixture::order), 1e-1);
    }
}
