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
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"

using namespace ddc;


template <class T>
class RungeKutta2DFixture;

template <std::size_t ORDER>
class RungeKutta2DFixture<std::tuple<std::integral_constant<std::size_t, ORDER>>>
    : public testing::Test
{
public:
    static int constexpr order = ORDER;

    struct X
    {
        static bool constexpr PERIODIC = false;
        static bool constexpr IS_COVARIANT = true;
        static bool constexpr IS_CONTRAVARIANT = true;
        using Dual = X;
    };

    struct Y
    {
        static bool constexpr PERIODIC = false;
        static bool constexpr IS_COVARIANT = true;
        static bool constexpr IS_CONTRAVARIANT = true;
        using Dual = Y;
    };
    using CoordX = Coordinate<X>;
    using CoordY = Coordinate<Y>;
    struct GridX : UniformGridBase<X>
    {
    };
    struct GridY : UniformGridBase<Y>
    {
    };
    using IdxX = Idx<GridX>;
    using IdxY = Idx<GridY>;
    using IdxStepX = IdxStep<GridX>;
    using IdxStepY = IdxStep<GridY>;
    using IdxRangeX = IdxRange<GridX>;
    using IdxRangeY = IdxRange<GridY>;
    using IdxXY = Idx<GridX, GridY>;
    using IdxRangeXY = IdxRange<GridX, GridY>;
    using AdvectionFieldMem = host_t<VectorFieldMem<double, IdxRangeXY, VectorIndexSet<X, Y>>>;
    using RungeKutta = std::conditional_t<
            ORDER == 2,
            RK2<AdvectionFieldMem, AdvectionFieldMem, Kokkos::DefaultHostExecutionSpace>,
            std::conditional_t<
                    ORDER == 3,
                    RK3<AdvectionFieldMem, AdvectionFieldMem, Kokkos::DefaultHostExecutionSpace>,
                    RK4<AdvectionFieldMem, AdvectionFieldMem, Kokkos::DefaultHostExecutionSpace>>>;
};

using runge_kutta_2d_types = testing::Types<
        std::tuple<std::integral_constant<std::size_t, 2>>,
        std::tuple<std::integral_constant<std::size_t, 3>>,
        std::tuple<std::integral_constant<std::size_t, 4>>>;

TYPED_TEST_SUITE(RungeKutta2DFixture, runge_kutta_2d_types);

TYPED_TEST(RungeKutta2DFixture, RungeKutta2DOrder)
{
    using X = typename TestFixture::X;
    using CoordX = typename TestFixture::CoordX;
    using GridX = typename TestFixture::GridX;
    using IdxX = typename TestFixture::IdxX;
    using IdxStepX = typename TestFixture::IdxStepX;
    using IdxRangeX = typename TestFixture::IdxRangeX;

    using Y = typename TestFixture::Y;
    using CoordY = typename TestFixture::CoordY;
    using GridY = typename TestFixture::GridY;
    using IdxY = typename TestFixture::IdxY;
    using IdxStepY = typename TestFixture::IdxStepY;
    using IdxRangeY = typename TestFixture::IdxRangeY;

    using IdxXY = typename TestFixture::IdxXY;
    using IdxRangeXY = typename TestFixture::IdxRangeXY;
    using RungeKutta = typename TestFixture::RungeKutta;
    using AdvectionFieldMem = typename TestFixture::AdvectionFieldMem;
    using AdvectionField = host_t<VectorField<double, IdxRangeXY, VectorIndexSet<X, Y>>>;
    using ConstAdvectionField = host_t<VectorConstField<double, IdxRangeXY, VectorIndexSet<X, Y>>>;

    CoordX x_min(-1.0);
    CoordX x_max(1.0);
    IdxStepX x_size(5);

    CoordY y_min(-1.0);
    CoordY y_max(1.0);
    IdxStepY y_size(5);

    IdxX start_x(0);
    IdxY start_y(0);

    int constexpr Ntests = 2;

    double const xc = 0.25;
    double const yc = 0.0;
    double const omega = 2 * M_PI;

    double dt(0.1);
    int Nt(1);

    std::array<double, Ntests> error;
    std::array<double, Ntests - 1> order;

    ddc::init_discrete_space<GridX>(GridX::init(x_min, x_max, x_size));
    IdxRangeX idx_range_x(start_x, x_size);
    ddc::init_discrete_space<GridY>(GridY::init(y_min, y_max, y_size));
    IdxRangeY idx_range_y(start_y, y_size);

    IdxRangeXY idx_range(idx_range_x, idx_range_y);

    RungeKutta runge_kutta(idx_range);

    AdvectionFieldMem vals(idx_range);
    AdvectionFieldMem result(idx_range);

    double cos_val = std::cos(omega * dt * Nt);
    double sin_val = std::sin(omega * dt * Nt);
    ddc::for_each(idx_range, [&](IdxXY ixy) {
        double const dist_x = (coordinate(select<GridX>(ixy)) - xc);
        double const dist_y = (coordinate(select<GridY>(ixy)) - yc);

        ddcHelper::get<X>(result)(ixy) = xc + dist_x * cos_val - dist_y * sin_val;
        ddcHelper::get<Y>(result)(ixy) = yc + dist_x * sin_val + dist_y * cos_val;
    });

    for (int j(0); j < Ntests; ++j) {
        ddc::for_each(idx_range, [&](IdxXY ixy) {
            ddcHelper::get<X>(vals)(ixy) = coordinate(select<GridX>(ixy));
            ddcHelper::get<Y>(vals)(ixy) = coordinate(select<GridY>(ixy));
        });

        for (int i(0); i < Nt; ++i) {
            runge_kutta.update(
                    Kokkos::DefaultHostExecutionSpace(),
                    get_field(vals),
                    dt,
                    [yc, xc, &idx_range, omega](AdvectionField dy, ConstAdvectionField y) {
                        ddc::for_each(idx_range, [&](IdxXY ixy) {
                            ddcHelper::get<X>(dy)(ixy) = omega * (yc - ddcHelper::get<Y>(y)(ixy));
                            ddcHelper::get<Y>(dy)(ixy) = omega * (ddcHelper::get<X>(y)(ixy) - xc);
                        });
                    },
                    [&idx_range](AdvectionField y, ConstAdvectionField dy, double dt) {
                        ddc::for_each(idx_range, [&](IdxXY ixy) {
                            ddcHelper::get<X>(y)(ixy) += ddcHelper::get<X>(dy)(ixy) * dt;
                            ddcHelper::get<Y>(y)(ixy) += ddcHelper::get<Y>(dy)(ixy) * dt;
                        });
                    });
        }

        double linf_err = 0.0;
        ddc::for_each(idx_range, [&](IdxXY ixy) {
            double const err_x = ddcHelper::get<X>(result)(ixy) - ddcHelper::get<X>(vals)(ixy);
            double const err_y = ddcHelper::get<Y>(result)(ixy) - ddcHelper::get<Y>(vals)(ixy);
            double const err = std::sqrt(err_x * err_x + err_y * err_y);
            linf_err = err > linf_err ? err : linf_err;
        });
        error[j] = linf_err;

        dt *= 0.5;
        Nt *= 2;
    }
    for (int j(0); j < Ntests - 1; ++j) {
        order[j] = log(error[j] / error[j + 1]) / log(2.0);
        EXPECT_NEAR(order[j], TestFixture::order, 1e-1);
    }
}
