// SPDX-License-Identifier: MIT
#include <tuple>
#include <type_traits>
#include <utility>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "crank_nicolson.hpp"
#include "directional_tag.hpp"
#include "utils_tools.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

using namespace ddc;


struct X
{
    bool PERIODIC = false;
};

struct Y
{
    bool PERIODIC = false;
};

struct GridX : UniformGridBase<X>
{
};

struct GridY : UniformGridBase<Y>
{
};

TEST(CrankNicolson2DFixtureMixedTypes, CrankNicolson2DOrderMixedTypes)
{
    using CoordX = Coord<X>;
    using CoordY = Coord<Y>;
    using CoordXY = Coord<X, Y>;
    using IdxX = Idx<GridX>;
    using IdxY = Idx<GridY>;
    using IdxStepX = IdxStep<GridX>;
    using IdxStepY = IdxStep<GridY>;
    using IdxRangeX = IdxRange<GridX>;
    using IdxRangeY = IdxRange<GridY>;
    using IdxXY = Idx<GridX, GridY>;
    using IdxRangeXY = IdxRange<GridX, GridY>;
    using AdvectionFieldMem = host_t<VectorFieldMem<double, IdxRangeXY, NDTag<X, Y>>>;
    using CFieldMemXY = host_t<FieldMem<CoordXY, IdxRangeXY>>;
    using Method = CrankNicolson<CFieldMemXY, AdvectionFieldMem, Kokkos::DefaultHostExecutionSpace>;
    using AdvectionField = host_t<VectorField<double, IdxRangeXY, NDTag<X, Y>>>;
    using ConstAdvectionField = host_t<VectorConstField<double, IdxRangeXY, NDTag<X, Y>>>;

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

    Method crank_nicolson(idx_range);

    CFieldMemXY vals(idx_range);
    host_t<CFieldMemXY> result(idx_range);

    double cos_val = std::cos(omega * dt * Nt);
    double sin_val = std::sin(omega * dt * Nt);
    ddc::for_each(idx_range, [&](IdxXY ixy) {
        double const dist_x = (coordinate(select<GridX>(ixy)) - xc);
        double const dist_y = (coordinate(select<GridY>(ixy)) - yc);

        ddc::get<X>(result(ixy)) = xc + dist_x * cos_val - dist_y * sin_val;
        ddc::get<Y>(result(ixy)) = yc + dist_x * sin_val + dist_y * cos_val;
    });

    for (int j(0); j < Ntests; ++j) {
        ddc::for_each(idx_range, [&](IdxXY ixy) { vals(ixy) = coordinate(ixy); });

        for (int i(0); i < Nt; ++i) {
            crank_nicolson.update(
                    Kokkos::DefaultHostExecutionSpace(),
                    vals,
                    dt,
                    [yc,
                     xc,
                     &idx_range,
                     omega](AdvectionField dy, host_t<ConstField<CoordXY, IdxRangeXY>> y) {
                        ddc::for_each(idx_range, [&](IdxXY ixy) {
                            ddcHelper::get<X>(dy)(ixy) = omega * (yc - ddc::get<Y>(y(ixy)));
                            ddcHelper::get<Y>(dy)(ixy) = omega * (ddc::get<X>(y(ixy)) - xc);
                        });
                    },
                    [&idx_range](
                            host_t<Field<CoordXY, IdxRangeXY>> y,
                            ConstAdvectionField dy,
                            double dt) {
                        ddc::for_each(idx_range, [&](IdxXY ixy) { y(ixy) += dt * dy(ixy); });
                    });
        }

        double linf_err = 0.0;
        ddc::for_each(idx_range, [&](IdxXY ixy) {
            double const err_x = ddc::get<X>(result(ixy) - vals(ixy));
            double const err_y = ddc::get<Y>(result(ixy) - vals(ixy));
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

void CrankNicolson2DOrderMixedTypesTest()
{
    using CoordX = Coord<X>;
    using CoordY = Coord<Y>;
    using CoordXY = Coord<X, Y>;
    using IdxX = Idx<GridX>;
    using IdxY = Idx<GridY>;
    using IdxStepX = IdxStep<GridX>;
    using IdxStepY = IdxStep<GridY>;
    using IdxRangeX = IdxRange<GridX>;
    using IdxRangeY = IdxRange<GridY>;
    using IdxXY = Idx<GridX, GridY>;
    using IdxRangeXY = IdxRange<GridX, GridY>;
    using AdvectionFieldMem = VectorFieldMem<double, IdxRangeXY, NDTag<X, Y>>;
    using CFieldMemXY = FieldMem<CoordXY, IdxRangeXY>;
    using CFieldXY = Field<CoordXY, IdxRangeXY>;
    using Method = CrankNicolson<CFieldMemXY, AdvectionFieldMem>;
    using AdvectionField = VectorField<double, IdxRangeXY, NDTag<X, Y>>;
    using ConstAdvectionField = VectorConstField<double, IdxRangeXY, NDTag<X, Y>>;

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

    Method crank_nicolson(idx_range);

    CFieldMemXY vals_alloc(idx_range);
    host_t<CFieldMemXY> result(idx_range);
    CFieldXY vals = get_field(vals_alloc);

    double cos_val = std::cos(omega * dt * Nt);
    double sin_val = std::sin(omega * dt * Nt);
    ddc::for_each(idx_range, [&](IdxXY ixy) {
        double const dist_x = (coordinate(select<GridX>(ixy)) - xc);
        double const dist_y = (coordinate(select<GridY>(ixy)) - yc);

        ddc::get<X>(result(ixy)) = xc + dist_x * cos_val - dist_y * sin_val;
        ddc::get<Y>(result(ixy)) = yc + dist_x * sin_val + dist_y * cos_val;
    });

    for (int j(0); j < Ntests; ++j) {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY ixy) { vals(ixy) = coordinate(ixy); });

        for (int i(0); i < Nt; ++i) {
            crank_nicolson.update(
                    Kokkos::DefaultExecutionSpace(),
                    vals,
                    dt,
                    [yc,
                     xc,
                     &idx_range,
                     omega](AdvectionField dy, ConstField<CoordXY, IdxRangeXY> y) {
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range,
                                KOKKOS_LAMBDA(IdxXY ixy) {
                                    ddcHelper::get<X>(dy)(ixy) = omega * (yc - ddc::get<Y>(y(ixy)));
                                    ddcHelper::get<Y>(dy)(ixy) = omega * (ddc::get<X>(y(ixy)) - xc);
                                });
                    },
                    [&idx_range](Field<CoordXY, IdxRangeXY> y, ConstAdvectionField dy, double dt) {
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range,
                                KOKKOS_LAMBDA(IdxXY ixy) { y(ixy) += dt * dy(ixy); });
                    });
        }

        auto vals_host = ddc::create_mirror_view_and_copy(vals);

        double linf_err = 0.0;
        ddc::for_each(idx_range, [&](IdxXY ixy) {
            double const err_x = ddc::get<X>(result(ixy) - vals_host(ixy));
            double const err_y = ddc::get<Y>(result(ixy) - vals_host(ixy));
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

TEST(CrankNicolson2DFixtureMixedTypes, CrankNicolson2DOrderMixedTypesGPU)
{
    CrankNicolson2DOrderMixedTypesTest();
}
