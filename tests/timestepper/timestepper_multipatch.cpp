// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_uniform.hpp"
#include "crank_nicolson.hpp"
#include "euler.hpp"
#include "multipatch_field.hpp"
#include "multipatch_field_mem.hpp"
#include "multipatch_type.hpp"
#include "rk2.hpp"
#include "rk3.hpp"
#include "rk4.hpp"

namespace {

using namespace onion_shape_uniform_2d_2patches;

using DMultipatchFieldMemR = MultipatchFieldMem<DFieldMem1OnPatch, Patch1, Patch2>;
using DMultipatchFieldR = typename DMultipatchFieldMemR::span_type;
using DMultipatchConstFieldR = typename DMultipatchFieldMemR::view_type;
using MultipatchIdxRangeR = MultipatchType<IdxRange1OnPatch, Patch1, Patch2>;

using CMultipatchFieldMemRTheta = MultipatchFieldMem<CoordFieldMemOnPatch, Patch1, Patch2>;
using CMultipatchFieldRTheta = typename CMultipatchFieldMemRTheta::span_type;
using CMultipatchConstFieldRTheta = typename CMultipatchFieldMemRTheta::view_type;
using MultipatchIdxRangeRTheta = MultipatchType<IdxRangeOnPatch, Patch1, Patch2>;
using CMultipatchVectorFieldMemRTheta = MultipatchFieldMem<DVectorFieldMemOnPatch, Patch1, Patch2>;
using CMultipatchVectorFieldRTheta = typename CMultipatchVectorFieldMemRTheta::span_type;
using CMultipatchVectorConstFieldRTheta = typename CMultipatchVectorFieldMemRTheta::view_type;

using CoordRTheta = Coord<R, Theta>;
using IdxR1 = Idx<GridR<1>>;
using IdxR2 = Idx<GridR<2>>;
using IdxTheta1 = Idx<GridTheta<1>>;
using IdxTheta2 = Idx<GridTheta<2>>;
using IdxRTheta1 = Idx<GridR<1>, GridTheta<1>>;
using IdxRTheta2 = Idx<GridR<2>, GridTheta<2>>;
using DFieldR1 = DField<IdxRange<GridR<1>>>;
using DFieldR2 = DField<IdxRange<GridR<2>>>;
using DConstFieldR1 = DConstField<IdxRange<GridR<1>>>;
using DConstFieldR2 = DConstField<IdxRange<GridR<2>>>;
using CFieldRTheta1 = Field<CoordRTheta, IdxRange<GridR<1>, GridTheta<1>>>;
using CFieldRTheta2 = Field<CoordRTheta, IdxRange<GridR<2>, GridTheta<2>>>;
using CConstFieldRTheta1 = ConstField<CoordRTheta, IdxRange<GridR<1>, GridTheta<1>>>;
using CConstFieldRTheta2 = ConstField<CoordRTheta, IdxRange<GridR<2>, GridTheta<2>>>;
using VectorFieldRTheta1 = DVectorFieldOnPatch<Patch1>;
using VectorFieldRTheta2 = DVectorFieldOnPatch<Patch2>;
using VectorConstFieldRTheta1 = DVectorConstFieldOnPatch<Patch1>;
using VectorConstFieldRTheta2 = DVectorConstFieldOnPatch<Patch2>;

template <class TimeStepper>
void multipatch_timestepper_test(double expected_order)
{
    int constexpr Ntests = 5;
    std::array<double, Ntests> error;
    std::array<double, Ntests - 1> order;

    Coord<R> r_min(0.0);
    Coord<R> r_mid(1.0);
    Coord<R> r_max(2.0);

    IdxStep<GridR<1>> ncells_r_patch_1(10);
    IdxStep<GridR<2>> ncells_r_patch_2(5);

    Idx<GridR<1>> idx_start1(0);
    Idx<GridR<2>> idx_start2(0);

    double dt(0.01);
    int Nt(1);

    ddc::init_discrete_space<GridR<1>>(GridR<1>::init(r_min, r_mid, ncells_r_patch_1));
    ddc::init_discrete_space<GridR<2>>(GridR<2>::init(r_mid, r_max, ncells_r_patch_2));

    IdxRange<GridR<1>> idx_range_1(idx_start1, ncells_r_patch_1);
    IdxRange<GridR<2>> idx_range_2(idx_start2, ncells_r_patch_2);
    MultipatchIdxRangeR idx_ranges(idx_range_1, idx_range_2);

    DMultipatchFieldMemR vals_alloc(idx_ranges);
    DMultipatchFieldR vals(vals_alloc);

    TimeStepper timestepper(idx_ranges);

    host_t<DFieldMem<IdxRange<GridR<1>>>> result_patch1(idx_range_1);
    host_t<DFieldMem<IdxRange<GridR<2>>>> result_patch2(idx_range_2);

    double exp_val = exp(5.0 * dt * Nt);
    ddc::for_each(idx_range_1, [&](IdxR1 idx) {
        double const C = (ddc::coordinate(idx) - 0.6);
        result_patch1(idx) = C * exp_val + 0.6;
    });
    ddc::for_each(idx_range_2, [&](IdxR2 idx) {
        double const C = (ddc::coordinate(idx) - 0.6);
        result_patch2(idx) = C * exp_val + 0.6;
    });

    for (int j(0); j < Ntests; ++j) {
        DField<IdxRange<GridR<1>>> vals_patch1 = vals.get<Patch1>();
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_1,
                KOKKOS_LAMBDA(IdxR1 idx) { vals_patch1(idx) = ddc::coordinate(idx); });
        DField<IdxRange<GridR<2>>> vals_patch2 = vals.get<Patch2>();
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_2,
                KOKKOS_LAMBDA(IdxR2 idx) { vals_patch2(idx) = ddc::coordinate(idx); });

        for (int i(0); i < Nt; ++i) {
            timestepper.update(
                    Kokkos::DefaultExecutionSpace(),
                    vals,
                    dt,
                    [&](DMultipatchFieldR dy, DMultipatchConstFieldR y) {
                        DFieldR1 dy_1 = dy.get<Patch1>();
                        DFieldR2 dy_2 = dy.get<Patch2>();
                        DConstFieldR1 y_1 = y.get<Patch1>();
                        DConstFieldR2 y_2 = y.get<Patch2>();
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range_1,
                                KOKKOS_LAMBDA(IdxR1 idx) { dy_1(idx) = 5.0 * y_1(idx) - 3.0; });
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range_2,
                                KOKKOS_LAMBDA(IdxR2 idx) { dy_2(idx) = 5.0 * y_2(idx) - 3.0; });
                    },
                    [&](DMultipatchFieldR y, DMultipatchConstFieldR dy, double dt) {
                        DFieldR1 y_1 = y.get<Patch1>();
                        DFieldR2 y_2 = y.get<Patch2>();
                        DConstFieldR1 dy_1 = dy.get<Patch1>();
                        DConstFieldR2 dy_2 = dy.get<Patch2>();
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range_1,
                                KOKKOS_LAMBDA(IdxR1 idx) { y_1(idx) += dt * dy_1(idx); });
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range_2,
                                KOKKOS_LAMBDA(IdxR2 idx) { y_2(idx) += dt * dy_2(idx); });
                    });
        }

        auto vals1_host = ddc::create_mirror_view_and_copy(vals.get<Patch1>());
        auto vals2_host = ddc::create_mirror_view_and_copy(vals.get<Patch2>());

        double linf_err = 0.0;
        ddc::for_each(idx_range_1, [&](Idx<GridR<1>> idx) {
            double const err = abs(result_patch1(idx) - vals1_host(idx));
            linf_err = err > linf_err ? err : linf_err;
        });
        ddc::for_each(idx_range_2, [&](Idx<GridR<2>> idx) {
            double const err = abs(result_patch2(idx) - vals2_host(idx));
            linf_err = err > linf_err ? err : linf_err;
        });
        error[j] = linf_err;

        dt *= 0.5;
        Nt *= 2;
    }
    for (int j(0); j < Ntests - 1; ++j) {
        order[j] = log(error[j] / error[j + 1]) / log(2.0);
        EXPECT_NEAR(order[j], expected_order, 1e-1);
    }
}

template <class TimeStepper>
void multipatch_timestepper_2D_test(double expected_order)
{
    int constexpr Ntests = 2;
    std::array<double, Ntests> error;
    std::array<double, Ntests - 1> order;

    Coord<R> r_min(0.0);
    Coord<R> r_mid(1.0);
    Coord<R> r_max(2.0);

    IdxStep<GridR<1>> ncells_r_patch_1(10);
    IdxStep<GridR<2>> ncells_r_patch_2(5);

    Idx<GridR<1>> idx_r_start1(0);
    Idx<GridR<2>> idx_r_start2(0);

    Coord<Theta> theta_min(0.0);
    Coord<Theta> theta_max(2 * M_PI);

    IdxStep<GridTheta<1>> ncells_theta_patch_1(10);
    IdxStep<GridTheta<2>> ncells_theta_patch_2(10);

    Idx<GridTheta<1>> idx_theta_start1(0);
    Idx<GridTheta<2>> idx_theta_start2(0);

    double const xc = 0.25;
    double const yc = 0.0;
    double const omega = 2 * M_PI;

    double dt(0.01);
    int Nt(1);

    ddc::init_discrete_space<GridR<1>>(GridR<1>::init(r_min, r_mid, ncells_r_patch_1));
    ddc::init_discrete_space<GridR<2>>(GridR<2>::init(r_mid, r_max, ncells_r_patch_2));
    ddc::init_discrete_space<GridTheta<1>>(
            GridTheta<1>::init(theta_min, theta_max, ncells_theta_patch_1));
    ddc::init_discrete_space<GridTheta<2>>(
            GridTheta<2>::init(theta_min, theta_max, ncells_theta_patch_2));

    IdxRange<GridR<1>> idx_range_r_1(idx_r_start1, ncells_r_patch_1);
    IdxRange<GridR<2>> idx_range_r_2(idx_r_start2, ncells_r_patch_2);
    IdxRange<GridTheta<1>> idx_range_theta_1(idx_theta_start1, ncells_theta_patch_1);
    IdxRange<GridTheta<2>> idx_range_theta_2(idx_theta_start2, ncells_theta_patch_2);
    IdxRange<GridR<1>, GridTheta<1>> idx_range_1(idx_range_r_1, idx_range_theta_1);
    IdxRange<GridR<2>, GridTheta<2>> idx_range_2(idx_range_r_2, idx_range_theta_2);
    MultipatchIdxRangeRTheta idx_ranges(idx_range_1, idx_range_2);

    CMultipatchFieldMemRTheta vals_alloc(idx_ranges);
    CMultipatchFieldRTheta vals(vals_alloc);

    TimeStepper timestepper(idx_ranges);

    host_t<FieldMem<CoordRTheta, IdxRange<GridR<1>, GridTheta<1>>>> result_patch1(idx_range_1);
    host_t<FieldMem<CoordRTheta, IdxRange<GridR<2>, GridTheta<2>>>> result_patch2(idx_range_2);

    double cos_val = std::cos(omega * dt * Nt);
    double sin_val = std::sin(omega * dt * Nt);
    ddc::for_each(idx_range_1, [&](IdxRTheta1 ixy) {
        double const dist_x = (ddc::coordinate(ddc::select<GridR<1>>(ixy)) - xc);
        double const dist_y = (ddc::coordinate(ddc::select<GridTheta<1>>(ixy)) - yc);

        ddc::get<R>(result_patch1(ixy)) = xc + dist_x * cos_val - dist_y * sin_val;
        ddc::get<Theta>(result_patch1(ixy)) = yc + dist_x * sin_val + dist_y * cos_val;
    });
    ddc::for_each(idx_range_2, [&](IdxRTheta2 ixy) {
        double const dist_x = (ddc::coordinate(ddc::select<GridR<2>>(ixy)) - xc);
        double const dist_y = (ddc::coordinate(ddc::select<GridTheta<2>>(ixy)) - yc);

        ddc::get<R>(result_patch2(ixy)) = xc + dist_x * cos_val - dist_y * sin_val;
        ddc::get<Theta>(result_patch2(ixy)) = yc + dist_x * sin_val + dist_y * cos_val;
    });

    for (int j(0); j < Ntests; ++j) {
        CFieldRTheta1 vals_patch1 = vals.get<Patch1>();
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_1,
                KOKKOS_LAMBDA(IdxRTheta1 idx) { vals_patch1(idx) = ddc::coordinate(idx); });
        CFieldRTheta2 vals_patch2 = vals.get<Patch2>();
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_2,
                KOKKOS_LAMBDA(IdxRTheta2 idx) { vals_patch2(idx) = ddc::coordinate(idx); });

        for (int i(0); i < Nt; ++i) {
            timestepper.update(
                    Kokkos::DefaultExecutionSpace(),
                    vals,
                    dt,
                    [yc,
                     xc,
                     &idx_range_2,
                     &idx_range_1,
                     omega](CMultipatchVectorFieldRTheta dy, CMultipatchConstFieldRTheta y) {
                        VectorFieldRTheta1 dy_1 = dy.get<Patch1>();
                        VectorFieldRTheta2 dy_2 = dy.get<Patch2>();
                        CConstFieldRTheta1 y_1 = y.get<Patch1>();
                        CConstFieldRTheta2 y_2 = y.get<Patch2>();
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range_1,
                                KOKKOS_LAMBDA(IdxRTheta1 ixy) {
                                    ddcHelper::get<R>(dy_1)(ixy)
                                            = omega * (yc - ddc::get<Theta>(y_1(ixy)));
                                    ddcHelper::get<Theta>(dy_1)(ixy)
                                            = omega * (ddc::get<R>(y_1(ixy)) - xc);
                                });
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range_2,
                                KOKKOS_LAMBDA(IdxRTheta2 ixy) {
                                    ddcHelper::get<R>(dy_2)(ixy)
                                            = omega * (yc - ddc::get<Theta>(y_2(ixy)));
                                    ddcHelper::get<Theta>(dy_2)(ixy)
                                            = omega * (ddc::get<R>(y_2(ixy)) - xc);
                                });
                    },
                    [&idx_range_1, &idx_range_2](
                            CMultipatchFieldRTheta y,
                            CMultipatchVectorConstFieldRTheta dy,
                            double dt) {
                        CFieldRTheta1 y_1 = y.get<Patch1>();
                        CFieldRTheta2 y_2 = y.get<Patch2>();
                        VectorConstFieldRTheta1 dy_1 = dy.get<Patch1>();
                        VectorConstFieldRTheta2 dy_2 = dy.get<Patch2>();
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range_1,
                                KOKKOS_LAMBDA(IdxRTheta1 ixy) { y_1(ixy) += dt * dy_1(ixy); });
                        ddc::parallel_for_each(
                                Kokkos::DefaultExecutionSpace(),
                                idx_range_2,
                                KOKKOS_LAMBDA(IdxRTheta2 ixy) { y_2(ixy) += dt * dy_2(ixy); });
                    });
        }

        auto vals1_host = ddc::create_mirror_view_and_copy(vals.get<Patch1>());
        auto vals2_host = ddc::create_mirror_view_and_copy(vals.get<Patch2>());

        double linf_err = 0.0;
        ddc::for_each(idx_range_1, [&](IdxRTheta1 idx) {
            double const err_x = ddc::get<R>(result_patch1(idx) - vals1_host(idx));
            double const err_y = ddc::get<Theta>(result_patch1(idx) - vals1_host(idx));
            double const err = std::sqrt(err_x * err_x + err_y * err_y);
            linf_err = err > linf_err ? err : linf_err;
        });
        ddc::for_each(idx_range_2, [&](IdxRTheta2 idx) {
            double const err_x = ddc::get<R>(result_patch2(idx) - vals2_host(idx));
            double const err_y = ddc::get<Theta>(result_patch2(idx) - vals2_host(idx));
            double const err = std::sqrt(err_x * err_x + err_y * err_y);
            linf_err = err > linf_err ? err : linf_err;
        });
        error[j] = linf_err;

        dt *= 0.5;
        Nt *= 2;
    }
    for (int j(0); j < Ntests - 1; ++j) {
        order[j] = log(error[j] / error[j + 1]) / log(2.0);
        EXPECT_NEAR(order[j], expected_order, 1e-1);
    }
}

TEST(EulerFixture, Multipatch)
{
    multipatch_timestepper_test<Euler<DMultipatchFieldMemR>>(1.0);
}

TEST(RK2Fixture, Multipatch)
{
    multipatch_timestepper_test<RK2<DMultipatchFieldMemR>>(2.0);
}

TEST(RK3Fixture, Multipatch)
{
    multipatch_timestepper_test<RK3<DMultipatchFieldMemR>>(3.0);
}

TEST(RK4Fixture, Multipatch)
{
    multipatch_timestepper_test<RK4<DMultipatchFieldMemR>>(4.0);
}

TEST(CrankNicolsonFixture, Multipatch)
{
    multipatch_timestepper_test<CrankNicolson<DMultipatchFieldMemR>>(2.0);
}

TEST(EulerFixture, Multipatch2D)
{
    multipatch_timestepper_2D_test<
            Euler<CMultipatchFieldMemRTheta, CMultipatchVectorFieldMemRTheta>>(1.0);
}

TEST(RK2Fixture, Multipatch2D)
{
    multipatch_timestepper_2D_test<RK2<CMultipatchFieldMemRTheta, CMultipatchVectorFieldMemRTheta>>(
            2.0);
}

TEST(RK3Fixture, Multipatch2D)
{
    multipatch_timestepper_2D_test<RK3<CMultipatchFieldMemRTheta, CMultipatchVectorFieldMemRTheta>>(
            3.0);
}

TEST(RK4Fixture, Multipatch2D)
{
    multipatch_timestepper_2D_test<RK4<CMultipatchFieldMemRTheta, CMultipatchVectorFieldMemRTheta>>(
            4.0);
}

TEST(CrankNicolsonFixture, Multipatch2D)
{
    multipatch_timestepper_2D_test<
            CrankNicolson<CMultipatchFieldMemRTheta, CMultipatchVectorFieldMemRTheta>>(2.0);
}

} // namespace
