// SPDX-License-Identifier: MIT
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "polar_bsplines.hpp"
#include "polar_spline.hpp"

struct R
{
    static constexpr bool PERIODIC = false;
};
struct Theta
{
    static constexpr bool PERIODIC = true;
};

static constexpr std::size_t spline_r_degree = 3;
static constexpr std::size_t spline_theta_degree = 3;
static constexpr int continuity = 1;



struct BSplinesR : ddc::NonUniformBSplines<R, spline_r_degree>
{
};
struct BSplinesTheta : ddc::UniformBSplines<Theta, spline_theta_degree>
{
};


using GrevillePointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using GrevillePointsTheta = ddc::GrevilleInterpolationPoints<
        BSplinesTheta,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;

struct GridR : GrevillePointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : GrevillePointsTheta::interpolation_discrete_dimension_type
{
};
struct BSplines : PolarBSplines<BSplinesR, BSplinesTheta, continuity>
{
};

using PolarCoord = ddc::Coordinate<R, Theta>;
using CoordR = ddc::Coordinate<R>;
using CoordTheta = ddc::Coordinate<Theta>;
using Spline_h = PolarSpline<BSplines, Kokkos::HostSpace>;
using Spline_d = PolarSpline<BSplines>;
using SplineSpan_h = PolarSplineSpan<BSplines, Kokkos::HostSpace>;
using SplineSpan_d = PolarSplineSpan<BSplines>;

template <class StartExecSpace, class EndExecSpace>
class MirrorTestClass
{
    using Spline = PolarSpline<BSplines, typename StartExecSpace::memory_space>;
    using BuilderRTheta = ddc::SplineBuilder2D<
            StartExecSpace,
            typename StartExecSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::LAPACK,
            GridR,
            GridTheta>;

    CoordR const r0 {0.};
    CoordR const rN {1.};
    CoordTheta const theta0 {0.};
    CoordTheta const thetaN {2. * M_PI};
    std::size_t const ncells;

public:
    MirrorTestClass(std::size_t n_cells) : ncells(n_cells)
    {
        ddc::DiscreteVector<GridR> npoints_r(ncells + 1);
        std::vector<CoordR> breaks_r(npoints_r);
        const double dr = (rN - r0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_r[i] = CoordR(r0 + i * dr);
        }
        ddc::init_discrete_space<BSplinesR>(breaks_r);
        ddc::init_discrete_space<BSplinesTheta>(theta0, thetaN, ncells);
        ddc::init_discrete_space<GridR>(GrevillePointsR::get_sampling<GridR>());
        ddc::init_discrete_space<GridTheta>(GrevillePointsTheta::get_sampling<GridTheta>());
    }

    Spline set_up()
    {
        ddc::DiscreteDomain<GridR, GridTheta> interpolation_idx_range(
                GrevillePointsR::get_domain<GridR>(),
                GrevillePointsTheta::get_domain<GridTheta>());
        BuilderRTheta builder_rtheta(interpolation_idx_range);
        Spline coef(builder_rtheta.spline_domain());
        ddc::parallel_fill(coef.singular_spline_coef, 42);
        ddc::parallel_fill(coef.spline_coef, 13);
        return coef;
    }

    bool check_mirror_same_attributes(
            Spline& coef,
            PolarSplineSpan<BSplines, typename EndExecSpace::memory_space> copy_test)
    {
        return (coef.spline_coef.domain() == copy_test.spline_coef.domain()
                && coef.singular_spline_coef.domain() == copy_test.singular_spline_coef.domain());
    }

    bool check_mirror_equals_copy(
            Spline& coef,
            PolarSplineSpan<BSplines, typename EndExecSpace::memory_space> copy_test)
    {
        check_mirror_same_attributes(coef, copy_test);
        auto copy_singular
                = ddc::create_mirror_and_copy(copy_test.singular_spline_coef.span_cview());
        auto coef_singular = ddc::create_mirror_and_copy(coef.singular_spline_coef.span_cview());
        auto copy_regular = ddc::create_mirror_and_copy(copy_test.spline_coef.span_cview());
        auto coef_regular = ddc::create_mirror_and_copy(coef.spline_coef.span_cview());

        bool const is_equal_singular_values = ddc::transform_reduce(
                coef.singular_spline_coef.domain(),
                true,
                ddc::reducer::land<bool>(),
                [&](auto idx) { return (coef_singular(idx) == copy_singular(idx)); });

        bool const is_equal_values = ddc::transform_reduce(
                coef.spline_coef.domain(),
                true,
                ddc::reducer::land<bool>(),
                [&](auto idx) { return (coef_regular(idx) == copy_regular(idx)); });

        return (is_equal_singular_values && is_equal_values);
    }
};

TEST(Create_mirror, Host_to_Host)
{
    MirrorTestClass<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace>
            host_to_host(20);
    Spline_h coef = host_to_host.set_up();
    auto copy_coef = create_mirror(coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), Spline_h>);
    ASSERT_TRUE(host_to_host.check_mirror_same_attributes(coef, copy_coef));
}

TEST(Create_mirror, Host_to_Device)
{
    MirrorTestClass<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>
            host_to_device(20);
    Spline_h coef = host_to_device.set_up();
    auto copy_coef = create_mirror(Kokkos::DefaultExecutionSpace(), coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), Spline_d>);
    ASSERT_TRUE(host_to_device.check_mirror_same_attributes(coef, copy_coef.span_view()));
}

TEST(Create_mirror, Device_to_Host)
{
    MirrorTestClass<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>
            device_to_host(20);
    Spline_d coef = device_to_host.set_up();
    auto copy_coef = create_mirror(coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), Spline_h>);
    ASSERT_TRUE(device_to_host.check_mirror_same_attributes(coef, copy_coef));
}

TEST(Create_mirror, Device_to_Device)
{
    MirrorTestClass<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace> device_to_device(
            20);
    Spline_d coef = device_to_device.set_up();
    auto copy_coef = create_mirror(Kokkos::DefaultExecutionSpace(), coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), Spline_d>);
    ASSERT_TRUE(device_to_device.check_mirror_same_attributes(coef, copy_coef));
}

TEST(Create_mirror_and_copy, Host_to_Host)
{
    MirrorTestClass<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace>
            host_to_host(20);
    Spline_h coef = host_to_host.set_up();
    auto copy_coef = create_mirror_and_copy(coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), Spline_h>);
    ASSERT_TRUE(host_to_host.check_mirror_equals_copy(coef, copy_coef));
}

TEST(Create_mirror_and_copy, Host_to_Device)
{
    MirrorTestClass<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>
            host_to_device(20);
    Spline_h coef = host_to_device.set_up();
    auto copy_coef = create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), Spline_d>);
    ASSERT_TRUE(host_to_device.check_mirror_equals_copy(coef, copy_coef));
}

TEST(Create_mirror_and_copy, Device_to_Host)
{
    MirrorTestClass<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>
            device_to_host(20);
    Spline_d coef = device_to_host.set_up();
    auto copy_coef = create_mirror_and_copy(coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), Spline_h>);
    ASSERT_TRUE(device_to_host.check_mirror_equals_copy(coef, copy_coef));
}

TEST(Create_mirror_and_copy, Device_to_Device)
{
    MirrorTestClass<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace> device_to_device(
            20);
    Spline_d coef = device_to_device.set_up();
    auto copy_coef = create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), Spline_d>);
    ASSERT_TRUE(device_to_device.check_mirror_equals_copy(coef, copy_coef));
}

TEST(Create_mirror_view, Host_to_Host)
{
    MirrorTestClass<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace>
            host_to_host(20);
    Spline_h coef = host_to_host.set_up();
    auto copy_coef = create_mirror_view(coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), SplineSpan_h>);
    ASSERT_TRUE(host_to_host.check_mirror_same_attributes(coef, copy_coef.span_view()));
}

TEST(Create_mirror_view, Host_to_Device)
{
    MirrorTestClass<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>
            host_to_device(20);
    Spline_h coef = host_to_device.set_up();
    auto copy_coef = create_mirror_view(Kokkos::DefaultExecutionSpace(), coef.span_view());
    static_assert(std::is_same_v<
                  decltype(copy_coef),
                  std::conditional_t<
                          std::is_same_v<
                                  Kokkos::DefaultExecutionSpace::memory_space,
                                  Kokkos::HostSpace>,
                          SplineSpan_d,
                          Spline_d>>);

    ASSERT_TRUE(host_to_device.check_mirror_same_attributes(coef, copy_coef.span_view()));
}

TEST(Create_mirror_view, Device_to_Host)
{
    MirrorTestClass<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>
            device_to_host(20);
    Spline_d coef = device_to_host.set_up();
    auto copy_coef = create_mirror_view(coef.span_view());
    static_assert(std::is_same_v<
                  decltype(copy_coef),
                  std::conditional_t<
                          std::is_same_v<
                                  Kokkos::DefaultExecutionSpace::memory_space,
                                  Kokkos::HostSpace>,
                          SplineSpan_h,
                          Spline_h>>);

    ASSERT_TRUE(device_to_host.check_mirror_same_attributes(coef, copy_coef.span_view()));
}

TEST(Create_mirror_view, Device_to_Device)
{
    MirrorTestClass<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace> device_to_device(
            20);
    Spline_d coef = device_to_device.set_up();
    auto copy_coef = create_mirror_view(Kokkos::DefaultExecutionSpace(), coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), SplineSpan_d>);
    ASSERT_TRUE(device_to_device.check_mirror_same_attributes(coef, copy_coef.span_view()));
}

TEST(Create_mirror_view_and_copy, Host_to_Host)
{
    MirrorTestClass<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace>
            host_to_host(20);
    Spline_h coef = host_to_host.set_up();
    auto copy_coef = create_mirror_view_and_copy(coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), SplineSpan_h>);
    ASSERT_TRUE(host_to_host.check_mirror_equals_copy(coef, copy_coef));
}

TEST(Create_mirror_view_and_copy, Host_to_Device)
{
    MirrorTestClass<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>
            host_to_device(20);
    Spline_h coef = host_to_device.set_up();
    auto copy_coef = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), coef.span_view());
    static_assert(std::is_same_v<
                  decltype(copy_coef),
                  std::conditional_t<
                          std::is_same_v<
                                  Kokkos::DefaultExecutionSpace::memory_space,
                                  Kokkos::HostSpace>,
                          SplineSpan_d,
                          Spline_d>>);

    ASSERT_TRUE(host_to_device.check_mirror_equals_copy(coef, copy_coef));
}

TEST(Create_mirror_view_and_copy, Device_to_Host)
{
    MirrorTestClass<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>
            device_to_host(20);
    Spline_d coef = device_to_host.set_up();
    auto copy_coef = create_mirror_view_and_copy(coef.span_view());

    static_assert(std::is_same_v<
                  decltype(copy_coef),
                  std::conditional_t<
                          std::is_same_v<
                                  Kokkos::DefaultExecutionSpace::memory_space,
                                  Kokkos::HostSpace>,
                          SplineSpan_h,
                          Spline_h>>);

    ASSERT_TRUE(device_to_host.check_mirror_equals_copy(coef, copy_coef));
}

TEST(Create_mirror_view_and_copy, Device_to_Device)
{
    MirrorTestClass<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace> device_to_device(
            20);
    Spline_d coef = device_to_device.set_up();
    auto copy_coef = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), coef.span_view());
    static_assert(std::is_same_v<decltype(copy_coef), SplineSpan_d>);
    ASSERT_TRUE(device_to_device.check_mirror_equals_copy(coef, copy_coef));
}
