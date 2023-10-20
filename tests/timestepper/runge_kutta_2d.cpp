// SPDX-License-Identifier: MIT
#include <tuple>
#include <type_traits>
#include <utility>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <directional_tag.hpp>
#include <rk2.hpp>
#include <rk3.hpp>
#include <rk4.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

using namespace ddc;


template <class T>
class RungeKutta2DFixture;

template <std::size_t ORDER>
class RungeKutta2DFixture<std::tuple<std::integral_constant<std::size_t, ORDER>>>
    : public testing::Test
{
public:
    static int constexpr order = ORDER;

    struct RDimX
    {
        static bool constexpr PERIODIC = false;
    };

    struct RDimY
    {
        static bool constexpr PERIODIC = false;
    };
    using CoordX = Coordinate<RDimX>;
    using CoordY = Coordinate<RDimY>;
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
    using RungeKutta = std::conditional_t<
            ORDER == 2,
            RK2<AdvectionField>,
            std::conditional_t<ORDER == 3, RK3<AdvectionField>, RK4<AdvectionField>>>;
};

using runge_kutta_2d_types = testing::Types<
        std::tuple<std::integral_constant<std::size_t, 2>>,
        std::tuple<std::integral_constant<std::size_t, 3>>,
        std::tuple<std::integral_constant<std::size_t, 4>>>;

TYPED_TEST_SUITE(RungeKutta2DFixture, runge_kutta_2d_types);

TYPED_TEST(RungeKutta2DFixture, RungeKutta2DOrder)
{
    using RDimX = typename TestFixture::RDimX;
    using CoordX = typename TestFixture::CoordX;
    using IDimX = typename TestFixture::IDimX;
    using IndexX = typename TestFixture::IndexX;
    using IVectX = typename TestFixture::IVectX;
    using IDomainX = typename TestFixture::IDomainX;

    using RDimY = typename TestFixture::RDimY;
    using CoordY = typename TestFixture::CoordY;
    using IDimY = typename TestFixture::IDimY;
    using IndexY = typename TestFixture::IndexY;
    using IVectY = typename TestFixture::IVectY;
    using IDomainY = typename TestFixture::IDomainY;

    using IndexXY = typename TestFixture::IndexXY;
    using IDomainXY = typename TestFixture::IDomainXY;
    using RungeKutta = typename TestFixture::RungeKutta;
    using AdvectionField = typename TestFixture::AdvectionField;
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

    RungeKutta runge_kutta(dom);

    AdvectionField vals(dom);
    AdvectionField result(dom);

    double cos_val = std::cos(omega * dt * Nt);
    double sin_val = std::sin(omega * dt * Nt);
    for_each(dom, [&](IndexXY ixy) {
        double const dist_x = (coordinate(select<IDimX>(ixy)) - xc);
        double const dist_y = (coordinate(select<IDimY>(ixy)) - yc);

        ddcHelper::get<RDimX>(result)(ixy) = xc + dist_x * cos_val - dist_y * sin_val;
        ddcHelper::get<RDimY>(result)(ixy) = yc + dist_x * sin_val + dist_y * cos_val;
    });

    for (int j(0); j < Ntests; ++j) {
        for_each(dom, [&](IndexXY ixy) {
            ddcHelper::get<RDimX>(vals)(ixy) = coordinate(select<IDimX>(ixy));
            ddcHelper::get<RDimY>(vals)(ixy) = coordinate(select<IDimY>(ixy));
        });

        for (int i(0); i < Nt; ++i) {
            runge_kutta.update(
                    vals,
                    dt,
                    [yc, xc, &dom, omega](AdvectionFieldSpan dy, AdvectionFieldView y) {
                        for_each(dom, [&](IndexXY ixy) {
                            ddcHelper::get<RDimX>(dy)(ixy)
                                    = omega * (yc - ddcHelper::get<RDimY>(y)(ixy));
                            ddcHelper::get<RDimY>(dy)(ixy)
                                    = omega * (ddcHelper::get<RDimX>(y)(ixy) - xc);
                        });
                    },
                    [&dom](AdvectionFieldSpan y, AdvectionFieldView dy, double dt) {
                        for_each(dom, [&](IndexXY ixy) {
                            ddcHelper::get<RDimX>(y)(ixy) += ddcHelper::get<RDimX>(dy)(ixy) * dt;
                            ddcHelper::get<RDimY>(y)(ixy) += ddcHelper::get<RDimY>(dy)(ixy) * dt;
                        });
                    });
        }

        double linf_err = 0.0;
        for_each(dom, [&](IndexXY ixy) {
            double const err_x
                    = ddcHelper::get<RDimX>(result)(ixy) - ddcHelper::get<RDimX>(vals)(ixy);
            double const err_y
                    = ddcHelper::get<RDimY>(result)(ixy) - ddcHelper::get<RDimY>(vals)(ixy);
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
