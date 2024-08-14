#include <fstream>
#include <random>
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Lagrange_interpolator.hpp"
#include "mesh_builder.hpp"
// Here Z stands to specify that we use ad hoc geometry
struct Z
{
    static bool constexpr PERIODIC = false;
};

using CoordZ = Coord<Z>;



template <std::size_t N>
class LagrangeTest
{
public:
    struct GridZ : NonUniformGridBase<Z>
    {
    };

private:
    using IdxStepZ = IdxStep<GridZ>;
    using IdxZ = Idx<GridZ>;
    const IdxStepZ x_ncells;
    const CoordZ x_min;
    const CoordZ x_max;
    const bool perturb;
    std::function<double(double, int)> const f_poly;
    std::vector<double> error;
    std::pair<double, double> err_out;


public:
    LagrangeTest() : x_ncells(10), x_min(-1), x_max(1), perturb(false) {};
    LagrangeTest(int nx, int deg, bool perturbation, std::function<double(double, int)> f)
        : x_ncells(nx)
        , x_min(-1)
        , x_max(1)
        , perturb(perturbation)
        , f_poly(f)
    {
        std::vector<CoordZ> point_sampling = build_uniform_break_points(x_min, x_max, x_ncells);
        ddc::init_discrete_space<GridZ>(point_sampling);
        Idx<GridZ> lbound(0);
        IdxRange<GridZ> interpolation_idx_range_x(lbound, x_ncells + 1);

        host_t<DFieldMem<IdxRange<GridZ>>> exact_host_alloc(interpolation_idx_range_x);
        host_t<DField<IdxRange<GridZ>>> exact_host = get_field(exact_host_alloc);
        ddc::for_each(get_idx_range(exact_host), [&](IdxZ const ix) {
            exact_host(ix) = f_poly(ddc::coordinate(ix).value(), deg);
        });
        int cpt = 1;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> unif_distr(-1, 1);
        host_t<FieldMem<Coord<Z>, IdxRange<GridZ>>> Sample_host(interpolation_idx_range_x);

        ddc::for_each(get_idx_range(Sample_host), [&](IdxZ const ix) {
            if (cpt % int(0.1 * x_ncells) == 0) {
                Sample_host(ix) = Coord<Z>(
                        ddc::coordinate(ix).value()
                        + (cpt > 0.1 * x_ncells) * (cpt < 0.9 * x_ncells) * perturb
                                  * unif_distr(generator) / x_ncells);
            } else {
                Sample_host(ix) = Coord<Z>(ddc::coordinate(ix).value());
            }
            cpt++;
        });
        IdxStep<GridZ> static constexpr gwx {0};

        LagrangeInterpolator<GridZ, BCond::DIRICHLET, BCond::DIRICHLET, GridZ>
                Test_Interpolator(deg, gwx);
        auto exact = ddc::
                create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), get_field(exact_host));
        auto Sample = ddc::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(Sample_host));
        Test_Interpolator(get_field(exact), get_const_field(Sample));
        ddc::parallel_deepcopy(exact_host, exact);
        ddc::parallel_deepcopy(Sample_host, Sample);
        ddc::for_each(interpolation_idx_range_x, [&](IdxZ const ix) {
            error.push_back(std::abs(exact_host(ix) - f_poly(Sample_host(ix), deg)));
        });
        double s = std::transform_reduce(
                std::begin(error),
                std::end(error),
                0.0,
                std::plus<double>(),
                [](double x) { return x * x; });
        err_out.first = *std::max_element(error.begin(), error.end()); // L_infinity error
        err_out.second = sqrt(s); // L² error
    }
    std::pair<double, double> get_errors()
    {
        return err_out;
    }
};

class LagrangeTestFixture : public ::testing::TestWithParam<int>
{
protected:
    LagrangeTest<100> lagrangetest;
};


TEST_P(LagrangeTestFixture, ExactNodes)
{
    auto deg = GetParam();
    LagrangeTest<20> Test_instance_twenty(20, deg, false, [](double x, int d) {
        return std::pow(x, d);
    });
    LagrangeTest<50> Test_instance_fifty(50, deg, false, [](double x, int d) {
        return std::pow(x, d);
    });
    LagrangeTest<100> Test_instance_hundred(100, deg, false, [](double x, int d) {
        return std::pow(x, d);
    });

    EXPECT_EQ(Test_instance_twenty.get_errors().first, 0.);
    EXPECT_EQ(Test_instance_fifty.get_errors().first, 0.);
    EXPECT_EQ(Test_instance_hundred.get_errors().first, 0.);
}

TEST_P(LagrangeTestFixture, ExactPolynomial)
{
    auto deg = GetParam();

    LagrangeTest<100> Test_instance_hundred(100, deg, true, [](double x, int d) {
        return std::pow(x, d);
    });
    double tol = 1e-14;

    EXPECT_LE(Test_instance_hundred.get_errors().first, tol);
    EXPECT_LE(Test_instance_hundred.get_errors().second, tol);
}


TEST_P(LagrangeTestFixture, InterpolationError)
{
    auto deg = GetParam();
    double Linf[3];

    LagrangeTest<25> Test_instance_twenty(25, deg, true, [](double x, int d) {
        return 1. / (1. + x * x);
    });
    LagrangeTest<50> Test_instance_fifty(50, deg, true, [](double x, int d) {
        return 1. / (1. + x * x);
    });
    LagrangeTest<100> Test_instance_hundred(100, deg, true, [](double x, int d) {
        return 1. / (1. + x * x);
    });
    Linf[0] = Test_instance_twenty.get_errors().first;
    Linf[1] = Test_instance_fifty.get_errors().first;
    Linf[2] = Test_instance_hundred.get_errors().first;
    double rate = (log(Linf[1]) - log(Linf[0])) / (log(50) - log(25));
    std ::cout << "Degree of Lagrange polynomial " << deg << std::endl
               << " Associated convergence rate " << rate << std::endl;
    EXPECT_LE(rate, -deg);
}

INSTANTIATE_TEST_SUITE_P(DegMesh, LagrangeTestFixture, ::testing::Values(3, 4, 5, 6, 7, 8, 9));
