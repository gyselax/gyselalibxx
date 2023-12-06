#include <fstream>
#include <random>
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Lagrange_interpolator.hpp"
// Here Z stands to specify that we use ad hoc geometry
struct RDimZ
{
    static bool constexpr PERIODIC = false;
};

template <std::size_t N>
struct CDimZ;
using CoordZ = ddc::Coordinate<RDimZ>;



template <std::size_t N>
class LagrangeTest
{
    using IDimZ = ddc::NonUniformPointSampling<CDimZ<N>>;
    using IVectZ = ddc::DiscreteVector<IDimZ>;
    using CDim = typename IDimZ::continuous_dimension_type;
    using DElemZ = ddc::DiscreteElement<IDimZ>;
    const IVectZ x_size;
    const CoordZ x_min;
    const CoordZ x_max;
    const bool perturb;
    std::function<double(double, int)> const fpol;
    std::vector<double> Erreur;
    std::pair<double, double> err_out;


public:
    LagrangeTest() : x_size(10), x_min(-1), x_max(1), perturb(false) {};
    LagrangeTest(int nx, int deg, bool perturbation, std::function<double(double, int)> f)
        : x_size(nx)
        , x_min(-1)
        , x_max(1)
        , perturb(perturbation)
        , fpol(f)
    {
        std::vector<double> point_sampling;

        double dx = (x_max - x_min) / x_size;
        // Create & intialize Uniform mesh
        for (int k = 0; k < x_size; k++) {
            point_sampling.push_back(x_min + k * dx);
        }
        ddc::init_discrete_space<IDimZ>(point_sampling);
        ddc::DiscreteElement<IDimZ> lbound(0);
        ddc::DiscreteVector<IDimZ> npoints(x_size);
        ddc::DiscreteDomain<IDimZ> interpolation_domain_x(lbound, npoints);

        ddc::Chunk Essai_alloc(interpolation_domain_x, ddc::HostAllocator<double>());
        ddc::ChunkSpan Essai = Essai_alloc.span_view();
        ddc::for_each(Essai.domain(), [&](DElemZ const ix) {
            Essai(ix) = fpol(ddc::coordinate(ix).value(), deg);
        });
        int cpt = 1;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> unif_distr(-1, 1);
        ddc::Chunk<ddc::Coordinate<CDim>, ddc::DiscreteDomain<IDimZ>> Sample(
                interpolation_domain_x);

        ddc::for_each(Sample.domain(), [&](DElemZ const ix) {
            if (cpt % int(0.1 * x_size) == 0) {
                Sample(ix) = ddc::Coordinate<CDim>(
                        ddc::coordinate(ix).value()
                        + (cpt > 0.1 * x_size) * (cpt < 0.9 * x_size) * perturb
                                  * unif_distr(generator) / x_size);
            } else {
                Sample(ix) = ddc::Coordinate<CDim>(ddc::coordinate(ix).value());
            }
            cpt++;
        });
        ddc::DiscreteVector<IDimZ> static constexpr gwx {0};

        LagrangeInterpolator<IDimZ, BCond::DIRICHLET, BCond::DIRICHLET>
                Test_Interpolator(deg, Essai.domain(), gwx);
        Test_Interpolator(Essai, Sample);
        ddc::for_each(interpolation_domain_x, [&](DElemZ const ix) {
            Erreur.push_back(std::abs(Essai(ix) - fpol(Sample(ix), deg)));
        });
        double s = std::transform_reduce(
                std::begin(Erreur),
                std::end(Erreur),
                0.0,
                std::plus<double>(),
                [](double x) { return x * x; });
        err_out.first = *std::max_element(Erreur.begin(), Erreur.end()); // L_infinity error
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
