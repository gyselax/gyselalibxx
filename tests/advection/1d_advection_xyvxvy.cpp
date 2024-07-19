
#include <cmath>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "bsl_advection_1d.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "rk2.hpp"
#include "spline_interpolator.hpp"
#include "vector_field_common.hpp"


namespace {
// Continuous dimensions
/// @brief A class which describes the real space in the first spatial direction X.
struct RDimX
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

/// @brief A class which describes the real space in the second spatial direction Y.
struct RDimY
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

/// @brief A class which describes the real space in the second velocity direction X.
struct RDimVx
{
};

/// @brief A class which describes the real space in the second velocity direction Y.
struct RDimVy
{
};


using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;
using CoordXY = ddc::Coordinate<RDimX, RDimY>;

using CoordVx = ddc::Coordinate<RDimVx>;
using CoordVy = ddc::Coordinate<RDimVy>;

// Splines
struct BSplinesX : ddc::UniformBSplines<RDimX, 3>
{
};
struct BSplinesY : ddc::UniformBSplines<RDimY, 3>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineYBoundary = ddc::BoundCond::PERIODIC;

// Discrete dimensions
struct IDimX : ddc::UniformPointSampling<RDimX>
{
};
struct IDimY : ddc::UniformPointSampling<RDimY>
{
};
struct IDimVx : ddc::UniformPointSampling<RDimVx>
{
};
struct IDimVy : ddc::UniformPointSampling<RDimVy>
{
};


using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;


using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;
using IndexVx = ddc::DiscreteElement<IDimVx>;
using IndexVy = ddc::DiscreteElement<IDimVy>;
using IndexXYVxVy = ddc::DiscreteElement<IDimX, IDimY, IDimVx, IDimVy>;


using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;
using IVectVx = ddc::DiscreteVector<IDimVx>;
using IVectVy = ddc::DiscreteVector<IDimVy>;


using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
using IDomainVx = ddc::DiscreteDomain<IDimVx>;
using IDomainVy = ddc::DiscreteDomain<IDimVy>;
using IDomainXYVxVy = ddc::DiscreteDomain<IDimX, IDimY, IDimVx, IDimVy>;



// Chunks, Spans and Views
template <class ElementType>
using FieldXY = device_t<ddc::Chunk<ElementType, IDomainXY>>;
using DFieldXY = FieldXY<double>;

template <class ElementType>
using FieldXYVxVy = device_t<ddc::Chunk<ElementType, IDomainXYVxVy>>;
using DFieldXYVxVy = FieldXYVxVy<double>;


template <class ElementType>
using SpanXY = device_t<ddc::ChunkSpan<ElementType, IDomainXY>>;
using DSpanXY = SpanXY<double>;

template <class ElementType>
using SpanXYVxVy = device_t<ddc::ChunkSpan<ElementType, IDomainXYVxVy>>;
using DSpanXYVxVy = SpanXYVxVy<double>;


// Operators
using SplineXBuilder_2d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY>;
using SplineXEvaluator_2d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimX,
        IDimY>;

using SplineYBuilder_2d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        IDimY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY>;
using SplineYEvaluator_2d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        IDimY,
        ddc::PeriodicExtrapolationRule<RDimY>,
        ddc::PeriodicExtrapolationRule<RDimY>,
        IDimX,
        IDimY>;


using SplineXBuilder_4d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineXEvaluator_4d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;

using SplineYBuilder_4d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        IDimY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineYEvaluator_4d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        IDimY,
        ddc::PeriodicExtrapolationRule<RDimY>,
        ddc::PeriodicExtrapolationRule<RDimY>,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;



class XYVxVyAdvection1DTest : public ::testing::Test
{
protected:
    static constexpr CoordX x_min = CoordX(-.5);
    static constexpr CoordX x_max = CoordX(.5);
    static constexpr IVectX x_size = IVectX(60);

    static constexpr CoordY y_min = CoordY(-.5);
    static constexpr CoordY y_max = CoordY(.5);
    static constexpr IVectY y_size = IVectY(60);

    static constexpr IndexVx idx0_vx = IndexVx(0);
    static constexpr IVectVx vx_size = IVectVx(2);

    static constexpr IndexVy idx0_vy = IndexVy(0);
    static constexpr IVectVy vy_size = IVectVy(2);


    IDomainX const interpolation_domain_x;
    IDomainY const interpolation_domain_y;
    IDomainXY const xy_grid;

    IDomainVx const domain_vx;
    IDomainVy const domain_vy;
    IDomainXYVxVy const xyvxvy_grid;

public:
    XYVxVyAdvection1DTest()
        : interpolation_domain_x(SplineInterpPointsX::get_domain<IDimX>())
        , interpolation_domain_y(SplineInterpPointsY::get_domain<IDimY>())
        , xy_grid(interpolation_domain_x, interpolation_domain_y)
        , domain_vx(idx0_vx, vx_size)
        , domain_vy(idx0_vy, vy_size)
        , xyvxvy_grid(interpolation_domain_x, interpolation_domain_y, domain_vx, domain_vy) {};

    ~XYVxVyAdvection1DTest() = default;


    static void SetUpTestSuite()
    {
        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<BSplinesY>(y_min, y_max, y_size);

        ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
        ddc::init_discrete_space<IDimY>(SplineInterpPointsY::get_sampling<IDimY>());
    }

    template <class AdvectionOperatorX, class AdvectionOperatorY>
    double AdvectionXY(AdvectionOperatorX const& advection_x, AdvectionOperatorY const& advection_y)
    {
        // TIME PARAMETERS ---------------------------------------------------------------------------
        double const dt = 0.05;
        double const final_t = 0.2;
        int const time_iter = int(final_t / dt);

        // INITIALISATION ----------------------------------------------------------------------------
        // Parameters
        double const xc = 0.1;
        double const yc = 0.2;
        double const a = 0.5;

        DFieldXYVxVy function_alloc(xyvxvy_grid);
        DSpanXYVxVy function = function_alloc.span_view();
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                xyvxvy_grid,
                KOKKOS_LAMBDA(IndexXYVxVy const idx) {
                    CoordXY coord_xy = CoordXY(ddc::coordinate(idx));
                    double const x = CoordX(coord_xy);
                    double const y = CoordY(coord_xy);

                    double const r1 = Kokkos::sqrt((x - xc) * (x - xc) + 8 * (y - yc) * (y - yc));
                    double const r2 = Kokkos::sqrt(8 * (x - xc) * (x - xc) + (y - yc) * (y - yc));
                    double const G1
                            = Kokkos::pow(Kokkos::cos(M_PI * r1 / 2. / a), 4) * (abs(r1) < a);
                    double const G2
                            = Kokkos::pow(Kokkos::cos(M_PI * r2 / 2. / a), 4) * (abs(r2) < a);
                    function(idx) = 0.5 * (G1 + G2);
                });

        /*
            The advection field is a Chunk of double instead of a Chunk 
            of coordinates, because inside the advection operator
            we use an interpolator which needs the DFieldXY format.
        */
        DFieldXY advection_field_x_alloc(xy_grid);
        DSpanXY advection_field_x = advection_field_x_alloc.span_view();

        DFieldXY advection_field_y_alloc(xy_grid);
        DSpanXY advection_field_y = advection_field_y_alloc.span_view();

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                xy_grid,
                KOKKOS_LAMBDA(IndexXY const idx) {
                    CoordXY coord_xy(ddc::coordinate(idx));
                    double const x = CoordX(coord_xy);
                    advection_field_x(idx) = Kokkos::sin(x * 2 * M_PI) / M_PI / 2.;

                    double const y = CoordY(coord_xy);
                    advection_field_y(idx) = Kokkos::sin(y * 2 * M_PI) / M_PI / 2.;
                });


        // EXACT ADVECTED FUNCTION -------------------------------------------------------------------
        host_t<DFieldXYVxVy> exact_function(xyvxvy_grid);
        ddc::for_each(xyvxvy_grid, [&](IndexXYVxVy const idx) {
            CoordXY coord_xy = CoordXY(ddc::coordinate(idx));
            double const x0 = CoordX(coord_xy);
            double const y0 = CoordY(coord_xy);

            double x = 2 * std::atan(std::tan(x0 * M_PI) * std::exp(-final_t)) / M_PI / 2.;
            double y = 2 * std::atan(std::tan(y0 * M_PI) * std::exp(-final_t)) / M_PI / 2.;

            // Replace inside the domain the feet if the dimension if periodic
            if (RDimX::PERIODIC) {
                x = fmod(x - double(x_min), double(x_max - x_min)) + double(x_min);
                x = x > double(x_min) ? x : x + double(x_max - x_min);
            }
            if (RDimY::PERIODIC) {
                y = fmod(y - double(y_min), double(y_max - y_min)) + double(y_min);
                y = y > double(y_min) ? y : y + double(y_max - y_min);
            }

            double const r1 = std::sqrt((x - xc) * (x - xc) + 8 * (y - yc) * (y - yc));
            double const r2 = std::sqrt(8 * (x - xc) * (x - xc) + (y - yc) * (y - yc));
            double const G1 = std::pow(std::cos(M_PI * r1 / 2. / a), 4) * (abs(r1) < a);
            double const G2 = std::pow(std::cos(M_PI * r2 / 2. / a), 4) * (abs(r2) < a);
            exact_function(idx) = 0.5 * (G1 + G2);
        });


        // SIMULATION --------------------------------------------------------------------------------
        for (int i(0); i < time_iter; i++) {
            advection_x(function, advection_field_x, dt / 2);
            advection_y(function, advection_field_y, dt);
            advection_x(function, advection_field_x, dt / 2);
        };


        // CHECK ERRORS ------------------------------------------------------------------------------
        /*
            Simulation launched on GPU but error checking on CPU. 
        */
        auto function_host
                = ddc::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), function);
        double max_relative_error = 0;
        ddc::for_each(xyvxvy_grid, [&](IndexXYVxVy const idx) {
            double const relative_error = abs(function_host(idx) - exact_function(idx));
            max_relative_error
                    = max_relative_error > relative_error ? max_relative_error : relative_error;
        });
        return max_relative_error;
    };
};

} // end namespace



TEST_F(XYVxVyAdvection1DTest, AdvectionXY)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    SplineXBuilder_2d const adv_field_builder_x(xy_grid);
    SplineYBuilder_2d const adv_field_builder_y(xy_grid);
    SplineXBuilder_4d const function_builder_x(xyvxvy_grid);
    SplineYBuilder_4d const function_builder_y(xyvxvy_grid);


    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
    SplineXEvaluator_2d const adv_field_spline_evaluator_x(bv_x_min, bv_x_max);
    SplineXEvaluator_4d const function_spline_evaluator_x(bv_x_min, bv_x_max);

    ddc::PeriodicExtrapolationRule<RDimY> bv_y_min;
    ddc::PeriodicExtrapolationRule<RDimY> bv_y_max;
    SplineYEvaluator_2d const adv_field_spline_evaluator_y(bv_y_min, bv_y_max);
    SplineYEvaluator_4d const function_spline_evaluator_y(bv_y_min, bv_y_max);


    PreallocatableSplineInterpolator const
            function_spline_x_interpolator(function_builder_x, function_spline_evaluator_x);
    PreallocatableSplineInterpolator const
            function_spline_y_interpolator(function_builder_y, function_spline_evaluator_y);


    RK2<FieldXY<CoordX>, DFieldXY> time_stepper_x(xy_grid);
    RK2<FieldXY<CoordY>, DFieldXY> time_stepper_y(xy_grid);

    BslAdvection1D<
            IDimX,
            IDomainXY,
            IDomainXYVxVy,
            SplineXBuilder_2d,
            SplineXEvaluator_2d,
            RK2<FieldXY<CoordX>, DFieldXY>> const
            advection_x(
                    function_spline_x_interpolator,
                    adv_field_builder_x,
                    adv_field_spline_evaluator_x,
                    time_stepper_x);
    BslAdvection1D<
            IDimY,
            IDomainXY,
            IDomainXYVxVy,
            SplineYBuilder_2d,
            SplineYEvaluator_2d,
            RK2<FieldXY<CoordY>, DFieldXY>> const
            advection_y(
                    function_spline_y_interpolator,
                    adv_field_builder_y,
                    adv_field_spline_evaluator_y,
                    time_stepper_y);


    double const max_relative_error = AdvectionXY(advection_x, advection_y);
    EXPECT_LE(max_relative_error, 7.e-2);
    std::cout << "Test on " << x_size << "x" << y_size << "x" << vx_size << "x" << vy_size
              << " grid: max relative error = " << max_relative_error << std::endl;
}
