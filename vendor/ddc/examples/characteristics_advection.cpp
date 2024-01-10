// SPDX-License-Identifier: MIT

//! [includes]
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <Kokkos_Core.hpp>
//! [includes]
static constexpr std::size_t s_degree_x = 3;

//! [X-dimension]
/// Our first continuous dimension
struct X
{
    static constexpr bool PERIODIC = true;
};
//! [X-dimension]

//! [X-discretization]
/// A uniform discretization of X
using BSplinesX = ddc::UniformBSplines<X, s_degree_x>;
using GrevillePoints = ddc::GrevilleInterpolationPoints<
        BSplinesX,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;
using DDimX = GrevillePoints::interpolation_mesh_type;
//! [X-discretization]

//! [Y-space]
// Our second continuous dimension
struct Y;
// Its uniform discretization
using DDimY = ddc::UniformPointSampling<Y>;
//! [Y-space]

//! [time-space]
// Our simulated time dimension
struct T;
// Its uniform discretization
using DDimT = ddc::UniformPointSampling<T>;
//! [time-space]

//! [display]
/** A function to pretty print the density 
 * @param time the time at which the output is made
 * @param density the density at this time-step
 */
template <class ChunkType>
void display(double time, ChunkType density)
{
    double const mean_density = ddc::transform_reduce(
                                        density.domain(),
                                        0.,
                                        ddc::reducer::sum<double>(),
                                        density)
                                / density.domain().size();
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "At t = " << time << ",\n";
    std::cout << "  * mean density  = " << mean_density << "\n";
    // take a slice in the middle of the box
    ddc::ChunkSpan density_slice = density
            [ddc::get_domain<DDimY>(density).front()
             + ddc::get_domain<DDimY>(density).size() / 2];
    std::cout << "  * density[y:"
              << ddc::get_domain<DDimY>(density).size() / 2 << "] = {";
    ddc::for_each(
            ddc::policies::serial_host,
            ddc::get_domain<DDimX>(density),
            [=](ddc::DiscreteElement<DDimX> const ix) {
                std::cout << std::setw(6) << density_slice(ix) << " ";
            });
    std::cout << " }" << std::endl;
}
//! [display]

//! [main-start]
int main(int argc, char** argv)
{
    ddc::ScopeGuard scope(argc, argv);

    // some parameters that would typically be read from some form of
    // configuration file in a more realistic code

    //! [parameters]
    // Start of the domain of interest in the X dimension
    double const x_start = -1.;
    // End of the domain of interest in the X dimension
    double const x_end = 1.;
    // Number of discretization points in the X dimension
    size_t const nb_x_points = 100;
    // Velocity along x dimension
    double const vx = .2;
    // Start of the domain of interest in the Y dimension
    double const y_start = -1.;
    // End of the domain of interest in the Y dimension
    double const y_end = 1.;
    // Number of discretization points in the Y dimension
    size_t const nb_y_points = 100;
    // Simulated time at which to start simulation
    double const start_time = 0.;
    // Simulated time to reach as target of the simulation
    double const end_time = 10.;
    // Number of time-steps between outputs
    ptrdiff_t const t_output_period = 10;
    // Maximum time-step
    ddc::Coordinate<T> const max_dt {0.1};
    //! [parameters]

    //! [main-start]

    //! [X-global-domain]
    // Initialization of the global domain in X
    ddc::init_discrete_space<BSplinesX>(
            ddc::Coordinate<X>(x_start),
            ddc::Coordinate<X>(x_end),
            nb_x_points);
    ddc::init_discrete_space<DDimX>(
            ddc::GrevilleInterpolationPoints<
                    BSplinesX,
                    ddc::BoundCond::PERIODIC,
                    ddc::BoundCond::PERIODIC>::get_sampling());

    auto const x_domain = ddc::GrevilleInterpolationPoints<
            BSplinesX,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC>::get_domain();
    //! [X-global-domain]
    // Initialization of the global domain in Y
    auto const y_domain = ddc::init_discrete_space(
            DDimY::
                    init(ddc::Coordinate<Y>(y_start),
                         ddc::Coordinate<Y>(y_end),
                         ddc::DiscreteVector<DDimY>(nb_y_points)));

    //! [time-domains]

    // number of time intervals required to reach the end time
    ddc::DiscreteVector<DDimT> const nb_time_steps {
            std::ceil((end_time - start_time) / max_dt) + .2};
    // Initialization of the global domain in time:
    // - the number of discrete time-points is equal to the number of
    //   steps + 1
    ddc::DiscreteDomain<DDimT> const time_domain
            = ddc::init_discrete_space(
                    DDimT::
                            init(ddc::Coordinate<T>(start_time),
                                 ddc::Coordinate<T>(end_time),
                                 nb_time_steps + 1));
    //! [time-domains]

    //! [data allocation]
    // Maps density into the full domain twice:
    // - once for the last fully computed time-step
    ddc::Chunk last_density_alloc(
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            ddc::DeviceAllocator<double>());

    // - once for time-step being computed
    ddc::Chunk next_density_alloc(
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            ddc::DeviceAllocator<double>());
    //! [data allocation]

    //! [initial-conditions]
    ddc::ChunkSpan const initial_density
            = last_density_alloc.span_view();
    // Initialize the density on the main domain
    ddc::DiscreteDomain<DDimX, DDimY> x_mesh
            = ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain);
    ddc::for_each(
            ddc::policies::parallel_device,
            x_mesh,
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                double const x
                        = ddc::coordinate(ddc::select<DDimX>(ixy));
                double const y
                        = ddc::coordinate(ddc::select<DDimY>(ixy));
                initial_density(ixy)
                        = 9.999
                          * Kokkos::exp(-(x * x + y * y) / 0.1 / 2);
                // initial_density(ixy) = 9.999 * ((x * x + y * y) < 0.25);
            });
    //! [initial-conditions]

    ddc::Chunk host_density_alloc(
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            ddc::HostAllocator<double>());


    //! [initial output]
    // display the initial data
    ddc::deepcopy(host_density_alloc, last_density_alloc);
    display(ddc::coordinate(time_domain.front()),
            host_density_alloc[x_domain][y_domain]);
    // time of the iteration where the last output happened
    ddc::DiscreteElement<DDimT> last_output = time_domain.front();
    //! [initial output]

    //! [instantiate solver]
    ddc::SplineBuilderBatched<
            ddc::SplineBuilder<
                    Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultExecutionSpace::memory_space,
                    BSplinesX,
                    DDimX,
                    ddc::BoundCond::PERIODIC,
                    ddc::BoundCond::PERIODIC>,
            DDimX,
            DDimY>
            spline_builder(x_mesh);
    ddc::SplineEvaluatorBatched<
            ddc::SplineEvaluator<
                    Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultExecutionSpace::memory_space,
                    BSplinesX,
                    DDimX>,
            DDimX,
            DDimY>
            spline_evaluator(
                    spline_builder.spline_domain(),
                    ddc::g_null_boundary<BSplinesX>,
                    ddc::g_null_boundary<BSplinesX>);
    //! [instantiate solver]

    //! [instantiate intermediate chunks]
    // Instantiate chunk of spline coefs to receive output of spline_builder
    ddc::Chunk coef_alloc(
            spline_builder.spline_domain(),
            ddc::DeviceAllocator<double>());
    ddc::ChunkSpan coef = coef_alloc.span_view();

    // Instantiate chunk to receive feet coords
    ddc::Chunk feet_coords_alloc(
            spline_builder.vals_domain(),
            ddc::DeviceAllocator<ddc::Coordinate<X>>());
    ddc::ChunkSpan feet_coords = feet_coords_alloc.span_view();
    //! [instantiate intermediate chunks]


    //! [time iteration]
    for (auto const iter :
         time_domain.remove_first(ddc::DiscreteVector<DDimT>(1))) {
        //! [time iteration]

        //! [manipulated views]
        // a span of the density at the time-step we
        // will build
        ddc::ChunkSpan const next_density {
                next_density_alloc.span_view()};
        // a read-only view of the density at the previous time-step
        ddc::ChunkSpan const last_density {
                last_density_alloc.span_view()};
        //! [manipulated views]

        //! [numerical scheme]
        // Stencil computation on the main domain
        // Find the coordinates of the characteristics feet
        ddc::for_each(
                ddc::policies::parallel_device,
                feet_coords.domain(),
                KOKKOS_LAMBDA(
                        ddc::DiscreteElement<DDimX, DDimY> const e) {
                    feet_coords(e)
                            = ddc::coordinate(ddc::select<DDimX>(e))
                              - ddc::Coordinate<X>(
                                      vx * ddc::step<DDimT>());
                });
        // Interpolate the values at feets on the grid
        spline_builder(coef, last_density);
        spline_evaluator(
                next_density,
                feet_coords.span_cview(),
                coef.span_cview());
        //! [numerical scheme]

        //! [output]
        if (iter - last_output >= t_output_period) {
            last_output = iter;
            ddc::deepcopy(host_density_alloc, last_density_alloc);
            display(ddc::coordinate(iter),
                    host_density_alloc[x_domain][y_domain]);
        }
        //! [output]

        //! [swap]
        // Swap our two buffers
        std::swap(last_density_alloc, next_density_alloc);
        //! [swap]
    }

    //! [final output]
    if (last_output < time_domain.back()) {
        ddc::deepcopy(host_density_alloc, last_density_alloc);
        display(ddc::coordinate(time_domain.back()),
                host_density_alloc[x_domain][y_domain]);
    }
    //! [final output]
}
