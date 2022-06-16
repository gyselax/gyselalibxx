// SPDX-License-Identifier: MIT

//! [includes]
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

#include <ddc/ddc.hpp>

#include <Kokkos_Core.hpp>
//! [includes]


//! [X-dimension]
/// Our first continuous dimension
struct X;
//! [X-dimension]

//! [X-discretization]
/// A uniform discretization of X
using DDimX = UniformDiscretization<X>;
//! [X-discretization]

//! [Y-space]
// Our second continuous dimension
struct Y;
// Its uniform discretization
using DDimY = UniformDiscretization<Y>;
//! [Y-space]

//! [time-space]
// Our simulated time dimension
struct T;
// Its uniform discretization
using DDimT = UniformDiscretization<T>;
//! [time-space]


//! [display]
/** A function to pretty print the temperature
 * @param time the time at which the output is made
 * @param temp the temperature at this time-step
 */
template <class ChunkType>
void display(double time, ChunkType temp)
{
    double const mean_temp = transform_reduce(
                                     temp.domain(),
                                     0.,
                                     reducer::sum<double>(),
                                     temp)
                             / temp.domain().size();
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "At t = " << time << ",\n";
    std::cout << "  * mean temperature  = " << mean_temp << "\n";
    // take a slice in the middle of the box
    ChunkSpan temp_slice
            = temp[get_domain<DDimY>(temp).front()
                   + get_domain<DDimY>(temp).size() / 2];
    std::cout << "  * temperature[y:"
              << get_domain<DDimY>(temp).size() / 2 << "] = {";
    for_each(
            policies::serial_host,
            get_domain<DDimX>(temp),
            [=](DiscreteCoordinate<DDimX> const ix) {
                std::cout << std::setw(6) << temp_slice(ix);
            });
    std::cout << " }" << std::endl;
}
//! [display]


//! [main-start]
int main(int argc, char** argv)
{
    ScopeGuard scope(argc, argv);

    // some parameters that would typically be read from some form of
    // configuration file in a more realistic code

    //! [parameters]
    // Start of the domain of interest in the X dimension
    double const x_start = -1.;
    // End of the domain of interest in the X dimension
    double const x_end = 1.;
    // Number of discretization points in the X dimension
    size_t const nb_x_points = 10;
    // Thermal diffusion coefficient
    double const kx = .01;
    // Start of the domain of interest in the Y dimension
    double const y_start = -1.;
    // End of the domain of interest in the Y dimension
    double const y_end = 1.;
    // Number of discretization points in the Y dimension
    size_t const nb_y_points = 100;
    // Thermal diffusion coefficient
    double const ky = .002;
    // Simulated time at which to start simulation
    double const start_time = 0.;
    // Simulated time to reach as target of the simulation
    double const end_time = 10.;
    // Number of time-steps between outputs
    size_t const t_output_period = 10;
    //! [parameters]

    //! [main-start]
    //! [X-parameters]
    // Number of ghost points to use on each side in X
    DiscreteVector<DDimX> static constexpr gwx {1};
    //! [X-parameters]

    //! [X-global-domain]
    // Initialization of the global domain in X with gwx ghost points on
    // each side
    auto const [x_domain, ghosted_x_domain, x_pre_ghost, x_post_ghost]
            = init_discretization(DDimX::init_ghosted(
                    Coordinate<X>(x_start),
                    Coordinate<X>(x_end),
                    DiscreteVector<DDimX>(nb_x_points),
                    gwx));
    //! [X-global-domain]

    //! [X-domains]
    // our zone at the start of the domain that will be mirrored to the
    // ghost
    DiscreteDomain const
            x_domain_begin(x_domain.front(), x_post_ghost.extents());
    // our zone at the end of the domain that will be mirrored to the
    // ghost
    DiscreteDomain const x_domain_end(
            x_domain.back() - x_pre_ghost.extents() + 1,
            x_pre_ghost.extents());
    //! [X-domains]

    //! [Y-domains]
    // Number of ghost points to use on each side in Y
    DiscreteVector<DDimY> static constexpr gwy {1};

    // Initialization of the global domain in Y with gwy ghost points on
    // each side
    auto const [y_domain, ghosted_y_domain, y_pre_ghost, y_post_ghost]
            = init_discretization(DDimY::init_ghosted(
                    Coordinate<Y>(y_start),
                    Coordinate<Y>(y_end),
                    DiscreteVector<DDimY>(nb_y_points),
                    gwy));

    // our zone at the start of the domain that will be mirrored to the
    // ghost
    DiscreteDomain const
            y_domain_begin(y_domain.front(), y_post_ghost.extents());
    // our zone at the end of the domain that will be mirrored to the
    // ghost
    DiscreteDomain const y_domain_end(
            y_domain.back() - y_pre_ghost.extents() + 1,
            y_pre_ghost.extents());
    //! [Y-domains]

    //! [time-domains]
    // max(1/dx^2)
    double const invdx2_max = transform_reduce(
            x_domain,
            0.,
            reducer::max<double>(),
            [](DiscreteCoordinate<DDimX> ix) {
                return 1.
                       / (distance_at_left(ix) * distance_at_right(ix));
            });
    // max(1/dy^2)
    double const invdy2_max = transform_reduce(
            y_domain,
            0.,
            reducer::max<double>(),
            [](DiscreteCoordinate<DDimY> iy) {
                return 1.
                       / (distance_at_left(iy) * distance_at_right(iy));
            });
    Coordinate<T> const max_dt {
            .5 / (kx * invdx2_max + ky * invdy2_max)};

    // number of time intervals required to reach the end time
    DiscreteVector<DDimT> const nb_time_steps {
            std::ceil((end_time - start_time) / max_dt) + .2};
    // Initialization of the global domain in time:
    // - the number of discrete time-points is equal to the number of
    //   steps + 1
    DiscreteDomain<DDimT> const time_domain
            = init_discretization(DDimT::
                                          init(Coordinate<T>(start_time),
                                               Coordinate<T>(end_time),
                                               nb_time_steps + 1));
    //! [time-domains]

    //! [data allocation]
    // Maps temperature into the full domain (including ghosts) twice:
    // - once for the last fully computed time-step
    Chunk ghosted_last_temp(
            DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            DeviceAllocator<double>());

    // - once for time-step being computed
    Chunk ghosted_next_temp(
            DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            DeviceAllocator<double>());
    //! [data allocation]

    //! [initial-conditions]
    ChunkSpan const ghosted_initial_temp = ghosted_last_temp.span_view();
    // Initialize the temperature on the main domain
    for_each(
            policies::parallel_device,
            DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            DDC_LAMBDA(DiscreteCoordinate<DDimX, DDimY> const ixy) {
                double const x = to_real(select<DDimX>(ixy));
                double const y = to_real(select<DDimY>(ixy));
                ghosted_initial_temp(ixy)
                        = 9.999 * ((x * x + y * y) < 0.25);
            });
    //! [initial-conditions]

    Chunk ghosted_temp(
            DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            HostAllocator<double>());


    //! [initial output]
    // display the initial data
    deepcopy(ghosted_temp, ghosted_last_temp);
    display(to_real(time_domain.front()),
            ghosted_temp[x_domain][y_domain]);
    // time of the iteration where the last output happened
    DiscreteCoordinate<DDimT> last_output = time_domain.front();
    //! [initial output]

    //! [time iteration]
    for (auto const iter :
         time_domain.remove_first(DiscreteVector<DDimT>(1))) {
        //! [time iteration]

        //! [boundary conditions]
        // Periodic boundary conditions
        deepcopy(
                ghosted_last_temp[x_pre_ghost][y_domain],
                ghosted_last_temp[y_domain][x_domain_end]);
        deepcopy(
                ghosted_last_temp[y_domain][x_post_ghost],
                ghosted_last_temp[y_domain][x_domain_begin]);
        deepcopy(
                ghosted_last_temp[x_domain][y_pre_ghost],
                ghosted_last_temp[x_domain][y_domain_end]);
        deepcopy(
                ghosted_last_temp[x_domain][y_post_ghost],
                ghosted_last_temp[x_domain][y_domain_begin]);
        //! [boundary conditions]

        //! [manipulated views]
        // a span excluding ghosts of the temperature at the time-step we
        // will build
        ChunkSpan const next_temp {
                ghosted_next_temp[x_domain][y_domain]};
        // a read-only view of the temperature at the previous time-step
        ChunkSpan const last_temp {ghosted_last_temp.span_view()};
        //! [manipulated views]

        //! [numerical scheme]
        // Stencil computation on the main domain
        for_each(
                policies::parallel_device,
                next_temp.domain(),
                DDC_LAMBDA(DiscreteCoordinate<DDimX, DDimY> const ixy) {
                    DiscreteCoordinate<DDimX> const ix
                            = select<DDimX>(ixy);
                    DiscreteCoordinate<DDimY> const iy
                            = select<DDimY>(ixy);
                    double const dx_l = distance_at_left(ix);
                    double const dx_r = distance_at_right(ix);
                    double const dx_m = 0.5 * (dx_l + dx_r);
                    double const dy_l = distance_at_left(iy);
                    double const dy_r = distance_at_right(iy);
                    double const dy_m = 0.5 * (dy_l + dy_r);
                    next_temp(ix, iy) = last_temp(ix, iy);
                    next_temp(ix, iy)
                            += kx * max_dt
                               * (dx_l * last_temp(ix + 1, iy)
                                  - 2.0 * dx_m * last_temp(ix, iy)
                                  + dx_r * last_temp(ix - 1, iy))
                               / (dx_l * dx_m * dx_r);
                    next_temp(ix, iy)
                            += ky * max_dt
                               * (dy_l * last_temp(ix, iy + 1)
                                  - 2.0 * dy_m * last_temp(ix, iy)
                                  + dy_r * last_temp(ix, iy - 1))
                               / (dy_l * dy_m * dy_r);
                });
        //! [numerical scheme]

        //! [output]
        if (iter - last_output >= t_output_period) {
            last_output = iter;
            deepcopy(ghosted_temp, ghosted_last_temp);
            display(to_real(iter), ghosted_temp[x_domain][y_domain]);
        }
        //! [output]

        //! [swap]
        // Swap our two buffers
        std::swap(ghosted_last_temp, ghosted_next_temp);
        //! [swap]
    }

    //! [final output]
    if (last_output < time_domain.back()) {
        deepcopy(ghosted_temp, ghosted_last_temp);
        display(to_real(time_domain.back()),
                ghosted_temp[x_domain][y_domain]);
    }
    //! [final output]
}
