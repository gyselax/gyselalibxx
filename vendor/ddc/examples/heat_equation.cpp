// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

//! [includes]
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

#include <ddc/ddc.hpp>
//! [includes]


//! [X-dimension]
/// Our first continuous dimension
struct X;
//! [X-dimension]

//! [X-discretization]
/// A uniform discretization of X
using DDimX = ddc::UniformPointSampling<X>;
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
/** A function to pretty print the temperature
 * @tparam ChunkType The type of chunk span. This way the template parameters are avoided,
 *                   should be deduced by the compiler.
 * @param time The time at which the output is made.
 * @param temp The temperature at this time-step.
 */
template <class ChunkType>
void display(double time, ChunkType temp)
{
    // For `Chunk`/`ChunkSpan`, the ()-operator is used to access stored values with a single `DiscreteElement` as
    // input. So it is used here as a function that maps indices of the temperature domain
    // to the temperature value at that point
    double const mean_temp = ddc::transform_reduce(
                                     temp.domain(),
                                     0.,
                                     ddc::reducer::sum<double>(),
                                     temp)
                             / temp.domain().size();
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "At t = " << time << ",\n";
    std::cout << "  * mean temperature  = " << mean_temp << "\n";
    // take a slice in the middle of the box
    ddc::ChunkSpan temp_slice
            = temp[ddc::get_domain<DDimY>(temp).front()
                   + ddc::get_domain<DDimY>(temp).size() / 2];
    std::cout << "  * temperature[y:"
              << ddc::get_domain<DDimY>(temp).size() / 2 << "] = {";
    ddc::for_each(
            ddc::get_domain<DDimX>(temp),
            [=](ddc::DiscreteElement<DDimX> const ix) {
                std::cout << std::setw(6) << temp_slice(ix);
            });
    std::cout << " }" << std::endl;

#ifdef DDC_BUILD_PDI_WRAPPER
    ddc::PdiEvent("display")
            .with("temp", temp)
            .and_with("mean_temp", mean_temp)
            .and_with("temp_slice", temp_slice);
#endif
}
//! [display]


//! [main-start]
int main(int argc, char** argv)
{
#ifdef DDC_BUILD_PDI_WRAPPER
    auto pdi_conf = PC_parse_string("");
    PDI_init(pdi_conf);
#endif
    Kokkos::ScopeGuard const kokkos_scope(argc, argv);
    ddc::ScopeGuard const ddc_scope(argc, argv);

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
    ptrdiff_t const t_output_period = 10;
    //! [parameters]

    //! [main-start]
    //! [X-parameters]
    // Number of ghost points to use on each side in X
    ddc::DiscreteVector<DDimX> static constexpr gwx {1};
    //! [X-parameters]

    //! [X-global-domain]
    // Initialization of the global domain in X with gwx ghost points on
    // each side
    auto const [x_domain, ghosted_x_domain, x_pre_ghost, x_post_ghost]
            = ddc::init_discrete_space(DDimX::init_ghosted(
                    ddc::Coordinate<X>(x_start),
                    ddc::Coordinate<X>(x_end),
                    ddc::DiscreteVector<DDimX>(nb_x_points),
                    gwx));
    //! [X-global-domain]

    //! [X-domains]
    // our zone at the start of the domain that will be mirrored to the
    // ghost
    ddc::DiscreteDomain<DDimX> const
            x_domain_begin(x_domain.front(), x_post_ghost.extents());
    // our zone at the end of the domain that will be mirrored to the
    // ghost
    ddc::DiscreteDomain<DDimX> const x_domain_end(
            x_domain.back() - x_pre_ghost.extents() + 1,
            x_pre_ghost.extents());
    //! [X-domains]

    //! [Y-domains]
    // Number of ghost points to use on each side in Y
    ddc::DiscreteVector<DDimY> static constexpr gwy {1};

    // Initialization of the global domain in Y with gwy ghost points on
    // each side
    auto const [y_domain, ghosted_y_domain, y_pre_ghost, y_post_ghost]
            = ddc::init_discrete_space(DDimY::init_ghosted(
                    ddc::Coordinate<Y>(y_start),
                    ddc::Coordinate<Y>(y_end),
                    ddc::DiscreteVector<DDimY>(nb_y_points),
                    gwy));

    // our zone at the start of the domain that will be mirrored to the
    // ghost
    ddc::DiscreteDomain<DDimY> const
            y_domain_begin(y_domain.front(), y_post_ghost.extents());
    // our zone at the end of the domain that will be mirrored to the
    // ghost
    ddc::DiscreteDomain<DDimY> const y_domain_end(
            y_domain.back() - y_pre_ghost.extents() + 1,
            y_pre_ghost.extents());
    //! [Y-domains]

    //! [time-domains]
    // max(1/dx^2)
    double const invdx2_max = ddc::transform_reduce(
            x_domain,
            0.,
            ddc::reducer::max<double>(),
            [](ddc::DiscreteElement<DDimX> ix) {
                return 1.
                       / (ddc::distance_at_left(ix)
                          * ddc::distance_at_right(ix));
            });
    // max(1/dy^2)
    double const invdy2_max = ddc::transform_reduce(
            y_domain,
            0.,
            ddc::reducer::max<double>(),
            [](ddc::DiscreteElement<DDimY> iy) {
                return 1.
                       / (ddc::distance_at_left(iy)
                          * ddc::distance_at_right(iy));
            });
    ddc::Coordinate<T> const max_dt {
            .5 / (kx * invdx2_max + ky * invdy2_max)};

    // number of time intervals required to reach the end time
    // The + .2 is used to make sure it is rounded to the correct number of steps
    ddc::DiscreteVector<DDimT> const nb_time_steps {
            std::ceil((end_time - start_time) / max_dt) + .2};
    // Initialization of the global domain in time:
    // - the number of discrete time-points is equal to the number of
    //   steps + 1
    // `init` takes required information to initialize the attributes of the dimension.
    ddc::DiscreteDomain<DDimT> const time_domain
            = ddc::init_discrete_space(
                    DDimT::
                            init(ddc::Coordinate<T>(start_time),
                                 ddc::Coordinate<T>(end_time),
                                 nb_time_steps + 1));
    //! [time-domains]

    //! [data allocation]
    // Maps temperature into the full domain (including ghosts) twice:
    // - once for the last fully computed time-step
    ddc::Chunk ghosted_last_temp(
            "ghosted_last_temp",
            ddc::DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            ddc::DeviceAllocator<double>());

    // - once for time-step being computed
    // The `DeviceAllocator` is responsible for allocating memory on the default memory space.
    ddc::Chunk ghosted_next_temp(
            "ghosted_next_temp",
            ddc::DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            ddc::DeviceAllocator<double>());
    //! [data allocation]

    //! [initial-conditions]
    // The const qualifier makes it clear that ghosted_initial_temp always references
    // the same chunk, `ghosted_last_temp` in this case
    ddc::ChunkSpan const ghosted_initial_temp
            = ghosted_last_temp.span_view();
    // Initialize the temperature on the main domain
    ddc::parallel_for_each(
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                double const x
                        = ddc::coordinate(ddc::select<DDimX>(ixy));
                double const y
                        = ddc::coordinate(ddc::select<DDimY>(ixy));
                ghosted_initial_temp(ixy)
                        = 9.999 * ((x * x + y * y) < 0.25);
            });
    //! [initial-conditions]

    ddc::Chunk ghosted_temp(
            "ghost_temp",
            ddc::DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            ddc::HostAllocator<double>());


    //! [initial output]
    // display the initial data
    ddc::parallel_deepcopy(ghosted_temp, ghosted_last_temp);
    display(ddc::coordinate(time_domain.front()),
            ghosted_temp[x_domain][y_domain]);
    // time of the iteration where the last output happened
    ddc::DiscreteElement<DDimT> last_output_iter = time_domain.front();
    //! [initial output]

    //! [time iteration]
    for (auto const iter :
         time_domain.remove_first(ddc::DiscreteVector<DDimT>(1))) {
        //! [time iteration]

        //! [boundary conditions]
        // Periodic boundary conditions
        ddc::parallel_deepcopy(
                ghosted_last_temp[x_pre_ghost][y_domain],
                ghosted_last_temp[y_domain][x_domain_end]);
        ddc::parallel_deepcopy(
                ghosted_last_temp[y_domain][x_post_ghost],
                ghosted_last_temp[y_domain][x_domain_begin]);
        ddc::parallel_deepcopy(
                ghosted_last_temp[x_domain][y_pre_ghost],
                ghosted_last_temp[x_domain][y_domain_end]);
        ddc::parallel_deepcopy(
                ghosted_last_temp[x_domain][y_post_ghost],
                ghosted_last_temp[x_domain][y_domain_begin]);
        //! [boundary conditions]

        //! [manipulated views]
        // a span excluding ghosts of the temperature at the time-step we
        // will build
        ddc::ChunkSpan const next_temp {
                ghosted_next_temp[x_domain][y_domain]};
        // a read-only view of the temperature at the previous time-step
        // span_cview returns a read-only ChunkSpan
        ddc::ChunkSpan const last_temp {ghosted_last_temp.span_cview()};
        //! [manipulated views]

        //! [numerical scheme]
        // Stencil computation on the main domain
        ddc::parallel_for_each(
                next_temp.domain(),
                KOKKOS_LAMBDA(
                        ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                    ddc::DiscreteElement<DDimX> const ix
                            = ddc::select<DDimX>(ixy);
                    ddc::DiscreteElement<DDimY> const iy
                            = ddc::select<DDimY>(ixy);
                    double const dx_l = ddc::distance_at_left(ix);
                    double const dx_r = ddc::distance_at_right(ix);
                    double const dx_m = 0.5 * (dx_l + dx_r);
                    double const dy_l = ddc::distance_at_left(iy);
                    double const dy_r = ddc::distance_at_right(iy);
                    double const dy_m = 0.5 * (dy_l + dy_r);
                    next_temp(ix, iy) = last_temp(ix, iy);
                    next_temp(ix, iy)
                            += kx * ddc::step<DDimT>()
                               * (dx_l * last_temp(ix + 1, iy)
                                  - 2.0 * dx_m * last_temp(ix, iy)
                                  + dx_r * last_temp(ix - 1, iy))
                               / (dx_l * dx_m * dx_r);
                    next_temp(ix, iy)
                            += ky * ddc::step<DDimT>()
                               * (dy_l * last_temp(ix, iy + 1)
                                  - 2.0 * dy_m * last_temp(ix, iy)
                                  + dy_r * last_temp(ix, iy - 1))
                               / (dy_l * dy_m * dy_r);
                });
        //! [numerical scheme]

        //! [output]
        if (iter - last_output_iter >= t_output_period) {
            last_output_iter = iter;
            ddc::parallel_deepcopy(ghosted_temp, ghosted_last_temp);
            display(ddc::coordinate(iter),
                    ghosted_temp[x_domain][y_domain]);
        }
        //! [output]

        //! [swap]
        // Swap our two buffers
        std::swap(ghosted_last_temp, ghosted_next_temp);
        //! [swap]
    }

    //! [final output]
    if (last_output_iter < time_domain.back()) {
        ddc::parallel_deepcopy(ghosted_temp, ghosted_last_temp);
        display(ddc::coordinate(time_domain.back()),
                ghosted_temp[x_domain][y_domain]);
    }
    //! [final output]

#ifdef DDC_BUILD_PDI_WRAPPER
    PDI_finalize();
    PC_tree_destroy(&pdi_conf);
#endif
}
