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
struct X;
//! [X-dimension]

//! [X-discretization]
struct DDimX : ddc::UniformPointSampling<X>
{
};
//! [X-discretization]

//! [Y-space]
struct Y;
struct DDimY : ddc::UniformPointSampling<Y>
{
};
//! [Y-space]

//! [time-space]
struct T;
struct DDimT : ddc::UniformPointSampling<T>
{
};
//! [time-space]

//! [display]

/** A function to pretty print the temperature
 * @tparam ChunkType The type of chunk span. This way the template parameters are avoided,
 *                   should be deduced by the compiler.
 * @param time The time at which the output is made.
 * @param temp The temperature at this time-step.
 */

//! [display]

template <class ChunkType>
void display(double time, ChunkType temp)
{
    double const mean_temp = ddc::transform_reduce(
                                     temp.domain(),
                                     0.,
                                     ddc::reducer::sum<double>(),
                                     temp)
                             / temp.domain().size();
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "At t = " << time << ",\n";
    std::cout << "  * mean temperature  = " << mean_temp << "\n";
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
//! [main-start-x-parameters]
int main(int argc, char** argv)
{
#ifdef DDC_BUILD_PDI_WRAPPER
    auto pdi_conf = PC_parse_string("");
    PDI_init(pdi_conf);
#endif
    Kokkos::ScopeGuard const kokkos_scope(argc, argv);
    ddc::ScopeGuard const ddc_scope(argc, argv);


    //! [parameters]
    double const x_start = -1.;
    double const x_end = 1.;
    std::size_t const nb_x_points = 10;
    double const kx = .01;
    //! [main-start-x-parameters]
    //! [main-start-y-parameters]
    double const y_start = -1.;
    double const y_end = 1.;
    std::size_t const nb_y_points = 100;
    double const ky = .002;
    //! [main-start-y-parameters]
    //! [main-start-t-parameters]
    double const start_time = 0.;
    double const end_time = 10.;
    //! [main-start-t-parameters]
    std::ptrdiff_t const t_output_period = 10;
    //! [parameters]

    //! [main-start]

    //! [X-parameters]
    ddc::DiscreteVector<DDimX> const gwx(1);
    //! [X-parameters]

    //! [X-global-domain]
    auto const [x_domain, ghosted_x_domain, x_pre_ghost, x_post_ghost]
            = ddc::init_discrete_space<DDimX>(DDimX::init_ghosted<DDimX>(
                    ddc::Coordinate<X>(x_start),
                    ddc::Coordinate<X>(x_end),
                    ddc::DiscreteVector<DDimX>(nb_x_points),
                    gwx));
    //! [X-global-domain]

    //! [X-domains]
    ddc::DiscreteDomain<DDimX> const
            x_domain_begin(x_domain.front(), x_post_ghost.extents());
    ddc::DiscreteDomain<DDimX> const x_domain_end(
            x_domain.back() - x_pre_ghost.extents() + 1,
            x_pre_ghost.extents());
    //! [X-domains]

    //! [Y-domains]
    ddc::DiscreteVector<DDimY> const gwy(1);

    auto const [y_domain, ghosted_y_domain, y_pre_ghost, y_post_ghost]
            = ddc::init_discrete_space<DDimY>(DDimY::init_ghosted<DDimY>(
                    ddc::Coordinate<Y>(y_start),
                    ddc::Coordinate<Y>(y_end),
                    ddc::DiscreteVector<DDimY>(nb_y_points),
                    gwy));

    ddc::DiscreteDomain<DDimY> const
            y_domain_begin(y_domain.front(), y_post_ghost.extents());

    ddc::DiscreteDomain<DDimY> const y_domain_end(
            y_domain.back() - y_pre_ghost.extents() + 1,
            y_pre_ghost.extents());
    //! [Y-domains]

    //! [CFL-condition]

    double const dx = ddc::step<DDimX>();
    double const dy = ddc::step<DDimY>();
    double const invdx2 = 1. / (dx * dx);
    double const invdy2 = 1. / (dy * dy);

    ddc::Coordinate<T> const dt(.5 / (kx * invdx2 + ky * invdy2));

    //! [CFL-condition]

    //! [time-domain]
    ddc::DiscreteVector<DDimT> const nb_time_steps(
            std::ceil((end_time - start_time) / dt) + .2);

    ddc::DiscreteDomain<DDimT> const time_domain
            = ddc::init_discrete_space<DDimT>(DDimT::init<DDimT>(
                    ddc::Coordinate<T>(start_time),
                    ddc::Coordinate<T>(end_time),
                    nb_time_steps + 1));

    //! [time-domain]

    //! [data allocation]
    ddc::Chunk ghosted_last_temp(
            "ghosted_last_temp",
            ddc::DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            ddc::DeviceAllocator<double>());

    ddc::Chunk ghosted_next_temp(
            "ghosted_next_temp",
            ddc::DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            ddc::DeviceAllocator<double>());
    //! [data allocation]

    //! [initial-chunkspan]
    ddc::ChunkSpan const ghosted_initial_temp
            = ghosted_last_temp.span_view();
    //! [initial-chunkspan]

    //! [fill-initial-chunkspan]
    ddc::parallel_for_each(
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                double const x = ddc::coordinate(
                        ddc::DiscreteElement<DDimX>(ixy));
                double const y = ddc::coordinate(
                        ddc::DiscreteElement<DDimY>(ixy));
                ghosted_initial_temp(ixy)
                        = 9.999 * ((x * x + y * y) < 0.25);
            });
    //! [fill-initial-chunkspan]

    //! [host-chunk]
    ddc::Chunk ghosted_temp(
            "ghost_temp",
            ddc::DiscreteDomain<
                    DDimX,
                    DDimY>(ghosted_x_domain, ghosted_y_domain),
            ddc::HostAllocator<double>());
    //! [host-chunk]

    //! [initial-deepcopy]
    ddc::parallel_deepcopy(ghosted_temp, ghosted_last_temp);
    //! [initial-deepcopy]

    //! [initial-display]
    display(ddc::coordinate(time_domain.front()),
            ghosted_temp[x_domain][y_domain]);
    //! [initial-display]

    //! [last-output-iter]
    ddc::DiscreteElement<DDimT> last_output_iter = time_domain.front();
    //! [last-output-iter]

    //! [time iteration]
    for (ddc::DiscreteElement<DDimT> const iter :
         time_domain.remove_first(ddc::DiscreteVector<DDimT>(1))) {
        //! [time iteration]

        //! [boundary conditions]
        ddc::parallel_deepcopy(
                ghosted_last_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(x_pre_ghost, y_domain)],
                ghosted_last_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(y_domain, x_domain_end)]);
        ddc::parallel_deepcopy(
                ghosted_last_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(y_domain, x_post_ghost)],
                ghosted_last_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(y_domain, x_domain_begin)]);
        ddc::parallel_deepcopy(
                ghosted_last_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(x_domain, y_pre_ghost)],
                ghosted_last_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(x_domain, y_domain_end)]);
        ddc::parallel_deepcopy(
                ghosted_last_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(x_domain, y_post_ghost)],
                ghosted_last_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(x_domain, y_domain_begin)]);
        //! [boundary conditions]

        //! [manipulated views]
        ddc::ChunkSpan const next_temp(
                ghosted_next_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(x_domain, y_domain)]);
        ddc::ChunkSpan const last_temp(ghosted_last_temp.span_cview());
        //! [manipulated views]

        //! [numerical scheme]
        ddc::parallel_for_each(
                next_temp.domain(),
                KOKKOS_LAMBDA(
                        ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                    ddc::DiscreteElement<DDimX> const ix(ixy);
                    ddc::DiscreteElement<DDimY> const iy(ixy);
                    double const dt = ddc::step<DDimT>();

                    next_temp(ix, iy) = last_temp(ix, iy);
                    next_temp(ix, iy) += kx * dt
                                         * (last_temp(ix + 1, iy)
                                            - 2.0 * last_temp(ix, iy)
                                            + last_temp(ix - 1, iy))
                                         * invdx2;

                    next_temp(ix, iy) += ky * dt
                                         * (last_temp(ix, iy + 1)
                                            - 2.0 * last_temp(ix, iy)
                                            + last_temp(ix, iy - 1))
                                         * invdy2;
                });
        //! [numerical scheme]

        //! [output]
        if (iter - last_output_iter >= t_output_period) {
            last_output_iter = iter;
            ddc::parallel_deepcopy(ghosted_temp, ghosted_last_temp);
            display(ddc::coordinate(iter),
                    ghosted_temp[ddc::DiscreteDomain<
                            DDimX,
                            DDimY>(x_domain, y_domain)]);
        }
        //! [output]

        //! [swap]
        std::swap(ghosted_last_temp, ghosted_next_temp);
        //! [swap]
    }


    //! [final output]
    if (last_output_iter < time_domain.back()) {
        ddc::parallel_deepcopy(ghosted_temp, ghosted_last_temp);
        display(ddc::coordinate(time_domain.back()),
                ghosted_temp[ddc::DiscreteDomain<
                        DDimX,
                        DDimY>(x_domain, y_domain)]);
    }
    //! [final output]

#ifdef DDC_BUILD_PDI_WRAPPER
    PDI_finalize();
    PC_tree_destroy(&pdi_conf);
#endif
}
