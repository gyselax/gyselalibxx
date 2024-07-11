// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

//! [includes]
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include <ddc/ddc.hpp>
//! [includes]

//! [vector_generator]
std::vector<double> generate_random_vector(
        int n,
        double lower_bound,
        double higher_bound)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double>
            dis(lower_bound, higher_bound);

    std::vector<double> vec(n);
    vec[0] = lower_bound;
    vec[n - 1] = higher_bound;

    for (int i = 1; i < vec.size() - 1; ++i) {
        vec[i] = dis(gen);
    }

    std::sort(vec.begin(), vec.end());
    return vec;
}
//! [vector_generator]

//! [X-dimension]
struct X;
//! [X-dimension]

//! [X-discretization]
struct DDimX : ddc::NonUniformPointSampling<X>
{
};
//! [X-discretization]

//! [Y-space]
struct Y;
struct DDimY : ddc::NonUniformPointSampling<Y>
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

    //! [iterator_main-domain]
    std::vector<double> x_domain_vect
            = generate_random_vector(nb_x_points, x_start, x_end);
    //! [iterator_main-domain]

    std::size_t size_x = x_domain_vect.size();

    //! [ghost_points_x]
    std::vector<double> x_pre_ghost_vect {
            x_domain_vect.front()
            - (x_domain_vect.back() - x_domain_vect[size_x - 2])};

    std::vector<double> x_post_ghost_vect {
            x_domain_vect.back()
            + (x_domain_vect[1] - x_domain_vect.front())};
    //! [ghost_points_x]

    //! [build-domains]
    auto const [x_domain, ghosted_x_domain, x_pre_ghost, x_post_ghost]
            = ddc::init_discrete_space<DDimX>(DDimX::init_ghosted<DDimX>(
                    x_domain_vect,
                    x_pre_ghost_vect,
                    x_post_ghost_vect));
    //! [build-domains]

    ddc::DiscreteDomain<DDimX> const
            x_domain_begin(x_domain.front(), x_post_ghost.extents());
    ddc::DiscreteDomain<DDimX> const x_domain_end(
            x_domain.back() - x_pre_ghost.extents() + 1,
            x_pre_ghost.extents());

    //! [Y-vectors]
    std::vector<double> y_domain_vect
            = generate_random_vector(nb_y_points, y_start, y_end);

    std::size_t size_y = y_domain_vect.size();

    //! [ghost_points_y]
    std::vector<double> y_pre_ghost_vect {
            y_domain_vect.front()
            - (y_domain_vect.back() - y_domain_vect[size_y - 2])};
    std::vector<double> y_post_ghost_vect {
            y_domain_vect.back()
            + (y_domain_vect[1] - y_domain_vect.front())};
    //! [ghost_points_y]

    //! [Y-vectors]

    //! [build-Y-domain]
    auto const [y_domain, ghosted_y_domain, y_pre_ghost, y_post_ghost]
            = ddc::init_discrete_space<DDimY>(DDimY::init_ghosted<DDimY>(
                    y_domain_vect,
                    y_pre_ghost_vect,
                    y_post_ghost_vect));
    //! [build-Y-domain]

    ddc::DiscreteDomain<DDimY> const
            y_domain_begin(y_domain.front(), y_post_ghost.extents());

    ddc::DiscreteDomain<DDimY> const y_domain_end(
            y_domain.back() - y_pre_ghost.extents() + 1,
            y_pre_ghost.extents());

    //! [CFL-condition]

    double const invdx2_max = ddc::transform_reduce(
            x_domain,
            0.,
            ddc::reducer::max<double>(),
            [](ddc::DiscreteElement<DDimX> ix) {
                return 1.
                       / (ddc::distance_at_left(ix)
                          * ddc::distance_at_right(ix));
            });

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

    //! [CFL-condition]

    //! [time-domain]
    ddc::DiscreteVector<DDimT> const nb_time_steps(
            std::ceil((end_time - start_time) / max_dt) + .2);

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
