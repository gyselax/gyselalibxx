// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

//! [includes]
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include <Kokkos_Core.hpp>
//! [includes]

//! [X-dimension]
/// Our first continuous dimension
struct X;
//! [X-dimension]

//! [X-discretization]
/// A uniform discretization of X
struct DDimX : ddc::UniformPointSampling<X>
{
};
//! [X-discretization]

struct DDimFx : ddc::PeriodicSampling<ddc::Fourier<X>>
{
};

//! [Y-space]
// Our second continuous dimension
struct Y;
// Its uniform discretization
struct DDimY : ddc::UniformPointSampling<Y>
{
};

struct DDimFy : ddc::PeriodicSampling<ddc::Fourier<Y>>
{
};
//! [Y-space]

//! [time-space]
// Our simulated time dimension
struct T;
// Its uniform discretization
struct DDimT : ddc::UniformPointSampling<T>
{
};
//! [time-space]

//! [display]
/** A function to pretty print the temperature
 * @param time the time at which the output is made
 * @param temp the temperature at this time-step
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
}
//! [display]

//! [main-start]
int main(int argc, char** argv)
{
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
    size_t const nb_x_points = 100;
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
    std::cout << "Using spectral method \n";

    //! [X-global-domain]
    // Initialization of the global domain in X with gwx ghost points on
    // each side
    auto const x_domain
            = ddc::init_discrete_space<DDimX>(DDimX::init<DDimX>(
                    ddc::Coordinate<X>(x_start),
                    ddc::Coordinate<X>(x_end),
                    ddc::DiscreteVector<DDimX>(nb_x_points)));
    //! [X-global-domain]

    // Initialization of the global domain in Y with gwy ghost points on
    // each side
    auto const y_domain
            = ddc::init_discrete_space<DDimY>(DDimY::init<DDimY>(
                    ddc::Coordinate<Y>(y_start),
                    ddc::Coordinate<Y>(y_end),
                    ddc::DiscreteVector<DDimY>(nb_y_points)));

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
            2. / Kokkos::pow(Kokkos::numbers::pi, 2)
            / (kx * invdx2_max
               + ky * invdy2_max)}; // Classical stability theory gives .5 but empirically we see that for FFT method we need .2

    // number of time intervals required to reach the end time
    ddc::DiscreteVector<DDimT> const nb_time_steps {
            std::ceil((end_time - start_time) / max_dt) + .2};
    // Initialization of the global domain in time:
    // - the number of discrete time-points is equal to the number of
    //   steps + 1
    ddc::DiscreteDomain<DDimT> const time_domain
            = ddc::init_discrete_space<DDimT>(DDimT::init<DDimT>(
                    ddc::Coordinate<T>(start_time),
                    ddc::Coordinate<T>(end_time),
                    nb_time_steps + 1));

    //! [data allocation]
    // Maps temperature into the full domain (including ghosts) twice:
    // - once for the last fully computed time-step
    ddc::Chunk _last_temp(
            "_last_temp",
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            ddc::DeviceAllocator<double>());

    // - once for time-step being computed
    ddc::Chunk _next_temp(
            "_next_temp",
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            ddc::DeviceAllocator<double>());
    //! [data allocation]

    //! [initial-conditions]
    ddc::ChunkSpan const initial_temp = _last_temp.span_view();
    // Initialize the temperature on the main domain
    ddc::DiscreteDomain<DDimX, DDimY> x_mesh
            = ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain);
    ddc::parallel_for_each(
            x_mesh,
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                double const x
                        = ddc::coordinate(ddc::select<DDimX>(ixy));
                double const y
                        = ddc::coordinate(ddc::select<DDimY>(ixy));
                initial_temp(ixy) = 9.999 * ((x * x + y * y) < 0.25);
            });
    //! [initial-conditions]

    ddc::Chunk _host_temp(
            "_host_temp",
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            ddc::HostAllocator<double>());


    //! [initial output]
    // display the initial data
    ddc::parallel_deepcopy(_host_temp, _last_temp);
    display(ddc::coordinate(time_domain.front()),
            _host_temp[x_domain][y_domain]);
    // time of the iteration where the last output happened
    ddc::DiscreteElement<DDimT> last_output = time_domain.front();
    //! [initial output]

    ddc::init_discrete_space<DDimFx>(ddc::init_fourier_space<DDimFx>(
            ddc::DiscreteDomain<DDimX>(initial_temp.domain())));
    ddc::init_discrete_space<DDimFy>(ddc::init_fourier_space<DDimFy>(
            ddc::DiscreteDomain<DDimY>(initial_temp.domain())));
    ddc::DiscreteDomain<DDimFx, DDimFy> const k_mesh = ddc::
            FourierMesh<DDimFx, DDimFy>(initial_temp.domain(), false);
    ddc::Chunk Ff_allocation = ddc::
            Chunk("Ff_allocation",
                  k_mesh,
                  ddc::DeviceAllocator<Kokkos::complex<double>>());
    ddc::ChunkSpan Ff = Ff_allocation.span_view();

    //! [time iteration]
    for (auto const iter :
         time_domain.remove_first(ddc::DiscreteVector<DDimT>(1))) {
        //! [time iteration]

        //! [boundary conditions]
        //! [boundary conditions]

        //! [manipulated views]
        // a span excluding ghosts of the temperature at the time-step we
        // will build
        ddc::ChunkSpan const next_temp {_next_temp.span_view()};
        // a read-only view of the temperature at the previous time-step
        ddc::ChunkSpan const last_temp {_last_temp.span_view()};
        //! [manipulated views]

        //! [numerical scheme]
        // Stencil computation on the main domain
        ddc::FFT_Normalization norm = ddc::FFT_Normalization::BACKWARD;
        ddc::fft(Kokkos::DefaultExecutionSpace(), Ff, last_temp, {norm});
        ddc::parallel_for_each(
                k_mesh,
                KOKKOS_LAMBDA(ddc::DiscreteElement<DDimFx, DDimFy> const
                                      ikxky) {
                    ddc::DiscreteElement<DDimFx> const ikx
                            = ddc::select<DDimFx>(ikxky);
                    ddc::DiscreteElement<DDimFy> const iky
                            = ddc::select<DDimFy>(ikxky);
                    Ff(ikx, iky)
                            = Ff(ikx, iky)
                              * (1
                                 - (coordinate(ikx) * coordinate(ikx)
                                            * kx
                                    + coordinate(iky) * coordinate(iky)
                                              * ky)
                                           * max_dt); // Ff(t+dt) = (1-D*k^2*dt)*Ff(t)
                });
        ddc::
                ifft(Kokkos::DefaultExecutionSpace(),
                     next_temp,
                     Ff,
                     {norm});
        //! [numerical scheme]

        //! [output]
        if (iter - last_output >= t_output_period) {
            last_output = iter;
            ddc::parallel_deepcopy(_host_temp, _last_temp);
            display(ddc::coordinate(iter),
                    _host_temp[x_domain][y_domain]);
        }
        //! [output]

        //! [swap]
        // Swap our two buffers
        std::swap(_last_temp, _next_temp);
        //! [swap]
    }

    //! [final output]
    if (last_output < time_domain.back()) {
        ddc::parallel_deepcopy(_host_temp, _last_temp);
        display(ddc::coordinate(time_domain.back()),
                _host_temp[x_domain][y_domain]);
    }
    //! [final output]
}
