// SPDX-License-Identifier: MIT

#include <cmath>

#include <ddc/Chunk>
#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/PdiEvent>
#include <ddc/UniformDiscretization>
#include <ddc/for_each>

// Name of the axis
struct X;
struct Y;

using DDimX = UniformDiscretization<X>;
using DDimY = UniformDiscretization<Y>;

static unsigned nt = 10;
static unsigned nx = 100;
static unsigned ny = 200;
static unsigned gw = 1;
static double kx = 100.;
static double ky = 1.;
static double cfl = 0.99;

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  ghostwidth: int
  iter : int

data:
  temperature_extents: { type: array, subtype: int64, size: 2 }
  temperature:
    type: array
    subtype: double
    size: [ '$temperature_extents[0]', '$temperature_extents[1]' ]
    start: [ '$ghostwidth', '$ghostwidth' ]
    subsize: [ '$temperature_extents[0]-2*$ghostwidth', '$temperature_extents[1]-2*$ghostwidth' ]

plugins:
  decl_hdf5:
    - file: 'temperature_${iter:04}.h5'
      on_event: temperature
      collision_policy: replace_and_warn
      write: [temperature]
  trace: ~
)PDI_CFG";

int main()
{
    //! [mesh]
    // Origin on X
    Coordinate<X> const min_x(-1.);

    // Sampling step on X
    Coordinate<X> const dx(0.02);

    // Actual mesh on X
    init_discretization<DDimX>(min_x, dx);

    // Origin on Y
    Coordinate<Y> const min_y(-1.);

    // Sampling step on Y
    Coordinate<Y> const dy(0.01);

    // Actual mesh on Y
    init_discretization<DDimY>(min_y, dy);

    // Two-dimensional mesh on X,Y
    //! [mesh]

    //! [domain]
    // Take (nx+2gw) x (ny+2gw) points of `mesh_xy` starting from (0,0)
    DiscreteDomain<DDimX, DDimY> const domain_xy(
            DiscreteVector<DDimX, DDimY>(nx + 2 * gw, ny + 2 * gw));

    // Take only the inner domain (i.e. without ghost zone)
    DiscreteDomain<DDimX, DDimY> const inner_xy(
            DiscreteCoordinate<DDimX, DDimY>(gw, gw),
            DiscreteVector<DDimX, DDimY>(nx, ny));
    //! [domain]

    // Allocate data located at each point of `domain_xy` (including ghost region)
    //! [memory allocation]
    Chunk<double, DiscreteDomain<DDimX, DDimY>> T_in(domain_xy);
    Chunk<double, DiscreteDomain<DDimX, DDimY>> T_out(domain_xy);
    //! [memory allocation]

    //! [subdomains]
    // Ghost borders
    ChunkSpan const temperature_g_x_left = T_in[DiscreteCoordinate<DDimX>(gw - 1)];
    ChunkSpan const temperature_g_x_right = T_in[DiscreteCoordinate<DDimX>(nx + 2 * gw - 1)];
    ChunkSpan const temperature_g_y_left = T_in[DiscreteCoordinate<DDimY>(gw - 1)];
    ChunkSpan const temperature_g_y_right = T_in[DiscreteCoordinate<DDimY>(ny + 2 * gw - 1)];

    // Inner borders
    ChunkSpan const temperature_i_x_left = std::as_const(T_in)[DiscreteCoordinate<DDimX>(gw)];
    ChunkSpan const temperature_i_x_right = std::as_const(T_in)[DiscreteCoordinate<DDimX>(nx + gw)];
    ChunkSpan const temperature_i_y_left = std::as_const(T_in)[DiscreteCoordinate<DDimY>(gw)];
    ChunkSpan const temperature_i_y_right = std::as_const(T_in)[DiscreteCoordinate<DDimY>(ny + gw)];
    //! [subdomains]

    // Initialize the whole domain
    for_each(domain_xy, [&](DiscreteCoordinate<DDimX, DDimY> const ixy) {
        double const x = to_real(select<DDimX>(ixy));
        double const y = to_real(select<DDimY>(ixy));
        T_in(ixy) = 0.75 * ((x * x + y * y) < 0.25);
    });

    PDI_init(PC_parse_string(PDI_CFG));
    PDI_expose("ghostwidth", &gw, PDI_OUT);

    // Some heuristic for the time step
    double invdx2_max = 0.0;
    for (auto ix : select<DDimX>(inner_xy)) {
        invdx2_max = std::fmax(invdx2_max, 1.0 / (distance_at_left(ix) * distance_at_right(ix)));
    }
    double invdy2_max = 0.0;
    for (auto iy : select<DDimY>(inner_xy)) {
        invdy2_max = std::fmax(invdy2_max, 1.0 / (distance_at_left(iy) * distance_at_right(iy)));
    }
    double const dt = 0.5 * cfl / (kx * invdx2_max + ky * invdy2_max);
    std::size_t iter = 0;
    for (; iter < nt; ++iter) {
        //! [io/pdi]
        PdiEvent("temperature").with("iter", iter).and_with("temperature", T_in);
        //! [io/pdi]

        //! [numerical scheme]
        // Periodic boundary conditions
        deepcopy(temperature_g_x_left, temperature_i_x_right);
        deepcopy(temperature_g_x_right, temperature_i_x_left);
        deepcopy(temperature_g_y_left, temperature_i_y_right);
        deepcopy(temperature_g_y_right, temperature_i_y_left);

        // Stencil computation on inner domain `inner_xy`
        for_each(inner_xy, [&](DiscreteCoordinate<DDimX, DDimY> const ixy) {
            DiscreteCoordinate<DDimX> const ix = select<DDimX>(ixy);
            DiscreteCoordinate<DDimY> const iy = select<DDimY>(ixy);
            double const dx_l = distance_at_left(ix);
            double const dx_r = distance_at_right(ix);
            double const dx_m = 0.5 * (dx_l + dx_r);
            double const dy_l = distance_at_left(iy);
            double const dy_r = distance_at_right(iy);
            double const dy_m = 0.5 * (dy_l + dy_r);
            T_out(ix, iy) = T_in(ix, iy);
            T_out(ix, iy) += kx * dt
                             * (dx_l * T_in(ix + 1, iy) - 2.0 * dx_m * T_in(ix, iy)
                                + dx_r * T_in(ix - 1, iy))
                             / (dx_l * dx_m * dx_r);
            T_out(ix, iy) += ky * dt
                             * (dy_l * T_in(ix, iy + 1) - 2.0 * dy_m * T_in(ix, iy)
                                + dy_r * T_in(ix, iy - 1))
                             / (dy_l * dy_m * dy_r);
        });
        //! [numerical scheme]

        // Copy buf2 into buf1, a swap could also do the job
        deepcopy(T_in, T_out);
    }

    PdiEvent("temperature").with("iter", iter).and_with("temperature", T_in);

    PDI_finalize();

    return 0;
}
