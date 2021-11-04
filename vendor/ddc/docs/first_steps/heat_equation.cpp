#include <ddc/Chunk>
#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/PdiEvent>
#include <ddc/UniformDiscretization>

// Name of the axis
struct X;
struct Y;

using DDimX = UniformDiscretization<X>;
using DDimY = UniformDiscretization<Y>;

static unsigned nt = 10;
static unsigned nx = 100;
static unsigned ny = 200;
static unsigned gw = 1;
static double k = 0.1;

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
    DDimX const ddim_x(min_x, dx);

    // Origin on Y
    Coordinate<Y> const min_y(-1.);

    // Sampling step on Y
    Coordinate<Y> const dy(0.01);

    // Actual mesh on Y
    DDimY const ddim_y(min_y, dy);

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
    for (DiscreteCoordinate<DDimX> const ix : select<DDimX>(domain_xy)) {
        double const x = ddim_x.to_real(ix);
        for (DiscreteCoordinate<DDimY> const iy : select<DDimY>(domain_xy)) {
            double const y = ddim_y.to_real(iy);
            T_in(ix, iy) = 0.75 * ((x * x + y * y) < 0.25);
        }
    }

    PDI_init(PC_parse_string(PDI_CFG));
    PDI_expose("ghostwidth", &gw, PDI_OUT);

    double const dt = 0.49 / (1.0 / (dx * dx) + 1.0 / (dy * dy));
    double const Cx = k * dt / (dx * dx);
    double const Cy = k * dt / (dy * dy);
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
        for (DiscreteCoordinate<DDimX> const ix : select<DDimX>(inner_xy)) {
            for (DiscreteCoordinate<DDimY> const iy : select<DDimY>(inner_xy)) {
                T_out(ix, iy) = T_in(ix, iy);
                T_out(ix, iy) += Cx * (T_in(ix + 1, iy) - 2.0 * T_in(ix, iy) + T_in(ix - 1, iy));
                T_out(ix, iy) += Cy * (T_in(ix, iy + 1) - 2.0 * T_in(ix, iy) + T_in(ix, iy - 1));
            }
        }
        //! [numerical scheme]

        // Copy buf2 into buf1, a swap could also do the job
        deepcopy(T_in, T_out);
    }

    PdiEvent("temperature").with("iter", iter).and_with("temperature", T_in);

    PDI_finalize();

    return 0;
}
