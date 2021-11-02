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

static int nt = 10;
static int nx = 100;
static int ny = 200;
static int gw = 1;
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

int main(int argc, char** argv)
{
    //! [mesh]
    // Origin on X
    Coordinate<X> min_x(-1.);

    // Sampling step on X
    Coordinate<X> dx(0.02);

    // Actual mesh on X
    DDimX ddim_x(min_x, dx);

    // Origin on Y
    Coordinate<Y> min_y(-1.);

    // Sampling step on Y
    Coordinate<Y> dy(0.01);

    // Actual mesh on Y
    DDimY ddim_y(min_y, dy);

    // Two-dimensional mesh on X,Y
    //! [mesh]

    //! [domain]
    // Take (nx+2gw) x (ny+2gw) points of `mesh_xy` starting from (0,0)
    DiscreteDomain<DDimX, DDimY> domain_xy(DiscreteVector<DDimX, DDimY>(nx + 2 * gw, ny + 2 * gw));

    // Take only the inner domain (i.e. without ghost zone)
    DiscreteDomain<DDimX, DDimY> inner_xy(
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
    auto temperature_g_x_left = T_in[DiscreteCoordinate<DDimX>(gw - 1)];
    auto temperature_g_x_right = T_in[DiscreteCoordinate<DDimX>(nx + 2 * gw - 1)];
    auto temperature_g_y_left = T_in[DiscreteCoordinate<DDimY>(gw - 1)];
    auto temperature_g_y_right = T_in[DiscreteCoordinate<DDimY>(ny + 2 * gw - 1)];

    // Inner borders
    auto temperature_i_x_left = T_in[DiscreteCoordinate<DDimX>(gw)];
    auto temperature_i_x_right = T_in[DiscreteCoordinate<DDimX>(nx + gw)];
    auto temperature_i_y_left = T_in[DiscreteCoordinate<DDimY>(gw)];
    auto temperature_i_y_right = T_in[DiscreteCoordinate<DDimY>(ny + gw)];
    //! [subdomains]

    // Initialize the whole domain
    for (DiscreteCoordinate<DDimX> ix : select<DDimX>(domain_xy)) {
        double x = ddim_x.to_real(ix);
        for (DiscreteCoordinate<DDimY> iy : select<DDimY>(domain_xy)) {
            double y = ddim_y.to_real(iy);
            T_in(ix, iy) = 0.75 * ((x * x + y * y) < 0.25);
        }
    }

    PDI_init(PC_parse_string(PDI_CFG));
    PDI_expose("ghostwidth", &gw, PDI_OUT);

    const double dt = 0.49 / (1.0 / (dx * dx) + 1.0 / (dy * dy));
    const double Cx = k * dt / (dx * dx);
    const double Cy = k * dt / (dy * dy);
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
        for (DiscreteCoordinate<DDimX> ix : select<DDimX>(inner_xy)) {
            const DiscreteCoordinate<DDimX> ixp(ix + 1);
            const DiscreteCoordinate<DDimX> ixm(ix - 1);
            for (DiscreteCoordinate<DDimY> iy : select<DDimY>(inner_xy)) {
                const DiscreteCoordinate<DDimY> iyp(iy + 1);
                const DiscreteCoordinate<DDimY> iym(iy - 1);
                T_out(ix, iy) = T_in(ix, iy);
                T_out(ix, iy) += Cx * (T_in(ixp, iy) - 2.0 * T_in(ix, iy) + T_in(ixm, iy));
                T_out(ix, iy) += Cy * (T_in(ix, iyp) - 2.0 * T_in(ix, iy) + T_in(ix, iym));
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
