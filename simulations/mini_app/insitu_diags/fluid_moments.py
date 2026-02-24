import xarray as xr

# -------------------------------------------------
# 1. Quadrature
# -------------------------------------------------
def get_trapezoid_weights_1d(coord):
    """
    Computes trapezoidal weights using Xarray vectorization.
    Formula: w_i = 0.5 * (x_{i+1} - x_{i-1})
    """
    forward = coord.shift({coord.dims[0]: -1})
    backward = coord.shift({coord.dims[0]: 1})

    weights = 0.5 * (forward - backward)

    first_dx = 0.5 * (coord.isel({coord.dims[0]: 1}) - coord.isel({coord.dims[0]: 0}))
    last_dx = 0.5 * (coord.isel({coord.dims[0]: -1}) - coord.isel({coord.dims[0]: -2}))

    return weights.fillna(0.0).where(coord != coord[0], first_dx).where(coord != coord[-1], last_dx)

# -------------------------------------------------
# 2. FluidMoments Class
# -------------------------------------------------
class FluidMoments:
    def __init__(self, vpar_coord, mu_coord):
        self.vpar = vpar_coord
        self.mu = mu_coord
        
        w_vpar = get_trapezoid_weights_1d(vpar_coord)
        w_mu = get_trapezoid_weights_1d(mu_coord)
        
        self.weights = w_vpar * w_mu

    def compute_density(self, f):
        return (f * self.weights).sum(dim=("vpar", "mu"))

    def compute_velocity(self, f, density):
        momentum = (f * self.vpar * self.weights).sum(dim=("vpar", "mu"))
        return momentum / density

    def compute_temperature(self, f, density, velocity):
        diff2 = (self.vpar - velocity)**2
        return (f * diff2 * self.weights).sum(dim=("vpar", "mu")) / density
