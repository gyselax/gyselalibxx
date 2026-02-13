import numpy as np
import xarray as xr

# -------------------------------------------------
# 1. Quadrature Coefficients
# -------------------------------------------------
def trapezoid_quadrature_coefficients_1d(grid):
    """
    Computes trapezoidal weights for a 1D non-uniform grid.
    Formula: w_i = 0.5 * (x_{i+1} - x_{i-1})
    """
    grid = np.asarray(grid)
    n = len(grid)
    coeffs = np.zeros(n)
    if n > 1:
        coeffs[0] = 0.5 * (grid[1] - grid[0])
        coeffs[-1] = 0.5 * (grid[-1] - grid[-2])
        coeffs[1:-1] = 0.5 * (grid[2:] - grid[:-2])
    return coeffs

def get_quadrature_weights(vpar_vals, mu_vals):
    """
    Generates the 2D weight matrix for (vpar, mu) integration.
    """
    w_vpar = trapezoid_quadrature_coefficients_1d(vpar_vals)
    w_mu = trapezoid_quadrature_coefficients_1d(mu_vals)
    weights = np.outer(w_vpar, w_mu)
    return weights

# -------------------------------------------------
# 2. FluidMoments Class
# -------------------------------------------------
class FluidMoments:
    def __init__(self, vpar_coord, mu_coord):
        self.vpar = xr.DataArray(vpar_coord, dims="vpar", coords={"vpar": vpar_coord})
        self.mu = xr.DataArray(mu_coord, dims="mu", coords={"mu": mu_coord})
        
        w_vpar = trapezoid_quadrature_coefficients_1d(vpar_coord)
        w_mu = trapezoid_quadrature_coefficients_1d(mu_coord)
        weights_2d = np.outer(w_vpar, w_mu)
        
        self.weights = xr.DataArray(
            weights_2d, 
            dims=("vpar", "mu"), 
            coords={"vpar": vpar_coord, "mu": mu_coord}
        )

    def compute_density(self, f):
        return (f * self.weights).sum(dim=("vpar", "mu"))

    def compute_velocity(self, f, density):
        momentum = (f * self.vpar * self.weights).sum(dim=("vpar", "mu"))
        return momentum / density

    def compute_temperature(self, f, density, velocity):
        diff2 = (self.vpar - velocity)**2
        return (f * diff2 * self.weights).sum(dim=("vpar", "mu")) / density
