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

# -------------------------------------------------
# 3. Fetch Data and Compute Fluid Moments
# -------------------------------------------------
f_path = "fdistribu_5D_output.h5"
mom_path = "fluid_moments.h5"

ds_f = xr.open_dataset(f_path, chunks="auto")
ds_ref = xr.open_dataset(mom_path)

ds_f = ds_f.rename({
    "phony_dim_0": "species",
    "phony_dim_1": "tor1",
    "phony_dim_2": "tor2",
    "phony_dim_3": "tor3",
    "phony_dim_4": "vpar",
    "phony_dim_5": "mu",
})

ds_ref = ds_ref.rename({
    "phony_dim_0": "species",
    "phony_dim_1": "tor1",
    "phony_dim_2": "tor2",
    "phony_dim_3": "tor3",
})

vpar_vals = ds_f["vpar"].values 
mu_vals = ds_f["mu"].values

f = ds_f["fdistribu_sptor3Dv2D"].assign_coords(vpar=vpar_vals, mu=mu_vals)

moments_calc = FluidMoments(vpar_vals, mu_vals)

# --- Compute Density ---
density = moments_calc.compute_density(f).compute()
# --- Compute Velocity ---
mean_velocity = moments_calc.compute_velocity(f, density).compute()
# --- Compute Temperature ---
temperature = moments_calc.compute_temperature(f, density, mean_velocity).compute()

# -------------------------------------------------
# 4. Verification
# -------------------------------------------------
density_ref = ds_ref["density"]
velocity_ref = ds_ref["mean_velocity"]
temperature_ref = ds_ref["temperature"]


def compare(name, computed, reference):
    diff = np.abs(computed - reference)
    max_abs = diff.max()
    max_rel = max_abs / np.max(np.abs(reference))
    print(f"\n{name}")
    print(f"  max abs error: {max_abs.values:.8e}")
    print(f"  max rel error: {max_rel.values:.8e}")

compare("Density", density, density_ref)
compare("Mean velocity", mean_velocity, velocity_ref)
compare("Temperature", temperature, temperature_ref)
ds_out = xr.Dataset({
'density': density,
'mean_velocity': mean_velocity,
'temperature': temperature,
})
ds_out.to_netcdf('fluid_moments.nc')

