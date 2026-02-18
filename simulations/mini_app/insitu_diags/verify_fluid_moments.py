import xarray as xr
from fluid_moments import FluidMoments

# -------------------------------------------------
# 3. Execution Script (Dask-Backed)
# -------------------------------------------------
f_path = "fdistribu_5D_output.h5"
mom_path = "fluid_moments.h5"

ds_f = xr.open_dataset(f_path, chunks={"phony_dim_1": 1, "phony_dim_2": -1})
ds_ref = xr.open_dataset(mom_path, chunks="auto")

rename_dict = {
    "phony_dim_0": "species", "phony_dim_1": "tor1", "phony_dim_2": "tor2",
    "phony_dim_3": "tor3", "phony_dim_4": "vpar", "phony_dim_5": "mu"
}
ds_f = ds_f.rename({k: v for k, v in rename_dict.items() if k in ds_f.dims})

ref_rename = {k: v for k, v in rename_dict.items() if k in ds_ref.dims}
ds_ref = ds_ref.rename(ref_rename)

moments_calc = FluidMoments(ds_f["vpar"], ds_f["mu"])

fdistribu = ds_f["fdistribu_sptor3Dv2D"]

# --- Computation ---
density = moments_calc.compute_density(fdistribu)
mean_velocity = moments_calc.compute_velocity(fdistribu, density)
temperature = moments_calc.compute_temperature(fdistribu, density, mean_velocity)

# --- Verification & Output ---
ds_out = xr.Dataset({
    "density": density,
    "mean_velocity": mean_velocity,
    "temperature": temperature,
}).compute()

diff = abs(ds_out - ds_ref).to_array()
print(f"Max Abs Error: {diff.max().values}")

ds_out.to_netcdf("fluid_moments.nc")
