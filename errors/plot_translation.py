"""Plot foot positions and values from translation.out on polar axes."""

import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Parse: (r, theta) (r_foot, theta_foot) val_xy val_rtheta exact
pattern = re.compile(
    r'\(\s*([^,]+),\s*([^)]+)\)\s+\(\s*([^,]+),\s*([^)]+)\)\s+(\S+)\s+(\S+)\s+(\S+)'
)

n_iter = 9
it = 1

r_coord, theta_coord = [], []
r_foot, theta_foot = [], []
val_xy, val_rtheta, exact = [], [], []

filename = sys.argv[1]
name = filename.removesuffix('.out')

with open(filename) as f:
    for line in f:
        m = pattern.match(line.strip())
        if m:
            r_coord.append(float(m.group(1)))
            theta_coord.append(float(m.group(2)))
            r_foot.append(float(m.group(3)))
            theta_foot.append(float(m.group(4)))
            val_rtheta.append(float(m.group(5)))
            val_xy.append(float(m.group(6)))
            exact.append(float(m.group(7)))

r_coord = np.array(r_coord)
theta_coord = np.array(theta_coord)
r_foot = np.array(r_foot)
theta_foot = np.array(theta_foot)
val_xy = np.array(val_xy)
val_rtheta = np.array(val_rtheta)
exact = np.array(exact)

# Infer grid from unique coordinate values (data ordered: r outer, theta inner)
r_vals = np.unique(r_coord)
theta_vals = np.unique(theta_coord)
n_r, n_theta = len(r_vals), len(theta_vals)
print(f"Grid: {n_r} r points x {n_theta} theta points ({len(r_coord)} total)")

def to_grid(arr):
    return arr.reshape(n_iter, n_r, n_theta)

R = to_grid(r_coord)
T = to_grid(theta_coord)
R_foot = to_grid(r_foot)
T_foot = to_grid(theta_foot)
VAL_XY = to_grid(val_xy)
VAL_RT = to_grid(val_rtheta)
EXACT = to_grid(exact)

fig = plt.figure(figsize=(18, 6))
fig.suptitle(name, fontsize=14)
ax_quiver = fig.add_subplot(131)
ax_xy = fig.add_subplot(131, projection='polar')
ax_rt = fig.add_subplot(132, projection='polar')
ax_ex = fig.add_subplot(133, projection='polar')

## --- Plot 1: quiver from grid point to foot in (theta, r) space ---
## Arrows show displacement directly in polar coordinates
#dtheta = T_foot[it] - T[it]
#dr = R_foot[it] - R[it]
#ax_quiver.quiver(
#    T[it].ravel(), R[it].ravel(),
#    dtheta.ravel(), dr.ravel(),
#    angles='xy', scale_units='xy', scale=1,
#    width=0.002, headwidth=4, headlength=5,
#)
#ax_quiver.set_xlabel('θ (rad)')
#ax_quiver.set_ylabel('r')
#ax_quiver.set_title('Feet (θ, r displacement)', pad=10)

# Extend theta/r by one cell to close the polar ring in pcolormesh
dtheta_grid = theta_vals[1] - theta_vals[0]
theta_edges = np.append(theta_vals, theta_vals[-1] + dtheta_grid)
dr_grid = r_vals[1] - r_vals[0]
r_edges = np.append(r_vals[0] - dr_grid / 2, r_vals + dr_grid / 2)

def polar_pcolormesh(ax, data, cmap, label, lognorm=False):
    vmin, vmax = data.min(), data.max()
    if lognorm and vmin <= 0:
        vmin = data[data > 0].min() if (data > 0).any() else 1e-16
    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax) if lognorm else None
    pc = ax.pcolormesh(theta_edges, r_edges, data, cmap=cmap, norm=norm, shading='flat')
    plt.colorbar(pc, ax=ax, label=label, shrink=0.75, pad=0.1)

# --- Plot 2: val_xy ---
polar_pcolormesh(ax_xy, VAL_XY[it], 'plasma', 'val_xy')
ax_xy.set_title('val_xy', pad=15)

# --- Plot 3: val_rtheta ---
polar_pcolormesh(ax_rt, VAL_RT[it], 'plasma', 'val_rθ')
ax_rt.set_title('val_rθ', pad=15)

# --- Plot 3: val_rtheta ---
polar_pcolormesh(ax_ex, EXACT[it], 'plasma', 'exact')
ax_ex.set_title('val_exact', pad=15)

print("Max error in x,y : ", np.max(np.abs(VAL_XY[it] - EXACT[it])))
print("Max error in r,theta : ", np.max(np.abs(VAL_RT[it] - EXACT[it])))
print("Max diff : ", np.max(np.abs(VAL_RT[it] - VAL_XY[it])))

plt.tight_layout()
plt.savefig(f'{name}_plots.png', dpi=150, bbox_inches='tight')
print(f"Saved {name}_plots.png")
plt.show()
