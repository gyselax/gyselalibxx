# SPDX-License-Identifier: MIT
""" File which allows the visualisation of basis splines
"""
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt

parser = ArgumentParser(description="View the shape of the basis-splines")
parser.add_argument('--degree', type = int, nargs = '?', default = 2, help = "The degree of the b-splines")
parser.add_argument('--ncells', type = int, nargs = '?', default = 5, help = "The number of cells in the domain")
parser.add_argument('--greville_bcs', action='store_true', help="Place all knots outside the domain at the boundary")
args = parser.parse_args()

nc = args.ncells # Number of cells
d = args.degree  # Degree of the bsplines
greville_boundary_knots = args.greville_bcs

if greville_boundary_knots:
    knots = np.array([0]*d+[*np.linspace(0,nc,nc+1)] + [nc]*d)
else:
    knots = np.linspace(-d,nc+d,nc+1+2*d)

def linear_bspline_eval(t,i):
    """
    Evaluate the i-th bspline of degree 1 at position t

    Parameters
    ----------
    t : float
        The position at which the bspline is evaluated
    i : int
        The index of the bspline
    """
    if knots[i] < t <= knots[i+1]:
        return 1
    else:
        return 0

linear_bspline_eval_vect = np.vectorize(linear_bspline_eval)

def bspline_eval(t,i,m):
    """
    Evaluate the i-th bspline of degree m at position t

    Parameters
    ----------
    t : float/array of floats
        The position(s) at which the bspline is evaluated
    i : int
        The index of the bspline
    m : int
        The degree of the bsplines
    """
    if m==0:
        return linear_bspline_eval_vect(t,i)
    else:
        part1 = (t-knots[i])/(knots[i+m]-knots[i])*bspline_eval(t,i,m-1) if knots[i]!=knots[i+m] else 0
        part2 = (knots[i+m+1]-t)/(knots[i+m+1]-knots[i+1]) * bspline_eval(t,i+1,m-1) if knots[i+m+1] != knots[i+1] else 0
        return part1 + part2

if greville_boundary_knots:
    x = np.linspace(0,nc,nc*100)
else:
    x = np.linspace(-d,nc+d,(nc+2*d)*20)
greville_pts = np.array( [np.sum(knots[i:i+d])/d for i in range(1,1+nc+d)] )
greville_pts = np.around(greville_pts, decimals = 15)

plt.figure(figsize=[12,6])
plt.rc('text', usetex=True)
plt.rc('font', size=20)
ax = plt.gca()

# Plot the knots
ax.axvline(knots[0], color='k', label='knots')
if greville_boundary_knots:
    for k in knots[d+1:nc+d]:
        ax.axvline(k, color='k')
else:
    for k in knots[1:]:
        ax.axvline(k, color='k')

# Plot the domain boundaries
ax.axvline(0, color='k', linewidth=4)
ax.axvline(nc, color='k', linewidth=4)

# Plot the bsplines
for i in range(nc+d):
    plt.plot(x,bspline_eval(x,i,d),label="$B_"+str(i)+"(x)$",color='C'+str(i),linewidth=4)

# Plot the greville points
for g in greville_pts:
    colour = 'k' if 0 <= g <= nc else 'r'
    plt.plot(g, 0, 'o', color = colour, markersize=10)
# Add greville labels
plt.plot(0, 2, 'ok', label='greville\n points', markersize=10)
plt.plot(0, 2, 'or', label='points which\n become\n boundary\n conditions', markersize=10)

plt.xticks(knots)
plt.title('B-Splines')
plt.ylim(-0.05,1.25)

ax = plt.gca()
# Shrink current axis by 10%
box = ax.get_position()
ax.set_position([box.x0/2, box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()
