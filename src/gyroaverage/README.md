# Gyroaverage operator (for circular geometry)

Calculate the gyroaveraging operator by averaging over sample points on circular trajectories.
The sampling points on the circular trajectories are constructed using Cartesian coordinates which
are then transformed into polar coordinates. The degree 3 spline interpolations are used on
polar coordinates to evaluate values on the sampling points. To handle a more general geometry,
this operator needs to be updated.

