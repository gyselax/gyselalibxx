# Multipatch geometries

This folder contains a variety of geometries that can be used to test different parts of the multipatch code. These geometries are created for testing purposes only so some may seem unnecessarily complex or describe domains that cannot be visualised in 2D.

The available geometries are:
- `2patches_2d_non_periodic_non_uniform` : A rectangular domain arbitrarily split into 2 patches with non-uniform grids.
- `2patches_2d_non_periodic_uniform` : A rectangular domain arbitrarily split into 2 patches with uniform grids.
- `2patches_2d_onion_shape_uniform` : A circular domain split into 2 patches along the radial direction with uniform grids.
- `3patches_2d_non_periodic_non_uniform` : A rectangular domain arbitrarily split into 3 patches with non-uniform grids.
- `3patches_2d_non_periodic_uniform` : A rectangular domain arbitrarily split into 3 patches with uniform grids.
- `5patches_figure_of_eight` : A complex domain describing the surface of a 3D object that can be used to test unusual interfaces (e.g. patches connected by more than one interface, continuous grid lines that cross the same domain twice). The domain is split into 5 patches which each have uniform grids.
- `9patches_2d_periodic_strips_uniform` : A rectangular domain split into 3x3 rectangular patches with uniform grids.
