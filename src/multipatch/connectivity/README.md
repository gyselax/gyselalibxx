# Multipatch connectivity 

## Patch tag
The tag Patch refers to a single 2D patch geometry. The tag contains aliases (or shortcuts) to the DDC geometry elements, such as: 

* Grids (or points sequence along one dimension) on both dimensions (`Grid1` and `Grid2`);  
* Associated continuous dimensions (`Dim1` and `Dim2`);  
* Dimensions for B-splines coefficients (`Bsplines1` and `Bsplines2`);  
* Coordinates of objects represented on the dimensions (`Coord1`, `Coord2` and `Coord12`);  
* The type which describes the index of a grid point (e.g. `Idx1`);  
* The type which describes a distance between grid points (e.g. `IdxStep1`); 
* The type which describes the index range of a grid (e.g. `IdxRange1`).
* The type which describes the index range of a spline coefficients grid (e.g. `BSIdxRange1`).

The domain defined on this patch is called *logical domain*. 