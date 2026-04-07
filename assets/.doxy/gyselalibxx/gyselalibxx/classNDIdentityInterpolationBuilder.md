

# Class NDIdentityInterpolationBuilder

**template &lt;class ExecSpace, class MemorySpace, class DataType, class IdxRangeInterpolation, class IdxRangeBasis&gt;**



[**ClassList**](annotated.md) **>** [**NDIdentityInterpolationBuilder**](classNDIdentityInterpolationBuilder.md)



_An ND builder that copies function values directly as interpolation coefficients._ [More...](#detailed-description)


































































## Detailed Description


An ND extension of [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md). No computation is required because the interpolation coefficients equal the function values at the grid nodes (e.g. for an ND Lagrange interpolation on a tensor-product grid of nodes).




**Template parameters:**


* `ExecSpace` The Kokkos execution space. 
* `MemorySpace` The Kokkos memory space. 
* `DataType` The data type of field values and coefficients. 
* `IdxRangeInterpolation` The ND index range for the interpolation mesh, of the form IdxRange&lt;Grid1, Grid2, ...&gt;. 
* `IdxRangeBasis` The ND index range for the basis types, one per interpolation dimension, in the same order as the grids in IdxRangeInterpolation. 




    

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/nd_identity_interpolation_builder.hpp`

