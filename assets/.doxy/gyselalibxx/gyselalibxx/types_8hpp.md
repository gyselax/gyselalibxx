

# File types.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**types.hpp**](types_8hpp.md)

[Go to the source code of this file](types_8hpp_source.md)



* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "vector_field.hpp"`
* `#include "vector_field_mem.hpp"`
* `#include "vector_index_tools.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Patch::BSplines1 | [**BSplines1OnPatch**](#typedef-bsplines1onpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: The BSplines on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef typename Patch::BSplines2 | [**BSplines2OnPatch**](#typedef-bsplines2onpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: The BSplines on the second of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef DConstField&lt; IdxRange&lt; ddc::Deriv&lt; typename Patch::Dim1 &gt;, ddc::Deriv&lt; typename Patch::Dim2 &gt; &gt; &gt; | [**ConstDeriv12\_OnPatch\_2D**](#typedef-constderiv12_onpatch_2d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant field of the derivatives_ _defined on the index ranges of valid n and m._ |
| typedef DConstField&lt; IdxRange&lt; ddc::Deriv&lt; typename Patch::Dim1 &gt;, typename Patch::Grid2 &gt; &gt; | [**ConstDeriv1\_OnPatch\_2D**](#typedef-constderiv1_onpatch_2d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant field of the n-th derivatives in the direction of_[_**Patch**_](structPatch.md) _'s first logical grid defined on the index range of n and the_[_**Patch**_](structPatch.md) _'s second logical grid._ |
| typedef DConstField&lt; IdxRange&lt; typename Patch::Grid1, ddc::Deriv&lt; typename Patch::Dim2 &gt; &gt; &gt; | [**ConstDeriv2\_OnPatch\_2D**](#typedef-constderiv2_onpatch_2d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant field of the n-th derivatives in the direction of_[_**Patch**_](structPatch.md) _'s second logical grid defined on the index range of n and the_[_**Patch**_](structPatch.md) _'s first logical grid._ |
| typedef DConstField&lt; typename Patch::IdxRangeBS1 &gt; | [**ConstSplineCoeff1OnPatch\_1D**](#typedef-constsplinecoeff1onpatch_1d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of spline coefficients for a non-batched spline defined on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef DConstField&lt; IdxRange&lt; typename Patch::BSplines1, typename Patch::Grid2 &gt; &gt; | [**ConstSplineCoeff1OnPatch\_2D**](#typedef-constsplinecoeff1onpatch_2d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of spline coefficients batched over the second of the_[_**Patch**_](structPatch.md) _'s logical dimensions for a spline defined on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef DConstField&lt; typename Patch::IdxRangeBS12 &gt; | [**ConstSplineCoeffOnPatch\_2D**](#typedef-constsplinecoeffonpatch_2d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of 2D spline coefficients for a non-batched spline defined on both of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef host\_t&lt; ConstSplineCoeffOnPatch\_2D&lt; [**Patch**](structPatch.md) &gt; &gt; | [**ConstSplineCoeffOnPatch\_2D\_host**](#typedef-constsplinecoeffonpatch_2d_host)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of 2D spline coefficients for a non-batched spline defined on both of the_[_**Patch**_](structPatch.md) _'s logical dimensions. Defined on host._ |
| typedef Field&lt; typename Patch::Coord1, typename Patch::IdxRange1 &gt; | [**Coord1Field1OnPatch\_1D**](#typedef-coord1field1onpatch_1d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: Field of 1D coordinates defined on the defined on domain of the_[_**Patch**_](structPatch.md) _'s first logical dimension._ |
| typedef ConstField&lt; typename Patch::Coord12, typename Patch::IdxRange12 &gt; | [**CoordConstFieldOnPatch**](#typedef-coordconstfieldonpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: ConstField of 2D coordinates defined on the 2D_[_**Patch**_](structPatch.md) _'s logical domain._ |
| typedef FieldMem&lt; typename Patch::Coord12, typename Patch::IdxRange12 &gt; | [**CoordFieldMemOnPatch**](#typedef-coordfieldmemonpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: Field of 2D coordinates defined on the 2D_[_**Patch**_](structPatch.md) _'s logical domain._ |
| typedef Field&lt; typename Patch::Coord12, typename Patch::IdxRange12 &gt; | [**CoordFieldOnPatch**](#typedef-coordfieldonpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: Field of 2D coordinates defined on the 2D_[_**Patch**_](structPatch.md) _'s logical domain._ |
| typedef DConstField&lt; typename Patch::IdxRange1 &gt; | [**DConstField1OnPatch**](#typedef-dconstfield1onpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant Field of doubles on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef DConstField&lt; typename Patch::IdxRange12 &gt; | [**DConstFieldOnPatch**](#typedef-dconstfieldonpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant Field defined on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._ |
| typedef DField&lt; typename Patch::IdxRange1 &gt; | [**DField1OnPatch**](#typedef-dfield1onpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A Field of doubles on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef DFieldMem&lt; typename Patch::IdxRange1 &gt; | [**DFieldMem1OnPatch**](#typedef-dfieldmem1onpatch)  <br>_A FieldMem defined on the_ [_**Patch**_](structPatch.md) _'s first logical domain._ |
| typedef DFieldMem&lt; typename Patch::IdxRange2 &gt; | [**DFieldMem2OnPatch**](#typedef-dfieldmem2onpatch)  <br>_A FieldMem defined on the_ [_**Patch**_](structPatch.md) _'s first logical domain._ |
| typedef DFieldMem&lt; typename Patch::IdxRange12 &gt; | [**DFieldMemOnPatch**](#typedef-dfieldmemonpatch)  <br>_A FieldMem defined on the_ [_**Patch**_](structPatch.md) _'s 2D logical domain._ |
| typedef DField&lt; typename Patch::IdxRange12 &gt; | [**DFieldOnPatch**](#typedef-dfieldonpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A Field defined on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._ |
| typedef host\_t&lt; DFieldOnPatch&lt; [**Patch**](structPatch.md) &gt; &gt; | [**DFieldOnPatch\_host**](#typedef-dfieldonpatch_host)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A Field defined on host on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._ |
| typedef [**VectorConstField**](classVectorField.md)&lt; double, typename Patch::IdxRange12, VectorIndexSet&lt; typename Patch::Dim1, typename Patch::Dim2 &gt; &gt; | [**DVectorConstFieldOnPatch**](#typedef-dvectorconstfieldonpatch)  <br>_A ConstVectorField defined on the_ [_**Patch**_](structPatch.md) _'s 2D logical domain._ |
| typedef [**VectorFieldMem**](classVectorFieldMem.md)&lt; double, typename Patch::IdxRange12, VectorIndexSet&lt; typename Patch::Dim1, typename Patch::Dim2 &gt; &gt; | [**DVectorFieldMemOnPatch**](#typedef-dvectorfieldmemonpatch)  <br>_A_ [_**VectorFieldMem**_](classVectorFieldMem.md) _defined on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._ |
| typedef [**VectorField**](classVectorField.md)&lt; double, typename Patch::IdxRange12, VectorIndexSet&lt; typename Patch::Dim1, typename Patch::Dim2 &gt; &gt; | [**DVectorFieldOnPatch**](#typedef-dvectorfieldonpatch)  <br>_A_ [_**VectorField**_](classVectorField.md) _defined on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._ |
| typedef typename Patch::Grid1 | [**Grid1OnPatch**](#typedef-grid1onpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: The grid on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef typename Patch::Grid2 | [**Grid2OnPatch**](#typedef-grid2onpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: The grid on the second of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef typename Patch::Idx1 | [**Idx1OnPatch**](#typedef-idx1onpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: An index for the grid on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef typename Patch::IdxRange1 | [**IdxRange1OnPatch**](#typedef-idxrange1onpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: An index range over the grids on the first_[_**Patch**_](structPatch.md) _'s 1D logical domain._ |
| typedef typename Patch::IdxRange12 | [**IdxRangeOnPatch**](#typedef-idxrangeonpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: An index range over the grids on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._ |
| typedef DField&lt; typename Patch::IdxRangeBS1 &gt; | [**SplineCoeff1OnPatch\_1D**](#typedef-splinecoeff1onpatch_1d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of spline coefficients for a non-batched spline defined on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef DField&lt; IdxRange&lt; typename Patch::BSplines1, typename Patch::Grid2 &gt; &gt; | [**SplineCoeff1OnPatch\_2D**](#typedef-splinecoeff1onpatch_2d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of spline coefficients batched over the second of the_[_**Patch**_](structPatch.md) _'s logical dimensions for a spline defined on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
| typedef DField&lt; typename Patch::IdxRangeBS12 &gt; | [**SplineCoeffOnPatch\_2D**](#typedef-splinecoeffonpatch_2d)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of 2D spline coefficients for a non-batched spline defined on both of the_[_**Patch**_](structPatch.md) _'s logical dimensions._ |
















































## Public Types Documentation




### typedef BSplines1OnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: The BSplines on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using BSplines1OnPatch =  typename Patch::BSplines1;
```




<hr>



### typedef BSplines2OnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: The BSplines on the second of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using BSplines2OnPatch =  typename Patch::BSplines2;
```




<hr>



### typedef ConstDeriv12\_OnPatch\_2D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant field of the derivatives_ _defined on the index ranges of valid n and m._
```C++
using ConstDeriv12_OnPatch_2D =  DConstField<IdxRange<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2> >>;
```




<hr>



### typedef ConstDeriv1\_OnPatch\_2D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant field of the n-th derivatives in the direction of_[_**Patch**_](structPatch.md) _'s first logical grid defined on the index range of n and the_[_**Patch**_](structPatch.md) _'s second logical grid._
```C++
using ConstDeriv1_OnPatch_2D =  DConstField<IdxRange<ddc::Deriv<typename Patch::Dim1>, typename Patch::Grid2> >;
```




<hr>



### typedef ConstDeriv2\_OnPatch\_2D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant field of the n-th derivatives in the direction of_[_**Patch**_](structPatch.md) _'s second logical grid defined on the index range of n and the_[_**Patch**_](structPatch.md) _'s first logical grid._
```C++
using ConstDeriv2_OnPatch_2D =  DConstField<IdxRange<typename Patch::Grid1, ddc::Deriv<typename Patch::Dim2> >>;
```




<hr>



### typedef ConstSplineCoeff1OnPatch\_1D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of spline coefficients for a non-batched spline defined on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using ConstSplineCoeff1OnPatch_1D =  DConstField<typename Patch::IdxRangeBS1>;
```




<hr>



### typedef ConstSplineCoeff1OnPatch\_2D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of spline coefficients batched over the second of the_[_**Patch**_](structPatch.md) _'s logical dimensions for a spline defined on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using ConstSplineCoeff1OnPatch_2D =  DConstField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2> >;
```




<hr>



### typedef ConstSplineCoeffOnPatch\_2D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of 2D spline coefficients for a non-batched spline defined on both of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using ConstSplineCoeffOnPatch_2D =  DConstField<typename Patch::IdxRangeBS12>;
```




<hr>



### typedef ConstSplineCoeffOnPatch\_2D\_host 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of 2D spline coefficients for a non-batched spline defined on both of the_[_**Patch**_](structPatch.md) _'s logical dimensions. Defined on host._
```C++
using ConstSplineCoeffOnPatch_2D_host =  host_t<ConstSplineCoeffOnPatch_2D<Patch> >;
```




<hr>



### typedef Coord1Field1OnPatch\_1D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: Field of 1D coordinates defined on the defined on domain of the_[_**Patch**_](structPatch.md) _'s first logical dimension._
```C++
using Coord1Field1OnPatch_1D =  Field<typename Patch::Coord1, typename Patch::IdxRange1>;
```




<hr>



### typedef CoordConstFieldOnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: ConstField of 2D coordinates defined on the 2D_[_**Patch**_](structPatch.md) _'s logical domain._
```C++
using CoordConstFieldOnPatch =  ConstField<typename Patch::Coord12, typename Patch::IdxRange12>;
```




<hr>



### typedef CoordFieldMemOnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: Field of 2D coordinates defined on the 2D_[_**Patch**_](structPatch.md) _'s logical domain._
```C++
using CoordFieldMemOnPatch =  FieldMem<typename Patch::Coord12, typename Patch::IdxRange12>;
```




<hr>



### typedef CoordFieldOnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: Field of 2D coordinates defined on the 2D_[_**Patch**_](structPatch.md) _'s logical domain._
```C++
using CoordFieldOnPatch =  Field<typename Patch::Coord12, typename Patch::IdxRange12>;
```




<hr>



### typedef DConstField1OnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant Field of doubles on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using DConstField1OnPatch =  DConstField<typename Patch::IdxRange1>;
```




<hr>



### typedef DConstFieldOnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A constant Field defined on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._
```C++
using DConstFieldOnPatch =  DConstField<typename Patch::IdxRange12>;
```




<hr>



### typedef DField1OnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A Field of doubles on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using DField1OnPatch =  DField<typename Patch::IdxRange1>;
```




<hr>



### typedef DFieldMem1OnPatch 

_A FieldMem defined on the_ [_**Patch**_](structPatch.md) _'s first logical domain._
```C++
using DFieldMem1OnPatch =  DFieldMem<typename Patch::IdxRange1>;
```




<hr>



### typedef DFieldMem2OnPatch 

_A FieldMem defined on the_ [_**Patch**_](structPatch.md) _'s first logical domain._
```C++
using DFieldMem2OnPatch =  DFieldMem<typename Patch::IdxRange2>;
```




<hr>



### typedef DFieldMemOnPatch 

_A FieldMem defined on the_ [_**Patch**_](structPatch.md) _'s 2D logical domain._
```C++
using DFieldMemOnPatch =  DFieldMem<typename Patch::IdxRange12>;
```




<hr>



### typedef DFieldOnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A Field defined on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._
```C++
using DFieldOnPatch =  DField<typename Patch::IdxRange12>;
```




<hr>



### typedef DFieldOnPatch\_host 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A Field defined on host on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._
```C++
using DFieldOnPatch_host =  host_t<DFieldOnPatch<Patch> >;
```




<hr>



### typedef DVectorConstFieldOnPatch 

_A ConstVectorField defined on the_ [_**Patch**_](structPatch.md) _'s 2D logical domain._
```C++
using DVectorConstFieldOnPatch =  VectorConstField< double, typename Patch::IdxRange12, VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2> >;
```




<hr>



### typedef DVectorFieldMemOnPatch 

_A_ [_**VectorFieldMem**_](classVectorFieldMem.md) _defined on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._
```C++
using DVectorFieldMemOnPatch =  VectorFieldMem< double, typename Patch::IdxRange12, VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2> >;
```




<hr>



### typedef DVectorFieldOnPatch 

_A_ [_**VectorField**_](classVectorField.md) _defined on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._
```C++
using DVectorFieldOnPatch =  VectorField< double, typename Patch::IdxRange12, VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2> >;
```




<hr>



### typedef Grid1OnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: The grid on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using Grid1OnPatch =  typename Patch::Grid1;
```




<hr>



### typedef Grid2OnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: The grid on the second of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using Grid2OnPatch =  typename Patch::Grid2;
```




<hr>



### typedef Idx1OnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: An index for the grid on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using Idx1OnPatch =  typename Patch::Idx1;
```




<hr>



### typedef IdxRange1OnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: An index range over the grids on the first_[_**Patch**_](structPatch.md) _'s 1D logical domain._
```C++
using IdxRange1OnPatch =  typename Patch::IdxRange1;
```




<hr>



### typedef IdxRangeOnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: An index range over the grids on the_[_**Patch**_](structPatch.md) _'s 2D logical domain._
```C++
using IdxRangeOnPatch =  typename Patch::IdxRange12;
```




<hr>



### typedef SplineCoeff1OnPatch\_1D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of spline coefficients for a non-batched spline defined on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using SplineCoeff1OnPatch_1D =  DField<typename Patch::IdxRangeBS1>;
```




<hr>



### typedef SplineCoeff1OnPatch\_2D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of spline coefficients batched over the second of the_[_**Patch**_](structPatch.md) _'s logical dimensions for a spline defined on the first of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using SplineCoeff1OnPatch_2D =  DField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2> >;
```




<hr>



### typedef SplineCoeffOnPatch\_2D 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of 2D spline coefficients for a non-batched spline defined on both of the_[_**Patch**_](structPatch.md) _'s logical dimensions._
```C++
using SplineCoeffOnPatch_2D =  DField<typename Patch::IdxRangeBS12>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/data_types/types.hpp`

