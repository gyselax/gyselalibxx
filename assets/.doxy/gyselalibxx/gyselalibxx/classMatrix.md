

# Class Matrix



[**ClassList**](annotated.md) **>** [**Matrix**](classMatrix.md)



_The super class from which matrix classes should inherit. This class is used to solve matrix equations._ 

* `#include <matrix.hpp>`





Inherited by the following classes: [Matrix\_Banded](classMatrix__Banded.md),  [Matrix\_Corner\_Block](classMatrix__Corner__Block.md),  [Matrix\_Dense](classMatrix__Dense.md),  [Matrix\_PDS\_Tridiag](classMatrix__PDS__Tridiag.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Matrix**](#function-matrix) (int mat\_size) <br>_A constructor for the matrix._  |
| virtual void | [**factorise**](#function-factorise) () <br>_Factorise the matrix. This method prepares the matrix for a call to the solve methods. For most matrix types a call to factorise causes the LU decomposition to be calculated. The elements of the matrix should not be accessed once this method has been called._  |
| virtual double | [**get\_element**](#function-get_element) (int i, int j) const = 0<br>_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._ |
|  int | [**get\_size**](#function-get_size) () const<br>_Get the size of the n x n_ [_**Matrix**_](classMatrix.md) _._ |
| virtual void | [**set\_element**](#function-set_element) (int i, int j, double a\_ij) = 0<br>_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual DSpan1D | [**solve\_inplace**](#function-solve_inplace) (DSpan1D b) const<br>_Solve the matrix equation in place._  |
| virtual DSpan2D | [**solve\_multiple\_inplace**](#function-solve_multiple_inplace) (DSpan2D bx) const<br>_Solve multiple matrix equations in place._  |
| virtual DSpan1D | [**solve\_transpose\_inplace**](#function-solve_transpose_inplace) (DSpan1D b) const<br>_Solve the transposed matrix equation in place._  |
| virtual  | [**~Matrix**](#function-matrix) () = default<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_banded**](#function-make_new_banded) (int n, int kl, int ku, bool pds) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a banded matrix. This method returns the appropriate subclass to minimise memory and computations._ |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_block\_with\_banded\_region**](#function-make_new_block_with_banded_region) (int n, int kl, int ku, bool pds, int block1\_size, int block2\_size=0) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a matrix which has a banded region. This method returns the appropriate subclass to minimise memory and computations. This matrix must be able to be described by the following block matrices:_ |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_periodic\_banded**](#function-make_new_periodic_banded) (int n, int kl, int ku, bool pds) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a periodic banded matrix. This method returns the appropriate subclass to minimise memory and computations. A periodic banded matrix is like a banded matrix but additionally contains non- zero values in the corners. I.e it has the following sparsity pattern:_ |






## Protected Attributes

| Type | Name |
| ---: | :--- |
|  int const | [**n**](#variable-n)  <br>_The matrix size._  |
















## Protected Functions

| Type | Name |
| ---: | :--- |
| virtual int | [**factorise\_method**](#function-factorise_method) () = 0<br>_Call the factorisation method._  |
| virtual int | [**solve\_inplace\_method**](#function-solve_inplace_method) (double \* b, char transpose, int n\_equations) const = 0<br>_Call the LAPACK solve method._  |




## Public Functions Documentation




### function Matrix 

_A constructor for the matrix._ 
```C++
inline Matrix::Matrix (
    int mat_size
) 
```





**Parameters:**


* `mat_size` The size of the n x n [**Matrix**](classMatrix.md). 




        

<hr>



### function factorise 

_Factorise the matrix. This method prepares the matrix for a call to the solve methods. For most matrix types a call to factorise causes the LU decomposition to be calculated. The elements of the matrix should not be accessed once this method has been called._ 
```C++
virtual void Matrix::factorise () 
```




<hr>



### function get\_element 

_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual double Matrix::get_element (
    int i,
    int j
) const = 0
```





**Parameters:**


* `i` The row index. 
* `j` The column index. 



**Returns:**

The element at position (i,j). 





        

<hr>



### function get\_size 

_Get the size of the n x n_ [_**Matrix**_](classMatrix.md) _._
```C++
inline int Matrix::get_size () const
```





**Returns:**

The size of the n x n [**Matrix**](classMatrix.md). 





        

<hr>



### function set\_element 

_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual void Matrix::set_element (
    int i,
    int j,
    double a_ij
) = 0
```





**Parameters:**


* `i` The row index. 
* `j` The column index. 
* `a_ij` The value that the element should be set to. 




        

<hr>



### function solve\_inplace 

_Solve the matrix equation in place._ 
```C++
virtual DSpan1D Matrix::solve_inplace (
    DSpan1D b
) const
```



Solve the following matrix equation:  The result  is saved into the memory allocated for .




**Parameters:**


* `b` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        

<hr>



### function solve\_multiple\_inplace 

_Solve multiple matrix equations in place._ 
```C++
virtual DSpan2D Matrix::solve_multiple_inplace (
    DSpan2D bx
) const
```



Solve the following matrix equation:  for multiple values of  and . The first dimension is iterated over with each slice representing an equation to be solved. The result  is saved into the memory allocated for .




**Parameters:**


* `bx` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        

<hr>



### function solve\_transpose\_inplace 

_Solve the transposed matrix equation in place._ 
```C++
virtual DSpan1D Matrix::solve_transpose_inplace (
    DSpan1D b
) const
```



Solve the following matrix equation:  The result  is saved into the memory allocated for .




**Parameters:**


* `b` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        

<hr>



### function ~Matrix 

```C++
virtual Matrix::~Matrix () = default
```




<hr>
## Public Static Functions Documentation




### function make\_new\_banded 

_Get a new_ [_**Matrix**_](classMatrix.md) _representing a banded matrix. This method returns the appropriate subclass to minimise memory and computations._
```C++
static std::unique_ptr< Matrix > Matrix::make_new_banded (
    int n,
    int kl,
    int ku,
    bool pds
) 
```





**Parameters:**


* `n` The size of the n x n [**Matrix**](classMatrix.md). 
* `kl` The number of lower diagonals. 
* `ku` The number of upper diagonals. 
* `pds` True if the matrix is positive-definite and symmetric.



**Returns:**

An instance of a [**Matrix**](classMatrix.md) class that can store the banded matrix. 





        

<hr>



### function make\_new\_block\_with\_banded\_region 

_Get a new_ [_**Matrix**_](classMatrix.md) _representing a matrix which has a banded region. This method returns the appropriate subclass to minimise memory and computations. This matrix must be able to be described by the following block matrices:_
```C++
static std::unique_ptr< Matrix > Matrix::make_new_block_with_banded_region (
    int n,
    int kl,
    int ku,
    bool pds,
    int block1_size,
    int block2_size=0
) 
```



 where E is a banded matrix, and A, E and I are square matrices.




**Parameters:**


* `n` The size of the n x n [**Matrix**](classMatrix.md). 
* `kl` The number of lower diagonals. 
* `ku` The number of upper diagonals. 
* `pds` True if the matrix is positive-definite and symmetric. 
* `block1_size` The size of the matrix A. 
* `block2_size` The size of the matrix I.



**Returns:**

An instance of a [**Matrix**](classMatrix.md) class that can store the banded matrix. 





        

<hr>



### function make\_new\_periodic\_banded 

_Get a new_ [_**Matrix**_](classMatrix.md) _representing a periodic banded matrix. This method returns the appropriate subclass to minimise memory and computations. A periodic banded matrix is like a banded matrix but additionally contains non- zero values in the corners. I.e it has the following sparsity pattern:_
```C++
static std::unique_ptr< Matrix > Matrix::make_new_periodic_banded (
    int n,
    int kl,
    int ku,
    bool pds
) 
```



 

**Parameters:**


* `n` The size of the n x n [**Matrix**](classMatrix.md). 
* `kl` The number of lower diagonals. 
* `ku` The number of upper diagonals. 
* `pds` True if the matrix is positive-definite and symmetric.



**Returns:**

An instance of a [**Matrix**](classMatrix.md) class that can store the banded matrix. 





        

<hr>
## Protected Attributes Documentation




### variable n 

_The matrix size._ 
```C++
int const Matrix::n;
```




<hr>
## Protected Functions Documentation




### function factorise\_method 

_Call the factorisation method._ 
```C++
virtual int Matrix::factorise_method () = 0
```





**Returns:**

The LAPACK error code. 





        

<hr>



### function solve\_inplace\_method 

_Call the LAPACK solve method._ 
```C++
virtual int Matrix::solve_inplace_method (
    double * b,
    char transpose,
    int n_equations
) const = 0
```





**Parameters:**


* `b` The data describing the right-hand side. 
* `transpose` A character ['N'/'[**T**](structT.md)'] describing whether the matrix or the transposed matrix appears in the matrix equation. 
* `n_equations` The number of equations being solved. 



**Returns:**

The LAPACK error code. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix.hpp`

