

# Class Matrix\_Dense



[**ClassList**](annotated.md) **>** [**Matrix\_Dense**](classMatrix__Dense.md)



_A class describing a dense matrix._ 

* `#include <matrix_dense.hpp>`



Inherits the following classes: [Matrix](classMatrix.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Matrix\_Dense**](#function-matrix_dense) (int n) <br>_A constructor for the matrix._  |
| virtual double | [**get\_element**](#function-get_element) (int i, int j) override const<br>_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual void | [**set\_element**](#function-set_element) (int i, int j, double a\_ij) override<br>_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._ |


## Public Functions inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
|   | [**Matrix**](classMatrix.md#function-matrix) (int mat\_size) <br>_A constructor for the matrix._  |
| virtual void | [**factorise**](classMatrix.md#function-factorise) () <br>_Factorise the matrix. This method prepares the matrix for a call to the solve methods. For most matrix types a call to factorise causes the LU decomposition to be calculated. The elements of the matrix should not be accessed once this method has been called._  |
| virtual double | [**get\_element**](classMatrix.md#function-get_element) (int i, int j) const = 0<br>_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._ |
|  int | [**get\_size**](classMatrix.md#function-get_size) () const<br>_Get the size of the n x n_ [_**Matrix**_](classMatrix.md) _._ |
| virtual void | [**set\_element**](classMatrix.md#function-set_element) (int i, int j, double a\_ij) = 0<br>_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual DSpan1D | [**solve\_inplace**](classMatrix.md#function-solve_inplace) (DSpan1D b) const<br>_Solve the matrix equation in place._  |
| virtual DSpan2D | [**solve\_multiple\_inplace**](classMatrix.md#function-solve_multiple_inplace) (DSpan2D bx) const<br>_Solve multiple matrix equations in place._  |
| virtual DSpan1D | [**solve\_transpose\_inplace**](classMatrix.md#function-solve_transpose_inplace) (DSpan1D b) const<br>_Solve the transposed matrix equation in place._  |
| virtual  | [**~Matrix**](classMatrix.md#function-matrix) () = default<br> |




## Public Static Functions inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_banded**](classMatrix.md#function-make_new_banded) (int n, int kl, int ku, bool pds) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a banded matrix. This method returns the appropriate subclass to minimise memory and computations._ |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_block\_with\_banded\_region**](classMatrix.md#function-make_new_block_with_banded_region) (int n, int kl, int ku, bool pds, int block1\_size, int block2\_size=0) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a matrix which has a banded region. This method returns the appropriate subclass to minimise memory and computations. This matrix must be able to be described by the following block matrices:_ |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_periodic\_banded**](classMatrix.md#function-make_new_periodic_banded) (int n, int kl, int ku, bool pds) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a periodic banded matrix. This method returns the appropriate subclass to minimise memory and computations. A periodic banded matrix is like a banded matrix but additionally contains non- zero values in the corners. I.e it has the following sparsity pattern:_ |












## Protected Attributes inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
|  int const | [**n**](classMatrix.md#variable-n)  <br>_The matrix size._  |
































## Protected Functions inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
| virtual int | [**factorise\_method**](classMatrix.md#function-factorise_method) () = 0<br>_Call the factorisation method._  |
| virtual int | [**solve\_inplace\_method**](classMatrix.md#function-solve_inplace_method) (double \* b, char transpose, int n\_equations) const = 0<br>_Call the LAPACK solve method._  |






## Public Functions Documentation




### function Matrix\_Dense 

_A constructor for the matrix._ 
```C++
explicit Matrix_Dense::Matrix_Dense (
    int n
) 
```





**Parameters:**


* `n` The size of the n x n [**Matrix**](classMatrix.md). 




        

<hr>



### function get\_element 

_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual double Matrix_Dense::get_element (
    int i,
    int j
) override const
```





**Parameters:**


* `i` The row index. 
* `j` The column index. 



**Returns:**

The element at position (i,j). 





        
Implements [*Matrix::get\_element*](classMatrix.md#function-get_element)


<hr>



### function set\_element 

_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual void Matrix_Dense::set_element (
    int i,
    int j,
    double a_ij
) override
```





**Parameters:**


* `i` The row index. 
* `j` The column index. 
* `a_ij` The value that the element should be set to. 




        
Implements [*Matrix::set\_element*](classMatrix.md#function-set_element)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_dense.hpp`

