

# Class Matrix\_Centre\_Block



[**ClassList**](annotated.md) **>** [**Matrix\_Centre\_Block**](classMatrix__Centre__Block.md)



_A_ [_**Matrix**_](classMatrix.md) _representing a matrix which has a banded region. This matrix must be able to be described by the following block matrices:_[More...](#detailed-description)

* `#include <matrix_centre_block.hpp>`



Inherits the following classes: [Matrix\_Corner\_Block](classMatrix__Corner__Block.md)










































































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Matrix\_Centre\_Block**](#function-matrix_centre_block) (int n, int top\_block\_size, int bottom\_block\_size, std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; q) <br>_A constructor for the matrix._  |
| virtual double | [**get\_element**](#function-get_element) (int i, int j) override const<br>_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual void | [**set\_element**](#function-set_element) (int i, int j, double a\_ij) override<br>_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual DSpan1D | [**solve\_inplace**](#function-solve_inplace) (DSpan1D b) override const<br>_Solve the matrix equation in place._  |
| virtual DSpan2D | [**solve\_multiple\_inplace**](#function-solve_multiple_inplace) (DSpan2D bx) override const<br>_Solve multiple matrix equations in place._  |
| virtual DSpan1D | [**solve\_transpose\_inplace**](#function-solve_transpose_inplace) (DSpan1D b) override const<br>_Solve the transposed matrix equation in place._  |


## Public Functions inherited from Matrix_Corner_Block

See [Matrix\_Corner\_Block](classMatrix__Corner__Block.md)

| Type | Name |
| ---: | :--- |
|   | [**Matrix\_Corner\_Block**](classMatrix__Corner__Block.md#function-matrix_corner_block-12) (int n, int k, std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; q) <br>_A constructor for the matrix._  |
| virtual void | [**factorise**](classMatrix__Corner__Block.md#function-factorise) () override<br>_Factorise the matrix. This method prepares the matrix for a call to the solve methods. For most matrix types a call to factorise causes the LU decomposition to be calculated. The elements of the matrix should not be accessed once this method has been called._  |
| virtual double | [**get\_element**](classMatrix__Corner__Block.md#function-get_element) (int i, int j) override const<br>_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual void | [**set\_element**](classMatrix__Corner__Block.md#function-set_element) (int i, int j, double a\_ij) override<br>_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual DSpan1D | [**solve\_inplace**](classMatrix__Corner__Block.md#function-solve_inplace) (DSpan1D b) override const<br>_Solve the matrix equation in place._  |
| virtual DSpan2D | [**solve\_multiple\_inplace**](classMatrix__Corner__Block.md#function-solve_multiple_inplace) (DSpan2D bx) override const<br>_Solve multiple matrix equations in place._  |
| virtual DSpan1D | [**solve\_transpose\_inplace**](classMatrix__Corner__Block.md#function-solve_transpose_inplace) (DSpan1D b) override const<br>_Solve the transposed matrix equation in place._  |


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














## Protected Attributes

| Type | Name |
| ---: | :--- |
|  int const | [**bottom\_block\_index**](#variable-bottom_block_index)  <br>_The index of the first element of the sub-matrix I (size(A) + size(E))._  |
|  int const | [**bottom\_block\_size**](#variable-bottom_block_size)  <br>_The size of the sub-matrix I._  |
|  std::unique\_ptr&lt; double[]&gt; | [**swap\_array**](#variable-swap_array)  <br>_A memory block of size top\_block\_size which is used during swap operations._  |
|  int const | [**top\_block\_size**](#variable-top_block_size)  <br>_The size of the sub-matrix A._  |


## Protected Attributes inherited from Matrix_Corner_Block

See [Matrix\_Corner\_Block](classMatrix__Corner__Block.md)

| Type | Name |
| ---: | :--- |
|  DSpan2D | [**Abm\_1\_gamma**](classMatrix__Corner__Block.md#variable-abm_1_gamma)  <br>_The upper-right sub-matrix. The sub-matrix gamma before the_ [_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called. The sub-matrix beta after the_[_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called._ |
|  std::unique\_ptr&lt; double[]&gt; | [**Abm\_1\_gamma\_ptr**](classMatrix__Corner__Block.md#variable-abm_1_gamma_ptr)  <br>_Data storage for the upper-right sub-matrix. This memory describes gamma before the_ [_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called. This memory describes beta after the_[_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called._ |
|  [**Matrix\_Dense**](classMatrix__Dense.md) | [**delta**](classMatrix__Corner__Block.md#variable-delta)  <br>_The sub-matrix delta which is a dense matrix._  |
|  int const | [**k**](classMatrix__Corner__Block.md#variable-k)  <br>_The size of the dense corner matrix delta._  |
|  DSpan2D | [**lambda**](classMatrix__Corner__Block.md#variable-lambda)  <br>_The sub-matrix lambda._  |
|  std::unique\_ptr&lt; double[]&gt; | [**lambda\_ptr**](classMatrix__Corner__Block.md#variable-lambda_ptr)  <br>_Data storage for the sub-matrix lambda._  |
|  int const | [**nb**](classMatrix__Corner__Block.md#variable-nb)  <br>_The size of the banded matrix._  |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**q\_block**](classMatrix__Corner__Block.md#variable-q_block)  <br>_The sub-matrix Q which is a banded matrix._  |


## Protected Attributes inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
|  int const | [**n**](classMatrix.md#variable-n)  <br>_The matrix size._  |












































## Protected Functions

| Type | Name |
| ---: | :--- |
|  void | [**adjust\_indexes**](#function-adjust_indexes) (int & i, int & j) const<br>_Change the indices so they index the matrix in the corner layout._  |
|  DSpan1D | [**swap\_array\_to\_centre**](#function-swap_array_to_centre-12) (DSpan1D bx) const<br>_Rearrange the elements of the solution of the equation calculated by the_ [_**Matrix\_Corner\_Block**_](classMatrix__Corner__Block.md) _superclass so they are aligned with the way in which the user accesses the data._ |
|  DSpan2D | [**swap\_array\_to\_centre**](#function-swap_array_to_centre-22) (DSpan2D bx) const<br>_Rearrange the elements of the solutions of the equations calculated by the_ [_**Matrix\_Corner\_Block**_](classMatrix__Corner__Block.md) _superclass so they are aligned with the way in which the user accesses the data._ |
|  DSpan1D | [**swap\_array\_to\_corner**](#function-swap_array_to_corner-12) (DSpan1D bx) const<br>_Rearrange the elements of the right-hand side of the equation so they are aligned with the way in which the matrix is stored._  |
|  DSpan2D | [**swap\_array\_to\_corner**](#function-swap_array_to_corner-22) (DSpan2D bx) const<br>_Rearrange the elements of the right-hand sides of the equations so they are aligned with the way in which the matrix is stored._  |


## Protected Functions inherited from Matrix_Corner_Block

See [Matrix\_Corner\_Block](classMatrix__Corner__Block.md)

| Type | Name |
| ---: | :--- |
|   | [**Matrix\_Corner\_Block**](classMatrix__Corner__Block.md#function-matrix_corner_block-22) (int n, int k, std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; q, int lambda\_size1, int lambda\_size2) <br>_A constructor for the matrix._  |
| virtual void | [**calculate\_delta\_to\_factorise**](classMatrix__Corner__Block.md#function-calculate_delta_to_factorise) () <br>_Calculate the contents of the dense matrix_  _that will be factorised. This is the element_ _of the blockwise LU decomposition of the matrix._ |
| virtual DSpan1D | [**solve\_gamma\_section**](classMatrix__Corner__Block.md#function-solve_gamma_section) (DSpan1D const u, DView1D const v) const<br>_Calculate the solution to the following equation:_  |
| virtual DSpan1D | [**solve\_gamma\_section\_transpose**](classMatrix__Corner__Block.md#function-solve_gamma_section_transpose) (DSpan1D const v, DView1D const u) const<br>_Calculate the solution to the following equation:_  |
| virtual DSpan1D | [**solve\_lambda\_section**](classMatrix__Corner__Block.md#function-solve_lambda_section) (DSpan1D v, DView1D u) const<br>_Calculate the solution to the following equation:_  |
| virtual DSpan1D | [**solve\_lambda\_section\_transpose**](classMatrix__Corner__Block.md#function-solve_lambda_section_transpose) (DSpan1D u, DView1D v) const<br>_Calculate the solution to the following equation:_  |


## Protected Functions inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
| virtual int | [**factorise\_method**](classMatrix.md#function-factorise_method) () = 0<br>_Call the factorisation method._  |
| virtual int | [**solve\_inplace\_method**](classMatrix.md#function-solve_inplace_method) (double \* b, char transpose, int n\_equations) const = 0<br>_Call the LAPACK solve method._  |








## Detailed Description


 where E is a banded matrix, and A, E and I are square matrices.


This matrix is solved by rearranging it to:  This new matrix is a corner matrix of the form:  with  and  Internally the matrix is saved in the corner format. 


    
## Public Functions Documentation




### function Matrix\_Centre\_Block 

_A constructor for the matrix._ 
```C++
Matrix_Centre_Block::Matrix_Centre_Block (
    int n,
    int top_block_size,
    int bottom_block_size,
    std::unique_ptr< Matrix > q
) 
```





**Parameters:**


* `n` The size of the n x n [**Matrix**](classMatrix.md). 
* `top_block_size` The size of the matrix A. 
* `bottom_block_size` The size of the matrix I. 
* `q` The banded section of the matrix. 




        

<hr>



### function get\_element 

_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual double Matrix_Centre_Block::get_element (
    int i,
    int j
) override const
```





**Parameters:**


* `i` The row index. 
* `j` The column index. 



**Returns:**

The element at position (i,j). 





        
Implements [*Matrix\_Corner\_Block::get\_element*](classMatrix__Corner__Block.md#function-get_element)


<hr>



### function set\_element 

_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual void Matrix_Centre_Block::set_element (
    int i,
    int j,
    double a_ij
) override
```





**Parameters:**


* `i` The row index. 
* `j` The column index. 
* `a_ij` The value that the element should be set to. 




        
Implements [*Matrix\_Corner\_Block::set\_element*](classMatrix__Corner__Block.md#function-set_element)


<hr>



### function solve\_inplace 

_Solve the matrix equation in place._ 
```C++
virtual DSpan1D Matrix_Centre_Block::solve_inplace (
    DSpan1D b
) override const
```



Solve the following matrix equation:  The result  is saved into the memory allocated for .




**Parameters:**


* `b` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        
Implements [*Matrix\_Corner\_Block::solve\_inplace*](classMatrix__Corner__Block.md#function-solve_inplace)


<hr>



### function solve\_multiple\_inplace 

_Solve multiple matrix equations in place._ 
```C++
virtual DSpan2D Matrix_Centre_Block::solve_multiple_inplace (
    DSpan2D bx
) override const
```



Solve the following matrix equation:  for multiple values of  and . The first dimension is iterated over with each slice representing an equation to be solved. The result  is saved into the memory allocated for .




**Parameters:**


* `bx` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        
Implements [*Matrix\_Corner\_Block::solve\_multiple\_inplace*](classMatrix__Corner__Block.md#function-solve_multiple_inplace)


<hr>



### function solve\_transpose\_inplace 

_Solve the transposed matrix equation in place._ 
```C++
virtual DSpan1D Matrix_Centre_Block::solve_transpose_inplace (
    DSpan1D b
) override const
```



Solve the following matrix equation:  The result  is saved into the memory allocated for .




**Parameters:**


* `b` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        
Implements [*Matrix\_Corner\_Block::solve\_transpose\_inplace*](classMatrix__Corner__Block.md#function-solve_transpose_inplace)


<hr>
## Protected Attributes Documentation




### variable bottom\_block\_index 

_The index of the first element of the sub-matrix I (size(A) + size(E))._ 
```C++
int const Matrix_Centre_Block::bottom_block_index;
```




<hr>



### variable bottom\_block\_size 

_The size of the sub-matrix I._ 
```C++
int const Matrix_Centre_Block::bottom_block_size;
```




<hr>



### variable swap\_array 

_A memory block of size top\_block\_size which is used during swap operations._ 
```C++
std::unique_ptr<double[]> Matrix_Centre_Block::swap_array;
```




<hr>



### variable top\_block\_size 

_The size of the sub-matrix A._ 
```C++
int const Matrix_Centre_Block::top_block_size;
```




<hr>
## Protected Functions Documentation




### function adjust\_indexes 

_Change the indices so they index the matrix in the corner layout._ 
```C++
void Matrix_Centre_Block::adjust_indexes (
    int & i,
    int & j
) const
```





**Parameters:**


* `i` On input: the row index of the matrix with the banded region in the middle. On output: the row index of the matrix with the banded region in the corner. 
* `j` On input: the column index of the matrix with the banded region in the middle. On output: the column index of the matrix with the banded region in the corner. 




        

<hr>



### function swap\_array\_to\_centre [1/2]

_Rearrange the elements of the solution of the equation calculated by the_ [_**Matrix\_Corner\_Block**_](classMatrix__Corner__Block.md) _superclass so they are aligned with the way in which the user accesses the data._
```C++
DSpan1D Matrix_Centre_Block::swap_array_to_centre (
    DSpan1D bx
) const
```



I.e. for  with  a vector of length nb, return  

**Parameters:**


* `bx` The solution of the rearranged matrix equation. 



**Returns:**

The solution of the matrix equation. 





        

<hr>



### function swap\_array\_to\_centre [2/2]

_Rearrange the elements of the solutions of the equations calculated by the_ [_**Matrix\_Corner\_Block**_](classMatrix__Corner__Block.md) _superclass so they are aligned with the way in which the user accesses the data._
```C++
DSpan2D Matrix_Centre_Block::swap_array_to_centre (
    DSpan2D bx
) const
```



I.e. for  with  a vector of length nb, return  

**Parameters:**


* `bx` The solutions of the rearranged matrix equations. 



**Returns:**

The solutions of the matrix equations. 





        

<hr>



### function swap\_array\_to\_corner [1/2]

_Rearrange the elements of the right-hand side of the equation so they are aligned with the way in which the matrix is stored._ 
```C++
DSpan1D Matrix_Centre_Block::swap_array_to_corner (
    DSpan1D bx
) const
```



I.e. for  with  a vector of length nb, return  

**Parameters:**


* `bx` The right-hand side of the matrix equation to be solved. 



**Returns:**

The right-hand side of the rearranged matrix equation. 





        

<hr>



### function swap\_array\_to\_corner [2/2]

_Rearrange the elements of the right-hand sides of the equations so they are aligned with the way in which the matrix is stored._ 
```C++
DSpan2D Matrix_Centre_Block::swap_array_to_corner (
    DSpan2D bx
) const
```



I.e. for  with  a vector of length nb, return  

**Parameters:**


* `bx` The right-hand sides of the matrix equations to be solved. 



**Returns:**

The right-hand sides of the rearranged matrix equations. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_centre_block.hpp`

