

# Class Matrix\_Periodic\_Banded



[**ClassList**](annotated.md) **>** [**Matrix\_Periodic\_Banded**](classMatrix__Periodic__Banded.md)



_A class representing a periodic banded matrix. A periodic banded matrix is like a banded matrix but additionally contains non- zero values in the corners. I.e it has the following sparsity pattern:_ [More...](#detailed-description)

* `#include <matrix_periodic_banded.hpp>`



Inherits the following classes: [Matrix\_Corner\_Block](classMatrix__Corner__Block.md)










































































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Matrix\_Periodic\_Banded**](#function-matrix_periodic_banded) (int n, int kl, int ku, std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; q) <br>_A constructor for the matrix._  |
| virtual double | [**get\_element**](#function-get_element) (int i, int j) override const<br>_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual void | [**set\_element**](#function-set_element) (int i, int j, double a\_ij) override<br>_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._ |


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
|  int const | [**kl**](#variable-kl)  <br>_Number of subdiagonals._  |
|  int const | [**ku**](#variable-ku)  <br>_Number of superdiagonals._  |


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
| virtual void | [**calculate\_delta\_to\_factorise**](#function-calculate_delta_to_factorise) () override<br>_Calculate the contents of the dense matrix_  _that will be factorised. This is the element_ _of the blockwise LU decomposition of the matrix._ |
| virtual DSpan1D | [**solve\_lambda\_section**](#function-solve_lambda_section) (DSpan1D v, DView1D u) override const<br>_Calculate the solution to the following equation:_  |
| virtual DSpan1D | [**solve\_lambda\_section\_transpose**](#function-solve_lambda_section_transpose) (DSpan1D u, DView1D v) override const<br>_Calculate the solution to the following equation:_  |


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


 


    
## Public Functions Documentation




### function Matrix\_Periodic\_Banded 

_A constructor for the matrix._ 
```C++
Matrix_Periodic_Banded::Matrix_Periodic_Banded (
    int n,
    int kl,
    int ku,
    std::unique_ptr< Matrix > q
) 
```





**Parameters:**


* `n` The size of the n x n matrix. 
* `kl` The number of lower diagonals. 
* `ku` The number of upper diagonals. 
* `q` The sub-matrix Q describing the banded sub-matrix. 




        

<hr>



### function get\_element 

_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual double Matrix_Periodic_Banded::get_element (
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
virtual void Matrix_Periodic_Banded::set_element (
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
## Protected Attributes Documentation




### variable kl 

_Number of subdiagonals._ 
```C++
int const Matrix_Periodic_Banded::kl;
```




<hr>



### variable ku 

_Number of superdiagonals._ 
```C++
int const Matrix_Periodic_Banded::ku;
```




<hr>
## Protected Functions Documentation




### function calculate\_delta\_to\_factorise 

_Calculate the contents of the dense matrix_  _that will be factorised. This is the element_ _of the blockwise LU decomposition of the matrix._
```C++
virtual void Matrix_Periodic_Banded::calculate_delta_to_factorise () override
```



 It is defined as:





With  the solution of the equation:





which should already be stored in the variable Abm\_1\_gamma when this method is called. 


        
Implements [*Matrix\_Corner\_Block::calculate\_delta\_to\_factorise*](classMatrix__Corner__Block.md#function-calculate_delta_to_factorise)


<hr>



### function solve\_lambda\_section 

_Calculate the solution to the following equation:_ 
```C++
virtual DSpan1D Matrix_Periodic_Banded::solve_lambda_section (
    DSpan1D v,
    DView1D u
) override const
```






The result is saved into the argument . This calculation is one of the steps necessary to solve the matrix equation.




**Parameters:**


* `v` A vector of length k describing the final values in the right-hand side of the matrix to be solved. 
* `u` A vector of length (n-k). This is an intermediate result of the calculation which solves the matrix equation.



**Returns:**

The vector v solution of the equation. 





        
Implements [*Matrix\_Corner\_Block::solve\_lambda\_section*](classMatrix__Corner__Block.md#function-solve_lambda_section)


<hr>



### function solve\_lambda\_section\_transpose 

_Calculate the solution to the following equation:_ 
```C++
virtual DSpan1D Matrix_Periodic_Banded::solve_lambda_section_transpose (
    DSpan1D u,
    DView1D v
) override const
```






The result is saved into the argument . This calculation is one of the steps necessary to solve the matrix equation.




**Parameters:**


* `u` A vector of length (n-k). This is an intermediate result of the calculation which solves the transposed matrix equation. 
* `v` A vector of length k. This is an intermediate result of the calculation which solves the transposed matrix equation.



**Returns:**

The vector u solution of the equation. 





        
Implements [*Matrix\_Corner\_Block::solve\_lambda\_section\_transpose*](classMatrix__Corner__Block.md#function-solve_lambda_section_transpose)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_periodic_banded.hpp`

