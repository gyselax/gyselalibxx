# Multipatch connectivity 

This directory defines structures and methods to describe patches and how they are connected to one another.

The following sections decribe the structures and methods implemented: 

-  [Patch](#src_multipatch_connectivity__Patch) - Definition of Patch tag.
-  [Interfaces](#src_multipatch_connectivity__Interfaces)
    - [Sticking of Two Edges](#src_multipatch_connectivity__Interfaces__Sticking_of_Two_Edges)
        - [Multi-patch domain](#src_multipatch_connectivity__Interfaces__Multi-patch_domain) - Mathematical definition of multi-patch domain.
        - [Edges](#src_multipatch_connectivity__Interfaces__Edges) - Mathematical defintion of edges. 
        - [Sticking and Coordinate Transformation](#src_multipatch_connectivity__Interfaces__Sticking_and_Coordinate_Transformation) - Mathematical defintion of  coordinate transformation and links to the implemented class (Interface, Edge and EdgeTransformation). 
        - [Index transformation](#src_multipatch_connectivity__Interfaces__Index_transformation) - Algorithm of the index transformation in EdgeTransformation. 
        - [Conformity of the meshes](#src_multipatch_connectivity__Interfaces__Conformity_of_the_meshes) - Definition of `UniformGridIdxMatching`. 
- [References](#src_multipatch_connectivity__References) - References. 
- [Contents](#src_multipatch_connectivity__Contents) - List of files in the folder. 


## Patch
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


## Interfaces 
### Sticking of Two Edges

We follow the idea to represent the topology of a multipatch domain using
the idea of *domain manifolds* (see [1]). This means that we have several
independent tensor-product logical domains and define their sticking via coordinate 
transformations.


**Remark:** In the referenced paper, the authors use affine isometries of the entire 
patch instead of coordinate transformations of the faces. This is equivalent to our approach in 
the sense that, when the scaling of the domains is removed and we consider only 
unit squares as logical domains, the coordinate transformations of the faces and the information 
about which edges are identified determine the 
isometries used in the paper and vice versa. This is easy to show.

#### Multi-patch domain

For simplicity, we will constrain ourselves to the 2D case, but this approach
could be generalized to arbitrary dimensions.

Let $\Omega$ be the domain of interest and assume that we have patches $\Omega^{(i)}$, $i=1,...,K$, which are disjoint, s.t.
```math
\Omega = \dot{\bigcup}_{i=1}^K \Omega^{(i)}.
```

For every patch $\Omega^{(i)}$ we have a mapping 
```math
F^{(i)}: [a_x^{(i)}, b_x^{(i)}]\times[a_y^{(i)} b_y^{(i)}] \rightarrow \Omega^{(i)}, 
```

where the rectangular domain $`[a_x^{(i)}, b_x^{(i)}]\times[a_y^{(i)} b_y^{(i)}]`$ defines 
the logical coordinates for the patch. We call this logical coordinate domain 
the *logical patch*.

#### Edges

So we obtain four logical edges
```math 
[a_x^{(i)}, b_x^{(i)}] \times \{ a_y^{(i)} \}, \quad
[a_x^{(i)}, b_x^{(i)}] \times \{ b_y^{(i)} \}, \quad
\{ a_x^{(i)} \} \times [a_y^{(i)}, b_y^{(i)}], \quad
\{ b_x^{(i)} \} \times [a_y^{(i)}, b_y^{(i)}] 
```

In the code, we define edges as follows. 
Every edge of a logical patch is identified via the patch it belongs to, the dimension 
and whether it is at the front or the back of the domain. 
So e.g.the edge $`[a_x^{(i)}, b_x^{(i)}] \times \{ a_y^{(i)} \}`$ would be identified with patch $i$, dimensions `Yi` and `FRONT`.
$`[a_x^{(i)}, b_x^{(i)}] \times \{ b_y^{(i)} \}`$ would be identified with patch $i$, dimensions `Yi` and `BACK` and 
$`\{ b_x^{(i)} \} \times [a_y^{(i)}, b_y^{(i)}]`$ would be identified with patch $i$, dimensions `Xi` and `BACK`.

#### Sticking and Coordinate Transformation

Any edge can then be associated with an edge on another patch. 
This corresponds to the 'sticking'. 
So for a different patch $`\Omega^{(j)}`$, we have the logical coordinate domain
$`[a_x^{(j)}, b_x^{(j)}] \times [a_y^{(j)} b_y^{(j)}]`$.
If we want to stick patch $i$ to patch $j$, we have to determine the edges that
are identified and how they are identified. 
The way that they are identified is mathematically determined via the coordinate transformation from one edge
to the other. 
Since the transformations are supposed to be affine and bijective, there are only two options: 
the transformation can be order preserving or 
order reversing (this corresponds to the orientation of the phyiscal edge where two parametrizations 
coming from the two patches can have either the same or the opposite orientation respectively).

So for example, we want to stick the edge $`\{ a_x^{(i)} \} \times [a_y^{(i)}, b_y^{(i)}]`$
on patch $i$ to the edge $`[a_x^{(j)}, b_x^{(j)}] \times \{ b_y^{(j)} \}`$ on patch $j$. 
If the transformation is order-preserving (i.e. the orientations of the parametrizations 
of the physical edge agree), then the transformation from the first edge to the second is 
```math
t \mapsto a_x^{(j)} + \frac{t - a_y^{(i)}}{b_y^{(i)} - a_y^{(i)}} \, (b_x^{(j)} - a_x^{(j)}).
```

If the transformation is order-reversing (i.e. the orientations of the parametrizations 
of the physical edge are opposite), then it is
```math
t \mapsto b_x^{(j)} - \frac{t - a_y^{(i)}}{b_y^{(i)} - a_y^{(i)}} \, (b_x^{(j)} - a_x^{(j)}).
```


In the code, the sticking of two edges is represented by an `Interface` structure which contains tags 
to the first patch and the second patch as well as the boolean member `orientations_agree`. 
The transformation from one edge to the other is done using the `EdgeTransformation` operator. 

This operator contains an `operator()` to transform a given coordinate on a coordinate on the other patch. 
It has also been extended for indices. 

#### Index transformation

All the patches are discretized into a finite number of points indexed by indices. 
So each index refers to a coordinate on a patch. 
For a given index on an edge of a patch at the interface, we can find an equivalent 
index on the other patch. This is always true for conforming meshes but not always the 
case for non-conforming meshes. 

For non-conforming meshes, the method `is_match_available()` is needed to identify if the given 
index has an equivalent or not. If `is_match_available()` returns true, we can then 
call the `operator()` to get the equivalent index on the other patch. 
Instead of calling both method, it is also possible to call `search_for_match()` which 
returns true if there is an equivalent index and updates it the given target index. 

**For uniform grids** $`\{x_i\}_{i=0, ..., N}`$ (current grid)  and $`\{x_j\}_{j=0, ..., M}`$ (target grid), 
finding an equivalent index can be done with the following steps:

* We compute the greatest common divisor of the number of cells $`N`$ and $`M`$ 
(e.g. $`\exists \ gcd \in\mathbb{N}, gcd>1; \ M \wedge N = gcd`$).

* The current index $`k`$ can have an equivalent index on the target patch. 

Depending on if $`k / p \in \mathbb{N}^*`$ with $`p = \frac{N}{gcd} \in \mathbb{N}^*`$,
this equivalent is given by multiplying by $`\frac{M}{N}`$ and adapt if the coordinate transformation is not order-preserving. 

* If the greatest common divisor is $`1`$, there is no equivalent index except for the first and last indices. 



**For non-uniform grids**, we cannot use this algorithm. We check all the indices of the target patch if 
one of them fit with the current one. 
We can use dichotomy method for the search. 



#### Conformity of the meshes
When there are non-conforming meshes, some operators need to loop on the "conforming indexes" at each Interface, 
i.e. the indexes with an equivalent index on the other patch of the interface. 

The `UniformIdxStepIndexMatching` storage class stores for a given Interface the "conforming indexes" in an IdxRangeSlice for each patch. 
Currently, the IdxRangeSlice stores indexes from an IdxRange with regular index steps. 
So the `UniformIdxStepIndexMatching` class does the same, and triggers assert when the conformity between two edges is not 
uniform. 

For example, `UniformIdxStepIndexMatching` is adapted to this type of interface
```
    0   1   2   3   4   5   6   7   8       <- Edge1
    |   |   |   |   |   |   |   |   |
       -2-     -2-     -2-     -2-          <- The index steps are regular. 
       -1-     -1-     -1-     -1-          <- The index steps are regular. 
    |       |       |       |       |
    0       1       2       3       4       <- Edge2
``` 

but fails for this type of interface
```
    0   1   2   3   4   5   6   7   8       <- Edge1
    |   |   |   |   |   |   |   |   |
     -1-     -3-       -2-     -2-          <- The index steps are NOT regular. 
     -1-     -1-       -1-     -1-          <- The index steps are regular. 
    |   |           |       |       |
    0   1           2       3       4       <- Edge2
``` 

The grids can be uniform or not uniform, as long as the index steps between the conforming indices are regular. 

If the grids are uniform, we know that the index steps between the conforming indices are regular. 
Moreover the index steps are given by the number of cells divided by the greatest common divisor of the 
two numbers of cells. 
(E.g. if the first example was on uniform grids: $`\gcd(8,4) = 4`$. On the first edge, the index step is $`\frac{8}{4} = 2`$ and 
on the second edge, $`\frac{4}{4} = 1`$.) 



## References
[1] Buchegger, F., Jüttler, B., Mantzaflaris, A., 
*Adaptively refined multi-patch B-splines with enhanced smoothness*.
Applied Mathematics and Computation,
Volume 272, Part 1
(2016).
https://www.sciencedirect.com/science/article/pii/S009630031500836X.


## Contents
- `edge_transformation.hpp`: Defines the `EdgeTransformation` operator to transform the coordinate from one edge to the other (see above).
- `edge.hpp`: Defines the `Edge` structure.
- `interface.hpp`: Defines the `Interface` structure.
- `patch.hpp`: Defines the `Patch` tag.
