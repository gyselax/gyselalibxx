

# File mesh\_builder.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**mesh\_builder.hpp**](mesh__builder_8hpp.md)

[Go to the source code of this file](mesh__builder_8hpp_source.md)



* `#include <cstdlib>`
* `#include <ctime>`
* `#include <vector>`
* `#include "ddc_aliases.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  std::vector&lt; Coord&lt; Dim &gt; &gt; | [**build\_random\_non\_uniform\_break\_points**](#function-build_random_non_uniform_break_points) (Coord&lt; Dim &gt; min, Coord&lt; Dim &gt; max, IdxStep&lt; Grid1D &gt; n\_cells, double const non\_uniformity=1.) <br> |
|  std::vector&lt; Coord&lt; Dim &gt; &gt; | [**build\_uniform\_break\_points**](#function-build_uniform_break_points) (Coord&lt; Dim &gt; min, Coord&lt; Dim &gt; max, IdxStep&lt; Grid1D &gt; n\_cells) <br> |




























## Public Functions Documentation




### function build\_random\_non\_uniform\_break\_points 

```C++
template<class Dim, class Grid1D>
std::vector< Coord< Dim > > build_random_non_uniform_break_points (
    Coord< Dim > min,
    Coord< Dim > max,
    IdxStep< Grid1D > n_cells,
    double const non_uniformity=1.
) 
```




<hr>



### function build\_uniform\_break\_points 

```C++
template<class Dim, class Grid1D>
std::vector< Coord< Dim > > build_uniform_break_points (
    Coord< Dim > min,
    Coord< Dim > max,
    IdxStep< Grid1D > n_cells
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/mesh_builder.hpp`

