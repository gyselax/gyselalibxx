# The discrete domain computation library (DDC)
<!--
Copyright (C) The DDC development team, see COPYRIGHT.md file

SPDX-License-Identifier: MIT
-->

See https://ddc.mdls.fr/

[DDC](https://ddc.mdls.fr/), is a C++-17 library that aims to offer to the C++/MPI world an equivalent to the [`xarray.DataArray`](https://xarray.pydata.org/en/stable/generated/xarray.DataArray.html)/[`dask.Array`](https://docs.dask.org/en/stable/array.html) python environment.
Where these two libraries are based on [numpy](https://numpy.org/), DDC relies on [Kokkos](https://github.com/kokkos/kokkos) and [mdspan](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2020/p0009r10.html) to offer CPU/GPU performance-portable multi-dimensional arrays and iterators.

DDC aims to offer abstractions to represent:
* tagged continuous dimensions,
* discretizations of these (multiple kinds of meshes, function spaces, Fourier, etc.),
* domains in these discretized dimensions,
* associating data to the elements of these discrete domains,
* efficient distribution and iteration over this data.

All these abstractions are handled at compilation using C++ template meta-programming to ensure zero overhead at execution and enable developers to change their design choices (*e.g.* from a regular mesh to a non-uniform one) with performance portability.

[Find out more...](https://ddc.mdls.fr/)

If you like the project, please leave us a github star.

If you want to know more, join un on [Slack](https://join.slack.com/t/ddc-lib/shared_invite/zt-14b6rjcrn-AwSfM6_arEamAKk_VgQPhg)
