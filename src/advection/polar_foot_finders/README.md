# Polar foot finder implementation details

This folder contains the implementation classes used internally by
[`PolarFootFinder`](../polar_foot_finder.hpp). They are not part of the public API.

## Structure

Each combination of `FootFindingSpace` and `AdvectionFieldSpace` requires a different
algorithm. The compile-time dispatch is handled by `ElementwiseChoice`, a partially
specialised struct that maps the two enum template parameters to a concrete type.

| `FFSpace` | `AFSpace` | Mem class | GPU functor |
|-----------|-----------|-----------|-------------|
| `LOGICAL` | `LOGICAL` | `ElementwiseLogicalAdvLogicalFootFinderMem` | `ElementwiseLogicalAdvLogicalFootFinder` |
| `PSEUDO_PHYSICAL` | `LOGICAL` | `ElementwiseLogicalAdvPseudoPhysFootFinderMem` | `ElementwiseLogicalAdvPseudoPhysFootFinder` |
| `PHYSICAL` | `LOGICAL` | `ElementwiseLogicalAdvPseudoPhysFootFinderMem` (uses physical mapping directly) | `ElementwiseLogicalAdvPseudoPhysFootFinder` |
| `PSEUDO_PHYSICAL` | `PHYSICAL` | `ElementwisePhysicalAdvPseudoPhysicalFootFinderMem` | `ElementwisePhysicalAdvPseudoPhysicalFootFinder` |
| `PHYSICAL` | `PHYSICAL` | `ElementwisePhysicalAdvPhysicalFootFinderMem` | `ElementwisePhysicalAdvPhysicalFootFinder` |

## Two-class pattern

Each algorithm is split into a pair of classes:

- **`...Mem`** — owns the spline coefficient field and a preallocated time stepper.
  Lives on the host. Calling `operator()(dt)` produces a GPU-copyable functor.
- **GPU functor** — holds only non-owning views. Can be copied to the device and
  called elementwise inside a `ddc::parallel_for_each` kernel.
