# Coding Covariant & Contravariant Tensors

## Introduction

This page explains how core mathematical concepts such as covariant and contravariant vectors and tensors are expressed in Gyselalib++.
It focuses on how these concepts are implemented using C++ templates, static typing, and compile-time structures.
We assume that you are already familiar with the mathematical and physical theory behind these ideas.
If not, please read the [Mathematical and Physical Conventions documentation](./mathematical_and_physical_conventions.md) first.

In tensor calculus, the distinction between covariant (lower index) and contravariant (upper index) components plays a critical role in how tensors transform under coordinate changes. In Gyselalib++, that distinction is encoded directly in the type system. This enables:

- Static verification of tensor operations (e.g., contraction only occurs between dual index types),

- Efficient tensor algebra with no runtime penalty,

- Self-documenting code that closely reflects mathematical intent.

The variance of an object is enforced via its type (as explained below) but it can also be deduced from the naming convention. If nothing is specified, contravariant behaviour can be assumed.
Covariant types are specified with a suffix `_cov`.

## Defining Continuous Dimensions and Variance

In tensor analysis, a vector space is defined over a basis. For example, a vector expressed in a basis $`\{\boldsymbol{\eta}_1, \boldsymbol{\eta}_2\}`$ would be written as:

$$
v = v^1 \boldsymbol{\eta}_1 + v^2 \boldsymbol{\eta}_2
$$

In this context, $`\boldsymbol{\eta}_1`$ and $`\boldsymbol{\eta}_2`$ are contravariant basis elements of a continuous vector space. The associated covariant basis elements are denoted $`\{\boldsymbol{\eta}^1, \boldsymbol{\eta}^2\}`$.
These concepts have a direct equivalent in the code.

Every basis element that may appear in a tensor expression must define the following three pieces of compile-time information:

- Whether the element can be used as a **contravariant** basis element
- Whether the element can be used as a **covariant** basis element
- What the associated basis element is in the **dual space** (the contravariant element associated with the covariant element and vice-versa).

This leads to the following code:

```cpp
struct Eta1_cov;
struct Eta1 {
    static constexpr bool IS_CONTRAVARIANT = true;
    static constexpr bool IS_COVARIANT = false;
    using Dual = Eta1_cov;
};

struct Eta1_cov {
    static constexpr bool IS_CONTRAVARIANT = false;
    static constexpr bool IS_COVARIANT = true;
    using Dual = Eta1;
};
```

Here:

- `Eta1` represents the **contravariant** basis element $`\boldsymbol{\eta}_1`$,
- `Eta1_cov` represents the **covariant dual** basis element $\boldsymbol{\eta}^1$,
- Each type declares its **dual**, linking the **contravariant** and **covariant** basis elements.

### Naming Convention

By convention:

- Contravariant dimensions are named `Eta1`, `Eta2`, `X`, `R`, etc.
- Covariant dimensions use the `_cov` suffix, such as `Eta1_cov`, `Eta2_cov`, `X_cov`, `R_cov`, etc.

### Cartesian Coordinates

In Cartesian coordinate systems, basis elements $`\{\boldsymbol{e}_i\}`$ are **self-dual**. In other words, they can be used as both contravariant and covariant bases. In such cases, it is valid to write:
$$
v = v^1e_1 + v^2e_2 = v_1e_1 + v_2e_2
$$

The corresponding code Cartesian basis elements is therefore simpler:

```cpp
struct X {
    static constexpr bool IS_CONTRAVARIANT = true;
    static constexpr bool IS_COVARIANT = true;
    using Dual = X;
};
```

This reflects the fact that $`\boldsymbol{e}_1 = \boldsymbol{e}^1`$, and no separate type is needed for the dual.

## Grouping Dimensions into a Bases

Once you have defined your continuous dimensions (e.g., `Eta1`, `Eta2`, `X`, `R`, etc.), you must group them into a set of basis elements to describe the basis on which a vector is defined.
This is done using the `VectorIndexSet` construct.

A `VectorIndexSet` represents an **ordered list of basis directions**. It is used as a template parameter to vector and tensor types and describes:

- Which basis elements are involved (e.g., `Eta1`, `Eta2`),
- In what order.

### Example: Defining Bases

```cpp
using EtaBasis = VectorIndexSet<Eta1, Eta2>;         // Contravariant basis
using EtaBasis_cov = VectorIndexSet<Eta1_cov, Eta2_cov>; // Covariant basis
```

This tells the library:

- `EtaBasis` defines a 2D contravariant space with directions `Eta1` and `Eta2`,
- `EtaBasis_cov` defines the associated covariant dual space.

You will use these sets as the type parameters for `Tensor` later.

###  Compile-Time Verification

Several compile-time tools are available to help ensure correctness and introspection of basis sets:

```cpp
static_assert(is_contravariant_vector_index_set_v<UnknownBasis>); // Ensures all elements are contravariant
static_assert(is_covariant_vector_index_set_v<UnknownBasis>);     // Ensures all elements are covariant
using UnknownBasisDual = vector_index_set_dual_t<UnknownBasis>);  // Computes the dual basis
```

Such verifications are used throughout the codebase. They ensure that tensor operations like contraction and multiplication are only performed between compatible spaces—that is, between elements with matching dual variances.

## Defining Vectors and Tensors

With your basis defined via `VectorIndexSet`, you can now create actual vector and tensor objects. In Gyselalib++, these are represented by `Tensor`. A tensor is a multidimensional array of components defined over a tensor product of basis elements.
For example, consider a 2D tensor `M` written as:
$$
M = m^i_{\;j} \boldsymbol{\eta}_i \otimes boldsymbol{b}^j
$$

This object is defined over the tensor product of the contravariant basis $`\{\boldsymbol{\eta}_i\}`$ and the covariant basis $boldsymbol{b}^j$.
Supposing the tensor components are of type `double`, this can be declared in Gyselalib++ as:

```cpp
Tensor<double, EtaBasis, BBasis_cov> M;
```

Since most tensors use `double` as their value type, the shorthand template DTensor is provided:
```cpp
DTensor<EtaBasis, BBasis_cov> M;
```

### Vectors

A vector is simply a rank-1 tensor. The same logic applies, and several equivalent declarations are supported for convenience.
For a vector expressed on the contravariant basis $`\{\boldsymbol{\eta}_1, \boldsymbol{\eta}_2\}`$, any of the following are valid:

```cpp
Tensor<double, EtaBasis> v;
DTensor<EtaBasis> v;
Vector<double, Eta1, Eta2> v;
DVector<Eta1, Eta2> v;
```
These are all equivalent and represent a 1D contravariant vector.

### Assigning Components

Component-wise assignment uses ddcHelper::select to target specific directions:

```cpp
ddcHelper::select<Eta1>(v) = 2.0;
ddcHelper::select<Eta2>(v) = 3.0;

ddcHelper::select<Eta1, B1_cov>(M) = 1.0;
ddcHelper::select<Eta1, B2_cov>(M) = 0.0;
ddcHelper::select<Eta2, B1_cov>(M) = 0.0;
ddcHelper::select<Eta2, B2_cov>(M) = 1.0;
```

This system ensures type safety: only combinations of indices consistent with the tensor’s declared basis are allowed, enabling early compile-time error detection.

## Tensor Multiplication

One of the most powerful features of Gyselalib++ is its support for tensor operations with variance-aware contraction, using Einstein notation at compile time.
Tensor contraction, such as a matrix-vector or matrix-matrix multiplication, is written using the `tensor_mul` function and `index` labels, ensuring:

- **Variance correctness**: contraction only compiles when applied between compatible covariant and contravariant indices.

- **Readability**: expressions resemble their mathematical form.

- **Performance**: evaluation is resolved statically with no runtime overhead.

### Example: Matrix–Vector Multiplication

Consider the tensor product:

$$
w^i = M^i_{\;j} \cdot v^j
$$

This is the contraction of a rank-2 tensor with a vector expressed on a contravariant basis (i.e., matrix-vector multiplication in index notation). In Gyselalib++, this is written as:

```cpp
DVector<Eta1, Eta2> w = tensor_mul(index<'i','j'>(M), index<'j'>(v));
```

Explanation:

- `M` is a `DTensor<EtaBasis, EtaBasis_cov>`, we are interested in the elements $`M^i_{\;j}`$
  - 'i' in the first dimension, i.e. the dimension associated with the contravariant basis $`\{\boldsymbol{\eta}_i\}`$
  - 'j' in the second dimension, i.e. the dimension associated with the covariant basis $`\{\boldsymbol{\eta}^i\}`$
- `v` is a `DVector<EtaBasis>`, we are interested in the elements $v^j$
  - 'j' in the first (only) dimension, i.e. the dimension associated with the contravariant basis $`\{\boldsymbol{\eta}_i\}`$
- The result is a contravariant vector (`DVector<Eta1, Eta2>`) thanks to the Einstein summation we know that we are interested in the elements $w^i$.

The indices used in `index<'i','j'>` are **compile-time labels** that identify how dimensions are mapped and contracted.

#### Einstein Summation and Compile-Time Matching

This syntax uses Einstein summation convention, which means:

- Any repeated index (like 'j' here) implies a summation over that index.

- All unrepeated (free) indices appear in the result — here, 'i'.

In this case:

- The 'j' index appears once as covariant index (in M) and once as contravariant index (in v), forming a valid contraction pair.

- The resulting tensor w has the free index 'i', matching the contravariant structure of the output vector.

The system performs strict compile-time checks:

- Ensuring contraction only occurs between dual variance pairs (e.g., contravariant ↔ covariant),

- Rejecting operations like contracting two covariant indices, or producing mismatched output types.

This not only enforces mathematical correctness but ensures that operations reflect the physical and geometrical structure of the problem being modelled.

### Example: Matrix–Matrix Multiplication

Now consider multiplying two tensors:

$$
P_i^{\;k} = M_{ij} \cdot N^{jk}
$$

This represents a standard matrix–matrix multiplication via contraction over index $j$.

In code:

```cpp
DTensor<EtaBasis_cov, EtaBasis_cov> M;
DTensor<EtaBasis, EtaBasis> N;
DTensor<EtaBasis_cov, EtaBasis> P = tensor_mul(index<'i','j'>(M), index<'j','k'>(N));
```

#### Einstein Summation and Compile-Time Matching

As in the Vector-Matrix case Einstein summation is used. To summarise:

- `index<'i','j'>(M)` specifies how to label the dimensions of `M`
- `index<'i','j'>(M)` specifies how to label the dimensions of `N`
- Since `'j'` appears once in each tensor, and each occurrence involves **dual variance types**, a **valid contraction** is performed
- The resulting tensor `P` has dimensions `'i'` and `'k'`, both preserved from the inputs as they are not repeated
- Variance is inferred based on the types of the basis elements used in each dimension, here `'i'` represents elements of the covariant basis on both sides of the equation, while `'k'` represents elements of the contravariant basis.

## Vector Fields

In many physical applications (including plasma physics) fields are not just scalar quantities, but vector-valued functions defined over a spatial domain. For example, the electric field $\mathbf{E}(\mathbf{x})$ associates a vector to each point in space.

In Gyselalib++, this is represented using `VectorField`, a field whose values are vectors, defined at each point on a discretised grid.

You can think of a DVectorField as a Field where each element is a DVector.

:warning: If you’re not yet familiar with scalar `Field<T, Idx>` types, refer to the documentation about [DDC in Gyselalib++](../first_steps/DDC_in_gyselalibxx.md) for an introduction.

### Syntax

```cpp
DVectorField<
    IdxRange<GridR, GridTheta>, // Discretised grid as for a scalar field
    CartesianBasis              // The basis on which the vectors are defined
> E;
```

This corresponds, for example, to defining a contravariant electric field $`\boldsymbol{E}(\boldsymbol{r}) = E^1(r_i, \theta_j) \boldsymbol{e}_1 + E^2(r_i, \theta_j) \boldsymbol{\e}_2 = E_x(r_i, \theta_j) \hat{x} + E_y(r_i, \theta_j) \hat{y}`$ on a polar mesh.

### Component Access and Assignment

The scalar fields in the vector fields can be extracted using `ddcHelper::select`. This allows components to be assigned as for scalar fields :

```cpp
IdxEta1 i;
ddcHelper::select<X>(E)(i) = 4.2; // Assign to component E_x at index i without a temporary
DField<IdxRange<GridR, GridTheta>> E_x = ddcHelper::select<X>(E);
E_x(i) = -0.3;                    // Assign to component E_x at index i via a temporary field
```

:warning: Currently it is only possible to get a constant vector from a vector field so it is not possible to directly assign a vector to an element of a vector field. The method `assign_vector_field_element` automatically carries out the elementwise assignment.

## Changing Vector Space

In many simulations, you may need to **transform a vector from one basis to another** — for example:

- Switching between covariant and contravariant representations
- Applying a coordinate transformation (e.g., from Cartesian to curvilinear basis)

This is done using the `to_vector_space` function and a coordinate transformation class.

### Example: Converting Covariant to Contravariant

Suppose you have a **covariant vector field** (e.g., a gradient) and want to convert it into a **contravariant representation**. This calculation requires a metric tensor. The metric tensor is defined from a coordinate transformation to/from a Cartesian basis.

Elementwise the operation is carried out at a given coordinate (as the metric tensor may depend on the coordinate:

```cpp
using PolarBasis = VectorIndexSet<R, Theta>;
Coord<R, Theta> coord; // The coordinate where d phi was calculated
DVector<R_cov, Theta_cov> d_phi; // The derivative with respect to phi
CircularToCartesian<R, Theta, X, Y> coord_transform;
DVector<R, Theta> grad_phi = to_vector_space<EtaBasis>(coord_transform, coord, d_phi);
```

For a vector field the code is:

```cpp
DVectorField<IdxRange<GridR, GridTheta>, VectorIndexSet<R_cov, Theta_cov>> d_phi;
DVectorField<IdxRange<GridR, GridTheta>, VectorIndexSet<R, Theta>> grad_phi; 
copy_to_vector_space<PolarBasis>(Kokkos::DefaultExecutionSpace(), grad_phi, coord_transform, d_phi);
```

### Example: Converting Between Bases

Converting between bases works very similarly to converting between covariant and contravariant bases. This time though it is important that the coordinate transformation maps between the coordinates that we are switching from/to.

```cpp
using PolarBasis = VectorIndexSet<R, Theta>;
using EtaBasis = VectorIndexSet<Eta1, Eta2>;
Coord<R, Theta> coord; // The coordinate where the vector was calculated
DVector<R_cov, Theta_cov> E;
CircularToCartesian<R, Theta, Eta1, Eta2> coord_transform;
DVector<Eta1, Eta2> grad_phi = to_vector_space<EtaBasis>(coord_transform, coord, d_phi);
```

## Static Tensors and Compile-Time Optimisation

Some tensors, like the identity tensor or the Levi-Civita symbol, are independent of runtime coordinates and can be fully defined at compile time. These are available through implementations in `static_tensors.hpp`.

Using static tensors enables:

- Efficient algebra with no runtime allocation or lookup

- Compile-time elimination of zero terms (e.g., sparse tensors, skew-symmetry)

- Clean expression of tensor identities and symmetries

For example, a cross product can be expressed with the Levi-Civita tensor as:

```cpp
CartesianLeviCivitaTensor<double, X, Y, Z> eps;
DVector<X, Y, Z> v;
DVector<X, Y, Z> w;
DVector<X, Y, Z> z = tensor_mul(index<'i','j','k'>(eps), index<'i'>(v), index<'j'>(w));
```

Thanks to inlining and the Einstein summation notation, the compiler understands the last line to mean:

```cpp
ddcHelper::select<X>(z) = 
          0*ddcHelper::select<X>(v)*ddcHelper::select<X>(w)  // i=X, j=X, k=X
        + 0*ddcHelper::select<X>(v)*ddcHelper::select<Y>(w)  // i=X, j=X, k=Y
        + 0*ddcHelper::select<X>(v)*ddcHelper::select<Z>(w)  // i=X, j=X, k=Z
        + 0*ddcHelper::select<Y>(v)*ddcHelper::select<X>(w)  // i=X, j=Y, k=X
        + 0*ddcHelper::select<Y>(v)*ddcHelper::select<Y>(w)  // i=X, j=Y, k=Y
        + 1*ddcHelper::select<Y>(v)*ddcHelper::select<Z>(w)  // i=X, j=Y, k=Z
        + 0*ddcHelper::select<Z>(v)*ddcHelper::select<X>(w)  // i=X, j=Z, k=X
        - 1*ddcHelper::select<Z>(v)*ddcHelper::select<Y>(w)  // i=X, j=Z, k=Y
        + 0*ddcHelper::select<Z>(v)*ddcHelper::select<Z>(w); // i=X, j=Z, k=Z
ddcHelper::select<Y>(z) = ...;
ddcHelper::select<Z>(z) = ...;
```

As long as at least `-O1` optimistations are activated, the compiler will optimise away the multiplications by 0, effectively turning this into the standard form of a cross product:

```cpp
ddcHelper::select<X>(z) = 
          ddcHelper::select<Y>(v)*ddcHelper::select<Z>(w) // i=X, j=Y, k=Z
        - ddcHelper::select<Z>(v)*ddcHelper::select<Y>(w); // i=X, j=Z, k=Y
ddcHelper::select<Y>(z) = ...;
ddcHelper::select<Z>(z) = ...;
```

