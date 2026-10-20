# Mappings

This folder contains code describing tools for handling different coordinate systems.

The current coordinate transformations implemented are:

## Analytically invertible coordinate transformations

### 1D coordinate transformations

#### Linear coordinate transformation

```math
y(x) = \alpha x + \beta
```

### 2D coordinate transformations

#### Circular coordinate transformation

- Mapping (CircularToCartesian):

```math
\left\{
\begin{aligned}
 & x(r,\theta) = r \cos(\theta), \\
 & y(r,\theta) = r \sin(\theta).
\end{aligned}
\right.
```

```math
J (r,\theta) = 
\begin{bmatrix} 
 \cos(\theta) & -r\sin(\theta) \\
 \sin(\theta) &  r\cos(\theta)
\end{bmatrix}.
```

with $\det(J) = r$.

- Inverse mapping (CartesianToCircular):

```math
\left\{
\begin{aligned}
 & r(x,y) = \sqrt{x^2 + y^2}, \\
 & \theta(x,y) = \text{atan2}(y, x).
\end{aligned}
\right.
```

#### Czarny coordinate transformation

- Mapping (CzarnyToCartesian): Taking $\varepsilon \in ]0,1[$ and $e$ as parameters.

```math
\left\{
\begin{aligned}
 & x(r,\theta) = \frac{1}{\varepsilon} \left( 1 - \sqrt{1 + \varepsilon(\varepsilon + 2 r \cos(\theta))} \right), \\
 & y(r,\theta) = \frac{e\xi r \sin(\theta)}{2 -\sqrt{1 + \varepsilon(\varepsilon + 2 r \cos(\theta))}}.
\end{aligned}
\right.
```

with $\xi = 1/\sqrt{1 - \varepsilon^2 /4}$.

```math
J (r,\theta) 
%= 
% \begin{bmatrix} 
%  -\frac{\cos(\theta)}{\sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}} 
%  & 
%  \frac{r\sin(\theta)}{\sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}} 
%  \\
%  \frac{\cos(\theta)}{\sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}} 
%  \frac{e\varepsilon \xi r \sin(\theta)}{\left(2 - \sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}\right)^2}
%  + 
%  \frac{e \xi \sin(\theta)}{2 - \sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}}
%  & 
%  \frac{-r\sin(\theta)}{\sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}} 
%  \frac{e\varepsilon \xi r \sin(\theta)}{\left(2 - \sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}\right)^2}
%  + 
%  \frac{e \xi r\cos(\theta)}{2 - \sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}}
%  \\
% \end{bmatrix}
= 
\frac{1}{\sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}}
\begin{bmatrix} 
 -1
 & 
 0 
 \\
 \frac{e\varepsilon \xi r \sin(\theta)}{\left(2 - \sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}\right)^2}
 & 
 \frac{\sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}}{2 - \sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}}
 \\
\end{bmatrix}
\begin{bmatrix} 
 \cos(\theta) & -r\sin(\theta) \\
 \sin(\theta) &  r\cos(\theta) \\
\end{bmatrix}.
```

with $\det(J) = \frac{-r}{\sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}}\frac{e \xi}{2 - \sqrt{1 + \varepsilon(\varepsilon + 2r\cos(\theta))}}$.

- Inverse mapping (CartesianToCzarny): Taking $\varepsilon \in ]0,1[$ and $e$ as parameters.

```math
\left\{
\begin{aligned}
 & r(x,y) = \sqrt{\frac{y^2 (1+\varepsilon x)^2}{e^2\xi^2+0.25(\varepsilon x^2-2x-\varepsilon)^2}}, \\
 & \theta(x,y) = \text{atan2}(2y (1+\varepsilon x), e \xi (\varepsilon x^2 - 2x-\varepsilon)).
\end{aligned}
\right.
```

### 3D coordinate transformations

#### Cylindrical to Cartesian transformation

- Mapping (CylindricalToCartesian):

```math
\left\{
\begin{aligned}
 & x(R,Z,\zeta) = R \cos(\zeta), \\
 & y(R,Z,\zeta) = R \sin(\zeta), \\
 & z(R,Z,\zeta) = Z
\end{aligned}
\right.
```

```math
J (R, Z, \zeta) = 
\begin{bmatrix} 
 \cos(\zeta) & 0 & -R\sin(\zeta) \\
 \sin(\zeta) & 0 &  R\cos(\zeta)
  0          & 1 &       0 \\
\end{bmatrix}.
```

with $\det(J) = -R$.

- Inverse mapping (CartesianToCylindrical):

```math
\left\{
\begin{aligned}
 & R(x,y,z) = \sqrt{x^2 + y^2}, \\
 & \zeta(x,y,z) = \text{atan2}(y, x), \\
 & Z(x,y,z) = z.
\end{aligned}
\right.
```

#### Toroidal to Cylindrical transformation

The toroidal to cylindrical transformation combines a 2D transformation in the polar plane of the torus with an orthogonal transformation in the remaining dimension.

```math
\begin{aligned}
(\rho, \theta) = f(R, Z) \\
\phi = - \zeta
\end{aligned}
```

where $f$ is a 2D curvilinear to cartesian transformation.

```math
J (R, Z, \zeta) =
\begin{bmatrix}
 J_f\.^{\rho}_{\;R} & J_f\.^{\rho}_{\;Z} & 0 \\
 J_f\.^{\theta}_{\;R} & 0 & J_f\.^{\theta}_{\;Z} & 0
  0          & 0 &       -1 \\
\end{bmatrix}.
```

## Barycentric coordinate transformation

- Mapping (BarycentricToCartesian):

```math
 (c_1, c_2, c_3) \rightarrow (x, y).
```

## Discrete coordinate transformation defined on B-splines

- Mapping (DiscretePoloidalCSSplineMapping)

```math
\left\{
\begin{aligned}
 & x(r,\theta) = \sum_k c_{x,k} B_k(r,\theta) , \\
 & y(r,\theta) = \sum_k c_{y,k} B_k(r,\theta) .
\end{aligned}
\right.
```

## Combined coordinate transformation which combines two of the coordinate transformations above

The tools are:

- InverseJacobianMatrix : this tool calculates the inverse Jacobian matrix on the specified coordinate system.
- InvJacobianOPoint : this tool calculates the inverse Jacobian matrix at the O-point on the specified coordinate system.
- MetricTensorEvaluator : this tool calculates the metric tensor associated with a coordinate transformation.
- VectorMapper : this tool helps when converting vectors stored in a `VectorField` from one coordinate system to another.
- other static analysis tools found in `coord_transformation_tools.hpp`

## Orthogonal coordinate transformation

A coordinate transformation that is comprised of coordinate transformation components which are orthogonal to one another.

## Changing the vector space of a vector

The function `to_vector_space` (in `vector_mapper.hpp`) expresses a vector $v$ on a different vector space.
There are two kinds of change:

- A change of variance (contravariant ↔ covariant) in the same coordinate system.
  This takes a metric tensor $G$ (e.g. a `MetricTensorEvaluator`).
- A change of coordinate system with the same variance.
  This takes a mapping $\mathcal{F} : q \mapsto p$, where $q$ are the argument coordinates and $p$ are the result coordinates.
  $J$ denotes its Jacobian matrix.

| Case | Takes | Computation |
|---|---|---|
| Same vector space | Mapping or metric tensor | $v$ is returned unchanged |
| Contravariant → covariant, same coordinate system | Metric tensor | $G v$ |
| Covariant → contravariant, same coordinate system | Metric tensor | $G^{-1} v$ |
| Argument basis → result basis, contravariant | Mapping | $J v$ |
| Argument basis → result basis, covariant | Mapping | $J^{-T} v$ |
| Result basis → argument basis, contravariant | Mapping | $J^{-1} v$ |
| Result basis → argument basis, covariant | Mapping | $J^T v$ |

A Cartesian basis is both covariant and contravariant, so a change to or from the Cartesian space is always a change of coordinate system.
To change both the coordinate system and the variance, call `to_vector_space` twice: once with the mapping and once with the metric tensor.

The functions `copy_to_vector_space` and `create_mirror_view_and_copy_on_vector_space` apply the same operations to every element of a vector field.
They also take either a mapping or a metric tensor.
