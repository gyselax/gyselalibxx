# Mappings

This folder contains code describing tools for handling different coordinate systems.

The current coordinate transformations implemented are:

## Analytically invertible coordinate transformations

### Circular coordinate transformation

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

### Czarny coordinate transformation

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

### Barycentric coordinate transformation

- Mapping (BarycentricToCartesian):

```math
 (c_1, c_2, c_3) \rightarrow (x, y).
```

## Discrete coordinate transformation defined on B-splines

- Mapping (DiscreteToCartesian)

```math
\left\{
\begin{aligned}
 & x(r,\theta) = \sum_k c_{x,k} B_k(r,\theta) , \\
 & y(r,\theta) = \sum_k c_{y,k} B_k(r,\theta) .
\end{aligned}
\right.
```

## 3D coordinate transformations

### Cylindrical to Cartesian transformation

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

### Toroidal to Cylindrical transformation

The toroidal to cylindrical transformation combines a 2D transformation in the polar plane of the torus with an orthogonal transformation in the remaining dimension.

```math
\begin{aligned}
(\rho, \theta) = f(R, Z) \\
\phi = - \zeta
\end{aligned}
```

where $f$ is a 2D curvilinear to cartesian transformation.

## Combined coordinate transformation which combines two of the coordinate transformations above

The tools are:

- InverseJacobianMatrix : this tool calculates the inverse Jacobian matrix on the specified coordinate system.
- InvJacobianOPoint : this tool calculates the inverse Jacobian matrix at the O-point on the specified coordinate system.
- MetricTensorEvaluator : this tool calculates the metric tensor associated with a coordinate transformation.
- VectorMapper : this tool helps when converting vectors stored in a `VectorField` from one coordinate system to another.
- other static analysis tools found in `mapping_tools.hpp`
