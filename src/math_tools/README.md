# Utility Functions

This folder contains mathematical classes and functions.

## Derivative tools

Functions for calculating derivatives with different methods:

- `spline_1d_partial_derivative.hpp` : calculate derivatives using a 1d spline representation.
- `central_fdm_partial_derivatives.hpp` : calculate derivatives using Finite Difference Method.

### Finite difference derivatives

The method is as follow:
Take a function $f$. We want to approximate the value of $`f'(x_1)`$ knowing the value of $`f(x_1)`$, $`f(x_2)`$ and $`f(x_3)`$ for
$`x_3>x_2>x_1`$.
Denote $`\alpha :=|x_2-x_1|`$ and $`\beta :=|x_3-x_2|`$.
We call $`Df(x_1)`$ the approximate value of $`f'(x_1)`$. We want $`|f'(x_1)-Df(x_1)|=o\left(\max(\alpha,\beta\right)^2)`$.
To this end, write

```math
 Df(x_1)=c_1f(x_1)+c_2f(x_2)+c_3f(x_3),
```

and expand $f$ to the second order.

```math
    \begin{aligned}
        Df(x_1) & = c_1f(x_1)+c_2\left(f(x_1)+\alpha f'(x_1)+\frac{1}{2}\alpha^2f''(x_1)+o(\alpha^2)\right)               \\
                & +c_3\left(f(x_1)+(\alpha+\beta) f'(x_1)+\frac{1}{2}(\alpha+\beta)^2f''(x_1)+o((\alpha+\beta)^2)\right).
    \end{aligned}
```

A term-by-term comparison gives us the following system:

```math
    \left\{
    \begin{aligned}
        c_1+c_2+c_3                      & =0  \\
        \alpha c_2+(\alpha+\beta)c_3     & =1  \\
        \alpha^2 c_2+(\alpha+\beta)^2c_3 & =0.
    \end{aligned}
    \right.
```

The solution is therefore given by

```math
 c_1=-\frac{2\alpha+\beta}{\alpha^2+\alpha\beta}\qquad c_2=\frac{1}{\alpha}+\frac{1}{\beta}\qquad c_3=-\frac{\alpha}{\beta(\alpha+\beta)}.
```

The same computation can be made for the backward and the centred FDM scheme, and lead to these three formulas:

```math
    \begin{aligned}
        Df(x_1) & = -\frac{2\alpha+\beta}{\alpha(\alpha+\beta)}f(x_1)+\left(\frac{1}{\alpha}+\frac{1}{\beta}\right)f(x_2)-\frac{\alpha}{\beta(\alpha+\beta)}f(x_3) \\
        Df(x_2) & = -\frac{\beta}{\alpha(\alpha+\beta)}f(x_1)+\left(\frac{1}{\alpha}-\frac{1}{\beta}\right)f(x_2)+\frac{\alpha}{\alpha(\alpha+\beta)}f(x_3)        \\
        Df(x_3) & = \frac{2\alpha+\beta}{\alpha(\alpha+\beta)}f(x_1)-\left(\frac{1}{\alpha}+\frac{1}{\beta}\right)f(x_2)+\frac{\alpha}{\beta(\alpha+\beta)}f(x_3).
    \end{aligned}
```

One can check than in the uniform case ($\alpha=\beta$) we recover the well known coefficient $-1/2$ and $1/2$ for the centred case
and $-3/2$, $2$ and $-1/2$ for the decentred case.

## Utility tools

The l\_norm\_tools.hpp file contains functions computing the infinity norm. For now, it computes the infinity norm of

- a double: $`\Vert x \Vert_{\infty} = x`$;
- a coordinate: $`\Vert x \Vert_{\infty} = \max_{i} (|x_i|)`$.
