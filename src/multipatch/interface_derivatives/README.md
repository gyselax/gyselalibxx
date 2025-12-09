# Interface derivatives

On multi-patch geometries, we can impose a $`\mathcal{C}^1`$ regularity on the whole domain by fixing identical derivatives and function values on the interface.
The operators implemented in this folder propose methods to compute the derivatives at the patch interfaces.

The method applied here is detailed in Vidal et al. (2025)[^2].

**Assumptions:**

- the splines are cubic;
- :warning: every break point has to be an interpolation point for the splines.

:warning: **Warning:**
For the non-uniform meshes, please do not use Greville points as interpolation points but the
break points as interpolation points.
If you do not use derivatives as closure condition, two additional interpolation points have to be placed in
the boundary cells.

## Contents

- [Relation between derivatives on the boundaries of two connected patches](#relation-between-derivatives-on-the-boundaries-of-two-connected-patches): It documents the `SingleInterfaceDerivativesCalculator` operator.
  - [How to use the `SingleInterfaceDerivativesCalculator` operator?](#how-to-use-the-singleinterfacederivativescalculator-operator): It documents how to use the operator in the code.
  - [Formulae](#formulae): It details the formulae applies in the operator.
- [Relation between the interface derivatives along one direction](#relation-between-the-interface-derivatives-along-one-direction): It documents the `InterfaceExactDerivativeMatrix` operator.
  - [How to use the `InterfaceExactDerivativeMatrix` operator?](#how-to-use-the-InterfaceExactDerivativeMatrix-operator): It documents how to use the operator in the code.

## Relation between derivatives on the boundaries of two connected patches

We consider an interface between two patches.
The method computes the derivative at the interface of an equivalent global spline defined on the two patches of the interface.

Here is an illustration of the patch layout. For this example, we chose this patch connection but the method
and the implemented operator works for any patch connection.

![Illustration patch](../../../docs/images/interface_derivatives/fig1\_illustration\_patch.png "")

In the code, we are using the following relation between derivatives at consecutive interfaces,

```math
    s'(x_{i}) = c^i_{N^L,N^R} + a^i_{N^L,N^R} s'(x_{i+N^R}) + b^i_{N^L,N^R} s'(x_{i-N^L}).
```

with

- $`x_{i}`$ the interpolation point located on the interface,
- $`x_{i-N^L}`$ the interpolation point located on the interface to the left of the given interface (left bound on patch 1),
- $`x_{i+N^R}`$ the interpolation point located on the interface to the right of the given interface (right bound on patch 2),
- $`N^L`$ the number of cells on the left patch (patch 1),
- $`N^R`$ the number of cells on the right patch (patch 2),
- $`a^i_{N^L,N^R}`$, $`b^i_{N^L,N^R}`$  and $`c^i_{N^L,N^R}`$ coefficients to determine.

**Note:** $`a^i_{N^L,N^R}`$ and $`b^i_{N^L,N^R}`$ only depend on the interpolation points. The coefficient $`c^i_{N^L,N^R}`$
can be written as a linear combination of the function values $`\{f_{i+k}\}_{k = - N^L}^{N^R}`$ on the patches:

```math
    c^i_{N^L,N^R} = \sum_{k = - N^L}^{N^R} \omega_{k, N^L,N^R}^i f_{i+k},
```

with $`\{\omega_{k, N^L,N^R}^i\}_{k = - N^L}^{N^R}`$ depending only on the interpolation points.

**Approximation:** The formula gives the exact derivative of an equivalent global spline.
It is possible to use an approximation of this formula by using is smaller number of cells on each patch.
Instead of choosing the $`N^L`$ cells of patch 1, we can chose $`N^L_{reduc}`$ first cells from the interface on patch 1.
And we can do the same on patch 2 with $`N^R_{reduc}`$ first cells from the interface.
The approximation is then given by

```math
    s'(x_{i}) = \sum_{k = - N^L_{reduc}}^{N^R_{reduc}} \omega_{k, N^L_{reduc},N^R_{reduc}}^i f_{i+k}.
```

There are different ways to compute the coefficients $`a^i_{N^L,N^R}`$,
$`b^i_{N^L,N^R}`$, and $`c^i_{N^L,N^R}`$. There are described in the section
[Formulae](#formulae).

Firstly, we describe how to use the `SingleInterfaceDerivativesCalculator` operator.

### How to use the SingleInterfaceDerivativesCalculator operator?

When `SingleInterfaceDerivativesCalculator` is instantiated, it computes and stores the coefficients
$`a^i_{N^L,N^R}`$ and $`b^i_{N^L,N^R}`$ and the weights $`\{\omega_{k, N^L,N^R}^i\}_{k = - N^L}^{N^R}`$.
To instantiate it, we need the index ranges of the interpolation points. If we want to use the exact formula,
we need to provide the index ranges with all the points,

```cpp
SingleInterfaceDerivativesCalculator<Interface_12> derivatives_calculator (idx_range_patch_1, idx_range_patch_2);
```

If we want to use an approximation, we provide the index ranges with the interpolation points on the selected cells,

```cpp
Patch1::IdxRange1D idx_range_patch_1_reduced (idx_range_patch_1.take_first(N_L_reduc +1));
Patch2::IdxRange1D idx_range_patch_2_reduced (idx_range_patch_2.take_last(N_R_reduc +1));
SingleInterfaceDerivativesCalculator<Interface_12> derivatives_calculator (idx_range_patch_1_reduced, idx_range_patch_2_reduced);
```

> **Remark:** For interpolation with interpolation points as closure condition, a special treatment has to be carried out on the boundary cells
> (see [Additional interpolation point not on a break point](#additional-interpolation-point-not-on-a-break-point)).
> To apply this treatment, we specify it in the template parameter as follows.

```cpp
// If we want to apply the treatment on Patch 1 and Patch 2
SingleInterfaceDerivativesCalculator<Interface_12> 
    derivatives_calculator (idx_range_patch_1, idx_range_patch_2, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE);
// or if we want to apply the treatment only on Patch 1
SingleInterfaceDerivativesCalculator<Interface_12> 
    derivatives_calculator (idx_range_patch_1, idx_range_patch_2, ddc::BoundCond::GREVILLE, ddc::BoundCond::HERMITE);
// or if we want to apply the treatment only on Patch 2
SingleInterfaceDerivativesCalculator<Interface_12> 
    derivatives_calculator (idx_range_patch_1, idx_range_patch_2, ddc::BoundCond::HERMITE, ddc::BoundCond::GREVILLE);
```

> If we want to use an approximation where the boundary cells are not involved (even for interpolation points as closure condition on the global domain),
> the default template parameter can be used to avoid applying the treatment on inner cells.

The coefficients can be collected with the following functions:

- `derivatives_calculator.get_coeff_deriv_patch_1()` returns the coefficient $`b^i_{N^L,N^R}`$;
- `derivatives_calculator.get_coeff_deriv_patch_2()` returns the coefficient $`a^i_{N^L,N^R}`$;
- `derivatives_calculator.get_function_coefficients(function_1, function_2)` returns the coefficient $`c^i_{N^L,N^R}`$
 (or $`c^i_{N^L_{reduc},N^R_{reduc}}`$ for approximation).
- `derivatives_calculator.get_derivatives_coefficients(deriv_1, deriv_2)` is similar to `derivatives_calculator.get_function_coefficients(function_1, function_2)` but is recommended to compute the cross-derivatives in a 2D case. This class is also used for 2D cases where we need to compute the corner cross-derivatives. We apply the same method 
as for first derivatives but with first derivatives instead of function values. The difference in this operator is the change of sign of the first derivatives on patch 2 if the orientations of the patches of interface desagree. 

If we want to apply the exact formula, we need to sum these coefficients,

```cpp
double const coeff_deriv_left = derivatives_calculator.get_coeff_deriv_patch_1(); // coeff b
double const coeff_deriv_right = derivatives_calculator.get_coeff_deriv_patch_2(); // coeff a
double const sum_values = derivatives_calculator.get_function_coefficients(function_1, function_2); // coeff c

double const deriv_interface_right = 1e1; // a given value
double const deriv_interface_left = 1e1; // a given value

double const deriv_interface = sum_values + coeff_deriv_left * deriv_interface_left + coeff_deriv_right * deriv_interface_right;
```

If we want to apply an approximation of the formula, we only need $`c^i_{N^L_{reduc},N^R_{reduc}}`$,

```cpp
double const deriv_interface = derivatives_calculator.get_function_coefficients(function_1, function_2); // coeff c
```

**Remark:** It is also possible to use slices for the functions values,

```cpp
double const deriv_interface =
    derivatives_calculator.get_function_coefficients(
            function_1[idx_range_patch_1_reduced],
            function_2[idx_range_patch_2_reduced]); // coeff c
```

### Formulae

We use a method inspired by Crouseille et al. (2009)[^1], and detailed in Vidal et al.(2025)[^2].
To establish the formulae, we use the following relation between each consecutive interpolation point,

```math
    s'(x_i) = \gamma_i + \alpha_i  s'(x_{i+1}) + \beta_i  s'(x_{i-1})
```

with the *local coefficients* given by

```math
\left\{
\begin{aligned}
    & \gamma_i =
    \frac{3}{2} \frac{1}{\Delta x^R_{i} +\Delta x^L_{i}}
    \left[
        \frac{\Delta x^L_{i}}{\Delta x^R_{i}}  f^R_{i+1}
        + \left(
            \frac{\Delta x^R_{i}}{\Delta x^L_{i}}
            - \frac{\Delta x^L_{i}}{\Delta x^R_{i}}
        \right) f^R_{i}
        - \frac{\Delta x^R_{i}}{\Delta x^L_{i}}  f^L_{i-1}
    \right],
    \\
    & \alpha_i = -\frac{1}{2} \frac{\Delta x^L_{i}}{\Delta x^R_{i} +\Delta x^L_{i}}, \\
    & \beta_i =  -\frac{1}{2} \frac{\Delta x^R_{i}}{\Delta x^R_{i} +\Delta x^L_{i}}.
\end{aligned}
\right.
```

and $`\Delta x^R_{i} = x^R_{i+1} -  x^R_{i} \text{, } \Delta x^L_{i} = x^L_{i} -  x^L_{i-1}`$,
and $`\{f^{+/-}_{i+k}\}_k`$ the function values at the interpolation points on patch 1 (+) and
patch 2 (-).

**Remark:** As mentioned above, the coefficient $`c^i_{N^L,N^R}`$ can be written as a linear combination
of the interpolating function values. In the code, we store the weights in front of each interpolating function value
as they depend only on the grids. For the same reason, we store $`\gamma_i`$ as a vector,

```math
{\bf \Gamma_i}
    = \left[({\bf \Gamma_i})_0, ({\bf \Gamma_i})_1, ({\bf \Gamma_i})_2\right]
    = \frac{3}{2} \frac{1}{\Delta x^R_{i} +\Delta x^L_{i}}
        \left[
            - \frac{\Delta x^R_{i}}{\Delta x^L_{i}},
            \left(
                \frac{\Delta x^R_{i}}{\Delta x^L_{i}}
                - \frac{\Delta x^L_{i}}{\Delta x^R_{i}}
            \right),
            \frac{\Delta x^L_{i}}{\Delta x^R_{i}}
        \right].
```

#### Recursive formula

We start the recursion on patch 2. We initialise it with,

```math
\begin{aligned}
    & c^i_{1,1} = \gamma_i,
    \qquad
    && a^i_{1,1} = \alpha_i,
    \qquad
    && b^i_{1,1} = \beta_i,
    \\
    & c^i_{1,2} = \frac{1}{1 - \alpha_i  \beta_{i+1}} \left[ \gamma_i + \alpha_i \gamma_{i+1} \right],
    \qquad
    && a^i_{1,2} = \frac{\alpha_i \alpha_{i+1}}{1 - \alpha_i  \beta_{i+1}},
    \qquad
    && b^i_{1,2} = \frac{\beta_i}{1 - \alpha_i  \beta_{i+1}},
\end{aligned}
```

(i.e. for the weights,

```math
\begin{aligned}
    & \omega_{k, 1,1}^i = \gamma_{i,k},
    && \omega_{k, 1,2}^i = \frac{1}{1 - \alpha_i  \beta_{i+1}} \omega_{k, 1,1}^i
                        + \frac{\alpha_i}{1 - \alpha_i  \beta_{i+1}} \gamma_{i+1,k},
\end{aligned}
```

with $`\gamma_{i,i-1} = ({\bf \Gamma_i})_0, \ \gamma_{i,i} = ({\bf \Gamma_i})_1, \ \gamma_{i,i+1} = ({\bf \Gamma_i})_2`$
and null for different *k*.)

and for $`n = 2, ..., N^R-1`$,

```math
\left\{
\begin{aligned}
    & c^i_{1,n+1} = \frac{1}{1 -  \beta_{i+n} \frac{a^i_{1,n}}{a^i_{1,n-1}} }
        \left[
            c^i_{1,n}
            + a^i_{1,n} \gamma_{i+n}
            - \beta_{i+n} \frac{a^i_{1,n}}{a^i_{1,n-1}} c^i_{1,n-1}
        \right] \\
    & a^i_{1,n+1} = \frac{a^i_{1,n} \alpha_{i+n}}{1 -  \beta_{i+n} \frac{a^i_{1,n}}{a^i_{1,n-1}} }, \\
    & b^i_{1,n+1} = \frac{1}{1 -  \beta_{i+n} \frac{a^i_{1,n}}{a^i_{1,n-1}} }
        \left[
            b^i_{1,n}
            - \beta_{i+n} \frac{a^i_{1,n} }{a^i_{1,n-1}} b^i_{1,n-1}
        \right].
\end{aligned}
\right.
```

(i.e. for the weights,

```math
\begin{aligned}
    & \omega_{k, 1,n+1}^i = \frac{1}{1 -  \beta_{i+n} \frac{a^i_{1,n}}{a^i_{1,n-1}} }
        \left[
            \omega_{k, 1,n}^i
            - \beta_{i+n} \frac{a^i_{1,n}}{a^i_{1,n-1}} \omega_{k, 1,n-1}^i
        \right]
        + \frac{1}{1 -  \beta_{i+n} \frac{a^i_{1,n}}{a^i_{1,n-1}} }a^i_{1,n} \gamma_{i+n, k}. \\
\end{aligned}
```

)

It gives us a relation linking the derivatives at the three points on the following illustration,

![Illustration recursion 1](../../../docs/images/interface_derivatives/fig2\_illustration\_recursion\_1.png "")

We keep doing the recursion on patch 1,

```math
    c^i_{2,N^R} = \frac{1}{1 - b^{i}_{1,N^R}  \alpha_{i-1}}
        \left[
            c^{i}_{1,N^R} + b^{i}_{1,N^R}  \gamma_{i-1}
        \right],
    \qquad
    a^i_{2,N^R} = \frac{a^{i}_{1,N^R}}{1 - b^{i}_{1,N^R} \alpha_{i-1}},
    \qquad
    b^i_{2,N^R} = \frac{b^{i}_{1,N^R}   \beta_{i-1}}{1 - b^{i}_{1,N^R}  \alpha_{i-1}},
```

(i.e. for the weights,

```math
\begin{aligned}
    & \omega_{k, 2,N^R}^i = \frac{1}{1 - b^{i}_{1,N^R}  \alpha_{i-1}} \omega_{k, 1,N^R}^i
                            + \frac{1}{1 - b^{i}_{1,N^R}  \alpha_{i-1}} b^{i}_{1,N^R}  \gamma_{i-1, k}.
\end{aligned}
```

)

and for $`m = 2, ..., N^L-1`$,

```math
\left\{
\begin{aligned}
        & c^{i}_{m+1,N^R} = \frac{1}{1 - \alpha_{i-m} \frac{b^{i}_{m,N^R}}{b^{i}_{m-1,N^R}}}
            \left[
                c^{i}_{m,N^R}
                + b^{i}_{m,N^R}  \gamma_{i-m}
                - \alpha_{i-m} \frac{b^{i}_{m,N^R}}{b^{i}_{m-1,N^R}} c^{i}_{m-1,N^R}
            \right],
        \\
        & a^{i}_{m+1,N^R} = \frac{1}{1 - \alpha_{i-m} \frac{b^{i}_{m,N^R}}{b^{i}_{m-1,N^R}}}
            \left[
                a^{i}_{m,N^R}
                - \alpha_{i-m} \frac{b^{i}_{m,N^R}}{b^{i}_{m-1,N^R}} a^{i}_{m-1,N^R}
            \right],
        \\
        & b^{i}_{m+1,N^R} = \frac{b^{i}_{m,N^R} \beta_{i-m} }
            {1 - \alpha_{i-m} \frac{b^{i}_{m,N^R}}{b^{i}_{m-1,N^R}}}.
\end{aligned}
\right.
```

(i.e. for the weights,

```math
\begin{aligned}
    & \omega_{k, m+1,N^R}^i =
    \frac{1}{1 - \alpha_{i-m} \frac{b^{i}_{m,N^R}}{b^{i}_{m-1,N^R}}}
            \left[
                \omega_{k, m,N^R}^i
                - \alpha_{i-m} \frac{b^{i}_{m,N^R}}{b^{i}_{m-1,N^R}} \omega_{k, m-1,N^R}^i
            \right]
    \frac{1}{1 - \alpha_{i-m} \frac{b^{i}_{m,N^R}}{b^{i}_{m-1,N^R}}} b^{i}_{m,N^R}  \gamma_{i-m, k}.
    \\
\end{aligned}
```

)

The two recursions gives us a relation linking the derivatives at the three consecutive interface
(see the following illustration),

![Illustration recursion 2](../../../docs/images/interface_derivatives/fig3\_illustration\_recursion\_2.png "")

#### Additional interpolation point not on a break point

If we use interpolation points as closure conditions for the equivalent global spline,
the last local coefficients must be modified.

:warning: **Warning:**
We consider here that all the break points are interpolation points. For cubic splines,
two additional interpolation points are added.
They are added in the boundary cells of the mesh where the equivalent global spline is defined.

![Illustration additional interpolation points.](../../../docs/images/interface_derivatives/fig4\_illustration\_additional\_interpolation\_pts.png "")

If there is an additional interpolation point in the left boundary cell, then the last step
of the second part of the recursion need to use these modified local coefficients,

```math
\left\{
\begin{aligned}
    & \gamma_i^{*,-} = \frac{1}{1 + \beta_i \frac{K_1^{*,-}}{K_0^{*,-}}}
        \left[
            \gamma_i
            + \beta_i \frac{1}{\Delta x^L_{i} K_0^{*,-}} s(x^L_{*})
            - \beta_i \frac{H_0^{*,-}}{\Delta x^L_{i} K_0^{*,-}} s(x^L_{i-1})
            - \beta_i \frac{H_1^{*,-}}{\Delta x^L_{i} K_0^{*,-}} s(x^L_{i})
        \right],
    \\
    & \alpha_i^{*,-} = \frac{\alpha_i}{1 + \beta_i \frac{K_1^{*,-}}{K_0^{*,-}}}, \\
    & \beta_i^{*,-} = 0
\end{aligned}
\right.
```

with

```math
    H_0(x) = (1-x)^2 (1+2x),
    \qquad
    H_1(x) = x^2 (3-2x),
    \qquad
    K_0(x) = (1-x)^2 x,
    \qquad
    K_1(x) = x^2 (x-1),
```

with $`H_0^{*,-} \text{, } H_1^{*,-} \text{, } K_0^{*,-} \text{ and } K_1^{*,-}`$ the evaluations of
$`H_0 \text{, } H_1 \text{, } K_0 \text{ and } K_1`$ at $`\frac{x^L_{*} - x^L_{i-1}}{\Delta x^L_{i}}`$,
with $`x^L_{*}`$ the additional interpolation point.

If there is an additional interpolation point in the right boundary cell, then the last step
of the first part of the recursion need to use these modified local coefficients,

```math
\left\{
\begin{aligned}
    & \gamma_i^{*,+} = \frac{1}{1 + \alpha_i \frac{K_0^{*,+}}{K_1^{*,+}}}
        \left[
            \gamma_i
            + \alpha_i \frac{1}{\Delta x^R_{i} K_1^{*,+}} s(x^R_{*})
            - \alpha_i \frac{H_0^{*,+}}{\Delta x^R_{i} K_1^{*,+}} s(x^R_{i})
            - \alpha_i \frac{H_1^{*,+}}{\Delta x^R_{i} K_1^{*,+}} s(x^R_{i+1})
        \right],
    \\
    & \alpha_i^{*,+} = 0, \\
    & \beta_i^{*,+} = \frac{\beta_i}{1 + \alpha_i \frac{K_0^{*,+}}{K_1^{*,+}}}.
\end{aligned}
\right.
```

and $`H_0^{*,+} \text{, } H_1^{*,+} \text{, } K_0^{*,+} \text{ and } K_1^{*,+}`$ the evaluations of
$`H_0 \text{, } H_1 \text{, } K_0 \text{ and } K_1`$ at $`\frac{x^R_{*} - x^R_{i}}{\Delta x^R_{i}}`$,
with $`x^R_{*}`$ the additional interpolation point.

#### Explicit formula

If the **interpolation points are uniform**
($`\Delta x^L_{i} = \Delta x^L  \text{ and } \Delta x^R_{i} = \Delta x^R \text{, } \forall i`$),
then the recursion formula is equivalent to the following explicit formula.

To lighten the notation, we note $`n_1 = N^L`$ and $`n_2 = N^R`$.

We introduce the following sequence to define the explicit formulae,

```math
    u_{k} = (2+\sqrt{3})^k - (2-\sqrt{3})^k,
    \qquad \forall k \in \mathbb{N}.
```

```math
\left\{
\begin{aligned}
    & a^I_{n_1,n_2} = \frac{(-1)^{n_2-1}u_{1}a^I_{1,1}u_{n_1}}{u_{n_1}u_{n_2} + u_{n_1}u_{n_2-1}a^I_{1,1} + u_{n_2} u_{n_1-1}b^I_{1,1}},
    \\
    & b^I_{n_1,n_2} = \frac{(-1)^{n_1-1} u_{1} b^I_{1,1} u_{n_2}}{u_{n_1}u_{n_2} + u_{n_1}u_{n_2-1}a^I_{1,1} + u_{n_2} u_{n_1-1}b^I_{1,1}},
    \\
    & c^I_{n_1,n_2} = \sum_{k=-n_1}^{n_2} \omega^I_{k,n_1,n_2} f_{i+k},
\end{aligned}
\right.
```

with the weights given by,

```math
\left\{
\begin{aligned}
    & \omega^I_{k,n_1,n_2} = 3(-1)^{k}
        \frac{\frac{a^I_{1,1}}{\Delta x^R} u_{n_1}u_{1}}
        {u_{n_1}u_{n_2} + u_{n_1}u_{n_2-1}a^I_{1,1} + u_{n_2} u_{n_1-1}b^I_{1,1}}
        , &k = n_2, \\
    & \omega^I_{k,n_1,n_2} = 3(-1)^{k}
        \frac{\frac{a^I_{1,1}}{\Delta x^R} u_{n_1}(u_{n_2-k+1} - u_{n_2-k-1})}
        {u_{n_1}u_{n_2} + u_{n_1}u_{n_2-1}a^I_{1,1} + u_{n_2} u_{n_1-1}b^I_{1,1}}
        , &k = 1, ..., n_2-1, \\
    & \omega^I_{k,n_1,n_2} = 3(-1)^k
        \frac{\frac{a^I_{1,1}}{\Delta x^R}u_{n_1} (u_{n_2} - u_{n_2-1}) - \frac{b^I_{1,1}}{\Delta x^L}u_{n_2}(u_{n_1} - u_{n_1-1})}
        {u_{n_1}u_{n_2} + u_{n_1}u_{n_2-1}a_{1,1} + u_{n_2} u_{n_1-1}b^I_{1,1}}
        , &k = 0, \\
    & \omega^I_{k,n_1,n_2} = 3(-1)^{k+1}
        \frac{\frac{b^I_{1,1}}{\Delta x^L} u_{n_2}(u_{n_1+k+1} - u_{n_1+k-1})}
        {u_{n_1}u_{n_2} + u_{n_1}u_{n_2-1}a^I_{1,1} + u_{n_2} u_{n_1-1}b^I_{1,1}}
        , &k = -(n_1-1), ..., -1, \\
    & \omega^I_{k,n_1,n_2} = 3(-1)^{k+1}
        \frac{\frac{b^I_{1,1}}{\Delta x^L} u_{n_2}u_{1}}
        {u_{n_1}u_{n_2} + u_{n_1}u_{n_2-1}a^I_{1,1} + u_{n_2} u_{n_1-1}b^I_{1,1}}
        , &k = -n_1,
\end{aligned}
\right.
```

and $`a^I_{1,1} = -\frac{1}{2} \frac{\Delta x^L}{\Delta x^R +\Delta x^L}`$
and $`b^I_{1,1} =  -\frac{1}{2} \frac{\Delta x^R}{\Delta x^R +\Delta x^L}`$.

## Relation between the interface derivatives along one direction

We consider now a set of patches connected via different interfaces.
We want to compute the derivatives at all the interfaces.

To compute all the interface derivatives along one direction, we collect all the relations between
three consecutive derivatives for each interface.
(For a given interface, the relation is defined in the previous section,
see [Relation between derivatives on the boundaries of two connected patches](#relation-between-derivatives-on-the-boundaries-of-two-connected-patches)).

All these relations can be stored in a matrix system as follows,

```math
\begin{bmatrix}
    s'(\mathcal{X}_1) \\
    \vdots \\
    s'(\mathcal{X}_I) \\
    \vdots \\
    s'(\mathcal{X}_{N_I-1}) \\
\end{bmatrix}
=
\begin{bmatrix}
    1 & -a^1 & \dots & 0\\
    -b^2 & 1 & -a^2 & \vdots \\
    \vdots & \ddots & \ddots & -a^{N_I -2}\\
    0 & \dots & -b^{N_I -1}& 1 \\
\end{bmatrix}^{-1}
\begin{bmatrix}
    c^1 + b^1 s'(\mathcal{X}_0) \\
    \vdots \\
    c^I \\
    \vdots \\
    c^{N_I-1} + a^{N_I-1} s'(\mathcal{X}_{N_I}) \\
\end{bmatrix},
```

with

- $`a^I, b^I \text{ and } c^I`$ the coefficients computed with `SingleInterfaceDerivativeCalculator`
for a given interface *I* (and given number of cells on the right and the left patches that we do not
precise here to lighten the notation).
- $`\{\mathcal{X}_I\}_I`$ a set of parallel interfaces.
- $`s'(\mathcal{X})`$ the derivative at the *I*th interface.

We simply refer it as

```math
    S = (\mathbb{I}-M)^{-1}C.
```

> **Remark:** This matrix is given in the case where the equivalent global spline passing through all
> the patches in the given direction use **Hermite boundary conditions**.
>
> If the equivalent global spline uses **additional interpolation points as closure condition**
> (see [Additional interpolation point not on a break point](#additional-interpolation-point-not-on-a-break-point)), then $`b^1 = 0 \text{ and } a^{N_I-1} = 0`$.
> So, the dependency on the boundary derivatives disappears in the vector *C*.
>
> If the equivalent global spline uses **periodic boundary conditions**, then we can add an additional
> relation in the matrix system to close the problem,

```math
\begin{bmatrix}
    s'(\mathcal{X}_0) \\
    s'(\mathcal{X}_1) \\
    \vdots \\
    s'(\mathcal{X}_I) \\
    \vdots \\
    s'(\mathcal{X}_{N_I-1}) \\
\end{bmatrix}
=
\begin{bmatrix}
    1 & -a^0 & 0 &\dots & -b^0\\
    -b^1 & 1 & -a^1 &  & 0\\
    0 & -b^2 & 1 & -a^2 & \vdots \\
    \vdots & & & & \\
    & & & \ddots & \\
    -a^{N_I -1} & 0 & \dots & -b^{N_I -1}& 1 \\
\end{bmatrix}^{-1}
\begin{bmatrix}
    c^0 \\
    c^1 \\
    \vdots \\
    c^I \\
    \vdots \\
    c^{N_I-1} \\
\end{bmatrix}.
```

### How to use the InterfaceExactDerivativeMatrix operator?

First, for each interface in the geometry, we instantiate a `SingleInterfaceDerivatorCalculator`
(see [How to use the SingleInterfaceDerivatorCalculator operator?](#how-to-use-the-singleinterfacederivativescalculator-operator)).
We store a constant reference of all the derivative calculator in a `SingleInterfaceDerivatorCalculatorCollection`.

```cpp
SingleInterfaceDerivatorCalculator<Interface_1> derivative_calculator_1(...);
SingleInterfaceDerivatorCalculator<Interface_2> derivative_calculator_2(...);
...

SingleInterfaceDerivatorCalculatorCollection derivative_calculators (derivative_calculator_1, derivative_calculator_2, ...);
```

We can then instantiate `InterfaceExactDerivativeMatrix` with the tuple of derivative calculator.

```cpp
InterfaceExactDerivativeMatrix<
        Connectivity,                               // MultipatchConnectivity class
        Grid1D,                                     // the given direction.
        ddc::detail::TypeSeq<Patch1, Patch2, ...>   // list of patches containing all the needed ones
        BoundCondGlobalLower,                       // boundary condition for the equivalent global spline
        BoundCondGlobalUpper,                       // boundary condition for the equivalent global spline
        SingleInterfaceDerivatorCalculatorCollection<Interface_1, Interface_2, ...>>
        matrix(idx_ranges, derivative_calculators);
```

with `idx_ranges` a `MultipatchType<IdxRangeonPatch, Patch1, Patch2, ...>` object.

During the instantiation, `InterfaceExactDerivativeMatrix` will allocate memory for the matrix $`(\mathbb{I} - M)`$,
for the right hand side vector *C* and for the solution vector *S*.
It also computes the matrix and stores it.

If we want to solve the system for a given function, we then use the operator `.solve_deriv()`.

```cpp
matrix.solve_deriv(
    functions_and_derivs,         // A MultipatchField collection of DerivField.
    );
```

During the call to `.solve_deriv()`, it computes the right hand side vector *C*.
It solves the matrix system to get the solution vector *S*.
It fill in the derivatives in th `functions_and_derivs` collection with the computed derivatives at the right place.

:warning: The derivatives in `functions_and_derivs` are overwritten during the last step.

#### A 2D case operator

This operator is actually implemented for 2D patches.
For 2D local splines with Hermite boundary conditions, we need the first derivatives along the first dimension,
first derivatives along the second dimension, and the cross-derivatives.

Let's take the following case,

![Illustration example](../../../docs/images/interface_derivatives/fig5\_example\_9\_patches.png "")

This geometry is composed of 9 patches forming periodic strips in the *x* direction.
We use additional interpolation points as closure condition for the equivalent global splines in
the *y* direction (i.e. `ddc::BoundCond::GREVILLE`).

So, we start by defining the `InterfaceExactDerivativeMatrix` matrices for each periodic directions
$`\vec{x_1}, \vec{x_4}, \text{ and } \vec{x_7}`$ (`GridX1`, `GridX4` and `GridX7`).

```cpp
InterfaceExactDerivativeMatrix<Connectivity, GridX1, 
        ddc::detail::TypeSeq<Patch1, Patch2, Patch3>, 
        ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC,
        SingleInterfaceDerivatorCalculatorCollection<Interface_1_2, Interface_2_3, Interface_3_1>>
        matrix_123(idx_ranges_123, derivative_calculators_123);

InterfaceExactDerivativeMatrix<Connectivity, GridX4, 
        ddc::detail::TypeSeq<Patch4, Patch5, Patch6>,
        ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC,
        SingleInterfaceDerivatorCalculatorCollection<Interface_4_5, Interface_5_6, Interface_6_4>>
        matrix_456(idx_ranges_456, derivative_calculators_456);
// ...
```

We also define `InterfaceExactDerivativeMatrix` matrices for each non-periodic directions
$`\vec{y_1}, \vec{y_2}, \text{ and } \vec{y_3}`$ (`GridY1`, `GridY2` and `GridY3`).

```cpp
InterfaceExactDerivativeMatrix<Connectivity, GridY1, 
        ddc::detail::TypeSeq<Patch1, Patch4, Patch7>,
        ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE,
        SingleInterfaceDerivatorCalculatorCollection<Interface_1_4, Interface_4_7>>
        matrix_147(idx_ranges_147, derivative_calculators_147);

InterfaceExactDerivativeMatrix<Connectivity, GridY2, 
        ddc::detail::TypeSeq<Patch2, Patch5, Patch8>,
        ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE,
        SingleInterfaceDerivatorCalculatorCollection<Interface_2_5, Interface_5_8>>
        matrix_258(idx_ranges_258, derivative_calculators_258);
// ...
```

> **Note:** the list of patches need to have all the patches on the given direction.
> In case of doubt, the full list of patches can be given. In this case, all the given collections
> have to be defined on the same patch set.

With these matrices, we can compute all the interface *x*-derivatives using the directions
$`\vec{y_1}, \vec{y_2}, \text{ and } \vec{y_3}`$,

```cpp
matrix_147.solve_deriv(functions_and_derivs_147);
matrix_258.solve_deriv(functions_and_derivs_258);
// ...
```

we can compute all the interface *y*-derivatives using the directions
$`\vec{x_1}, \vec{x_4}, \text{ and } \vec{x_7}`$,

```cpp
matrix_123.solve_deriv(functions_and_derivs_123);
matrix_456.solve_deriv(functions_and_derivs_456);
// ...
```

From the computed first derivatives, we can compute the cross-derivatives along one of the other direction.
E.g.

```cpp
matrix_123.solve_cross_deriv(functions_and_derivs_123); // along x direction using y-derivatives 
// ...
// or
matrix_147.solve_cross_deriv(functions_and_derivs_147); // along y direction using x-derivatives
// ...
```

Once the all derivatives computed on for every patches using all the interfaces,
all have the data to build local spline representations.
On conforming global meshes, the local splines are exactly pieces of an equivalent
global spline.


:warning: **Warning:** 
The cross-derivatives are computed from the first derivatives. So, the `.solve_deriv()` has to be applied before the `.solve_cross_deriv()`. 
It is recommended to applied `.solve_deriv()` on _every_ group of patches before. 

## References

[^1]: Crouseilles, N., Latu, G., Sonnendrücker, E.:
*A parallel vlasov solver based on local cubic spline interpolation on patches.*
Journal of Computational Physics 228(5), 1429–1446 (2009)

[^2]: Vidal, P., Bourne, E., Grandgirard, V., Mehrenberger, M., Sonnendrücker, E.,
*Local cubic spline interpolation for Vlasov-type equations on a multi-patch geometry.*
Journal of Scientific Computing, (2025) [ACCEPTED].
Available on arXiv: [https://arxiv.org/abs/2505.22078](https://arxiv.org/abs/2505.22078)
