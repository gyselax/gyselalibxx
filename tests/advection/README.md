# Tests on the templated advection operators

## Contents 

### Advection along x on (t, x)
* 1d\_advection\_x.cpp: test the BslAdvection1D operator for a 1D case. 

The test is the following: 
```math
    \partial_t f(t,x) + A(x) \partial_x f(t,x) = 0,
```

on $`\Omega = [-\pi, \pi]`$ with 

```math
\left\{
\begin{aligned}
    & f_0(x) = \sin(2x), \\
    & A(x) = \sin(x).
\end{aligned}
\right.
```

We test it on a grid $`N_x = 32`$  with $`dt = 0.05`$. The simulation runs on $`t\in[0,0.4]`$ and the relative error is expected to be below $` 5*10^{-3}`$. 


### Advection along x on (t, x, vx)
* 1d\_advection\_xvx.cpp: test the BslAdvection1D operator for a 1Dx1V case. 

The test is the following: 
```math
    \partial_t f(t,x, v_x) + v_x \partial_x f(t,x, v_x) = 0, \\
```

on $`\Omega = [-\pi, \pi]\times[-6, 6]`$ with 

```math
\left\{
\begin{aligned}
    & f_0(x, v_x) = \cos(x), \\
    & A(x, v_x) = v_x.
\end{aligned}
\right.
```

We test it on a grid $`N_x\times N_{v_x} = 100\times50`$  with $`dt = 0.1`$. The simulation runs on $`t\in[0,0.4]`$ and the relative error is expected to be below $` 5*10^{-7}`$. 


### Advection along x and y on (t, x, y, vx, vy)
* 1d\_advection\_xyvxvy.cpp: test the BslAdvection1D operator for a 2Dx2V case with a Strang splitting for the resolution. 

The test is the following: 
```math
    \partial_t f(t,x,y,v_x, v_y) + A(x,y) \cdot \nabla f(t,x,y, v_x, v_y) = 0,
```

on $`\Omega = [-0.5, 0.5]^2\times\{0,1\}^2`$ with 

```math
\left\{
\begin{aligned}
    & f_0(x) = \frac{1}{2}(G(r_1(x,y)) + G(r_2(x,y))), \\
    & A(x, y) = \frac{1}{2\pi}
    \begin{bmatrix}
        \sin(2\pi x) \\
        \sin(2\pi y) \\
    \end{bmatrix},
\end{aligned}
\right.
```

with 

```math
\left\{
\begin{aligned}
    & G (r)  = \cos(\frac{\pi r}{2a}) * \mathbb{1}_{r<a}(r), \\
    & r_1 (x,y) = \sqrt{x^2 + 8 y^2},  \\
    & r_2 (x,y) = \sqrt{8 x^2 + y^2},  \\
    & \omega = 2 \pi, \\
    & (x_c, y_c) = (0.25, 0). \\
\end{aligned}
\right.
```

We test it on a grid $`N_x \times N_y \times N_{v_x} \times N_{v_y} = 60 \times 60 \times 2 \times 2`$  with $`dt = 0.05`$. The simulation runs on $`t\in[0,0.2]`$ and the relative error is expected to be below $` 7*10^{-2}`$. 

### Same test cases as for BslAdvectionSpatial and BslAdvectionVelocity
* 1d\_spatial\_advection.cpp: test the BslAdvection1D operator for a 1Dx1V case.  
This test is similar to `tests/geometryXVx/spatialadvection.cpp` which tests the BslAdvectionSpatial operator. 

The test is the following: 
```math
    \partial_t f_s(t,x, v_x) - \sqrt{\frac{m_e}{m_s}} v_x \partial_x f_s(t,x, v_x) = 0, \\
```

on $`\Omega = [-\pi, \pi]\times[-6, 6]`$ with 

```math
\left\{
\begin{aligned}
    & f_0(x, v_x) = \cos(x), \\
    & A(x, v_x) = - \sqrt{\frac{m_e}{m_s}} v_x.
\end{aligned}
\right.
```

We test it on a grid $`N_x\times N_{v_x} = 100\times50`$  with $`dt = 0.1`$. The simulation runs on $`t\in[0,0.1]`$ and the relative error is expected to be below $` 1*10^{-6}`$. 


* 1d\_velocity\_advection.cpp: test the BslAdvection1D operator for a 1Dx1V case.  
This test is similar to `tests/geometryXVx/velocityadvection.cpp` which tests the BslAdvectionVelocity operator. 

The test is the following: 
```math
    \partial_t f_s(t,x, v_x) - q_s\sqrt{\frac{m_e}{m_s}} E \partial_{v_x} f_s(t,x, v_x) = 0, \\
```

on $`\Omega = [0, 2\pi]\times[-10, 10]`$ with 

```math
\left\{
\begin{aligned}
    & f_0(x, v_x) = \exp\left(-\frac{1}{2}v^2\right), \\
    & A(x, v_x) =  - q_s\sqrt{\frac{m_e}{m_s}} E.
\end{aligned}
\right.
```

We test it on a grid $`N_x\times N_{v_x} = 50\times100`$  with $`dt = 0.1`$. The simulation runs on $`t\in[0,0.1]`$ and the relative error is expected to be below $` 1*10^{-5}`$. 
