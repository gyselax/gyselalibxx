# Time integration

The time integration folder contains methods for solving a Vlasov-Poisson system of equations defined as:
```math
\partial_t f(t, x, v) + v_x \partial_x f + v_y \partial_y f + \frac{E_x}{m} \partial_{v_x} f + \frac{E_y}{m} \partial_{v_x} f = 0
-\Delta \phi(x) = \sum_s \rho_{q,s}(x,y)
\rho_{q,s}(x, y) = \int q_s f_s(x,y,v_x,v_y) dv_x v_y
E = - \nabla \phi
```

The implemented time integrators are:
- PredCorr : A predictor-corrector method

