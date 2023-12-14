# Tests on spline interpolator in polar coordinates

### Functions tests: 

* $(r,\theta) \mapsto r^p$ with $`0\leq p<4`$: interpolation with cubic bsplines supposed to be exact.
* $(r,\theta) \mapsto r^p$ with $p \geq 4$: interpolation with cubic bsplines supposed to be order 4.
* $(r,\theta) \mapsto r \cos(\theta)$ and $`(r,\theta) \mapsto \exp\left(- \frac{(x(r,\theta) - x_0)^2}{2\sigma_x^2} - \frac{(y(r,\theta) - y_0)^2}{2\sigma_y^2} \right)`$ : added and interpolation with cubic bsplines supposed to be order 4.

Test the interpolator on mesh points where the errors are supposed to be exact (machine error) and test the interpolator on "pseudo-random" points where the errors are supposed to be order 4 with a cubic bsplines basis. 
    
### Python scripts
    
* `test_convergence.py` : compute the convergence order in space for the tests where the errors are not supposed to be exact. 
    * launch the executable for a 64x64 grid and for a 128x128 grid to compare the errors. 
    * possibility of changing the default sizes by adding the two new sizes (N for NxN grid) at this end of the command (see `--help`).

    * Command to launch the test in this folder: `python3 test_convergence.py ../../../build/tests/geometryRTheta/2d_spline_interpolator/2d_spline_interpolator_tests` (`N0`) (`N1`)
    
* `display_curves.py` : display errors for several grid sizes and compute the slope of the errors for each function.

    * Command to launch the test in this folder: `python3 display_curves.py ./../../build/tests/geometryRTheta/2d_spline_interpolator/2d_spline_interpolator_tests`

 
