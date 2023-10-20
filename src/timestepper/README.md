# Time Stepping Methods

Time stepping methods are methods for calculating how a value (or dimensioned values) evolves over time. Such an evolution should be expressed as a system of equations of the form:

$$
dy/dt = f(t, y)\\
y(t_0) = y_0
$$

The implemented methods are:
- Second order Runge Kutta (RK2)
- Third order Runge Kutta (RK3)
- Fourth order Runge Kutta (RK4)

These classes all contain an `update` method which carries out one time step of the algorithm.
