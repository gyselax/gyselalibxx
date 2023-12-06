# Advection methods

The `advection/` folder gathers the backward semi lagrangian scheme classes. There are two main classes, AdvectionSpatial and AdvectionVelocity. Implementing these operators separately makes sense, since we are using time splitting.

Two implementations of the operators are available. For both cases the feet of the characteristic curves are used to interpolate the updated distribution function on mesh points. The sequential version uses classical interpolators and needs to iterate over a slice of interest in the phase space. The parallel one has the suffix "batched". It uses batched interpolators so the interpolation step is done over the whole distribution function.

## Spatial advection
Here the purpose is the advection along a direction on the physical space dimension of the phase space. 
The dynamics of the motion on the spatial dimension are governed by the following equation. 

$$ \frac{df_s}{dt}= \sqrt{\frac{m_e}{m_s}} v \frac{\partial f_s}{\partial x} $$

## Velocity advection
Here the purpose is the advection along a direction on the velocity space dimension of the phase space.
The dynamics of the motion on the velocity dimension are governed by the following equation, where E is the electric field.

$$ \frac{df_s}{dt}= q_s \sqrt{\frac{m_e}{m_s}} E \frac{\partial f_s}{\partial v} $$
