# Interface derivatives

This folder contains functions to estimate derivatives at the boundary between 2 patches. The patches may represent different geometries across a local multi-patch object, or similar geometries distributed across MPI ranks.

The following methods are implemented so far:
- DerivReconstructionCoeffs [^2]

[^2]: Vidal, P., Bourne, E., Grandgirard, V., Mehrenberger, M., Sonnendrücker, E.,
*Local cubic spline interpolation for Vlasov-type equations on a multi-patch geometry.*
Journal of Scientific Computing 106, 58 (2026). <https://doi.org/10.1007/s10915-025-03161-z>

