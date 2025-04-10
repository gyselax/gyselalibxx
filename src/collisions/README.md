# Collisions

Collision operator in $`(v_\parallel,\mu)`$ applied to all against all species at the same time. As such, the operator is $O(N^2)$ in number of species.

It was rewritten from Fortran into C++. The operator's name is KOkkos coLIsion OPerator (KOLIOP) and is provided as a submodule (see `/vendor`).

While the operator is in C++, the interface with Fortran and other C++ code is made through a C layer. This is also known as the hourglass pattern (conceptually, the narrow neck of the device represent the lowering to C).

To integrate Koliop into Gyselalib++, we wrap its functionalities into a DDC aware operator present in `collision_operator.hpp`. In Gyselalib++, operator are expected to support multiple if not all kind of geometries. But Koliop expect some data in layout right [sp, phi, theta, r, vpar, mu] instead of the [sp, phi, r, theta, vpar, mu] layout that is going to be favoured in Gyselalib++. We have some machinery that setup input configuration data depending on the geometry. These are in `collision_configuration_sprvparmu.hpp`, `collision_configuration_spvparmu.hpp`.

More information can be found in the [Gysela collision operator](../../docs/latex/collisions/Gysela_collision.pdf)
