# Geometry5D : 

The geometry folder contains a helper file `geometry.hpp` which provides shortcuts to the types needed to define the 5D geometry (3D in space + 2D in velocity).
- Concerning the 3D space, several coordinate systems are used:
  * (Cart1, Cart2, Cart3) : 3D cartesian geometry (x, y, z)
  * (Tor1, Tor2, Tor3) : 3D toric geometry (r, theta, phi)
  * (Cyl1, Cyl2, Cyl3) : 3D cylindric geometry (R, Z, phi)
- Concerning the 2D space, we consider:
  * Vpar : the parallel velocity 
  * Mu : the magnetic momentum (i.e Mu = Vperp^2/(2B))
  * V2D = (Vpar, Mu)




