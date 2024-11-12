# SPDX-License-Identifier: MIT

add_executable(test_collTokamAxi test_collTokamAxi.cpp)
target_link_libraries(test_collTokamAxi PUBLIC
  DDC::DDC
  DDC::PDI_Wrapper
  gslx::advection
  gslx::collisions
  gslx::geometry_tokamaxi
  gslx::interpolation
  gslx::io
  gslx::mpi_parallelisation
  gslx::paraconfpp
  gslx::quadrature
  gslx::speciesinfo
  gslx::timestepper
  gslx::utils

  PDI::pdi
  paraconf::paraconf
  )

install(TARGETS test_collTokamAxi)

