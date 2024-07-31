# SPDX-License-Identifier: MIT


add_executable(testcollisions testcollisions.cpp)
target_link_libraries(testcollisions
  PUBLIC
  DDC::DDC
  DDC::PDI_Wrapper
  gslx::geometry5D
  gslx::io
  gslx::paraconfpp
  gslx::quadrature
  gslx::collisions
  gslx::speciesinfo
  gslx::utils

  PDI::pdi
  paraconf::paraconf
  )

install(TARGETS testcollisions)

