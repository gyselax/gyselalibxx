# SPDX-License-Identifier: MIT

add_executable(test_collSpVparMu test_collSpVparMu.cpp)
target_link_libraries(test_collSpVparMu
  PUBLIC
  DDC::core
  DDC::pdi
  PDI::pdi
  gslx::collisions
  gslx::geometry_vparmu
  gslx::initialization_vparmu
  gslx::io
  gslx::quadrature
  gslx::speciesinfo
  gslx::utils

  paraconf::paraconf
  )

install(TARGETS test_collSpVparMu)

