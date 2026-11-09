# SPDX-License-Identifier: MIT

add_executable(test_collSpVparMu test_collSpVparMu.cpp)
target_link_libraries(test_collSpVparMu
  PUBLIC
  DDC::core
  DDC::pdi
  PDI::PDI_C
  gyselalibxx::collisions
  gyselalibxx::geometry_vparmu
  gyselalibxx::geometry_collisions_vparmu
  gyselalibxx::initialisation_vparmu
  gyselalibxx::io
  gyselalibxx::quadrature
  gyselalibxx::speciesinfo
  gyselalibxx::utils

  paraconf::paraconf
  )

install(TARGETS test_collSpVparMu)

