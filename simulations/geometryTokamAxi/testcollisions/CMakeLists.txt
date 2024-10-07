# SPDX-License-Identifier: MIT

add_executable(test_collTokamAxi test_collTokamAxi.cpp)
target_link_libraries(test_collTokamAxi PUBLIC
  DDC::DDC
  DDC::PDI_Wrapper
  gslx::geometry_tokamaxi
  gslx::io
  gslx::paraconfpp
  gslx::quadrature
  gslx::speciesinfo
  gslx::utils

  PDI::pdi
  paraconf::paraconf
  )

install(TARGETS test_collTokamAxi)

