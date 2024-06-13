# SPDX-License-Identifier: MIT


add_executable(testcollisions testcollisions.cpp)
target_compile_features(testcollisions PUBLIC cxx_std_17)
target_link_libraries(testcollisions
  PUBLIC
  DDC::DDC
  DDC::PDI_Wrapper
  gslx::geometry5D
  gslx::io
  gslx::paraconfpp
  gslx::quadrature
  gslx::collision
  gslx::speciesinfo
  PDI::pdi
  paraconf::paraconf
  )

install(TARGETS testcollisions)

