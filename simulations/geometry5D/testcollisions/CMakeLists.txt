# SPDX-License-Identifier: MIT


add_executable(testcollisions testcollisions.cpp)
target_compile_features(testcollisions PUBLIC cxx_std_17)
target_link_libraries(testcollisions
  PUBLIC
  DDC::DDC
  DDC::PDI_Wrapper
  gslx::speciesinfo
  gslx::geometry5D
  gslx::paraconfpp
  gslx::quadrature
  gslx::collisions_rhs_5D
  PDI::pdi
  paraconf::paraconf
  )

install(TARGETS testcollisions)

