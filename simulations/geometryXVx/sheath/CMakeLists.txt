# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(sheath_${GEOMETRY_VARIANT} sheath.cpp)
target_compile_features(sheath_${GEOMETRY_VARIANT} PUBLIC cxx_std_17)
target_link_libraries(sheath_${GEOMETRY_VARIANT}
  PUBLIC
  DDC::DDC
  DDC::PDI_Wrapper
  gslx::initialization_${GEOMETRY_VARIANT}
  gslx::interpolation
  gslx::poisson_${GEOMETRY_VARIANT}
  gslx::time_integration_${GEOMETRY_VARIANT}
  gslx::speciesinfo
  gslx::boltzmann_${GEOMETRY_VARIANT}
  gslx::advection
  gslx::paraconfpp
  gslx::rhs_${GEOMETRY_VARIANT}
  PDI::pdi
  paraconf::paraconf
  sll::splines
  )

install(TARGETS sheath_${GEOMETRY_VARIANT})

endforeach()
