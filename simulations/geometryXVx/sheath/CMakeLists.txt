# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(sheath_${GEOMETRY_VARIANT} sheath.cpp)
target_compile_features(sheath_${GEOMETRY_VARIANT} PUBLIC cxx_std_17)
target_link_libraries(sheath_${GEOMETRY_VARIANT}
  PUBLIC
  DDC::DDC
  DDC::PDI_Wrapper
  vcx::initialization_${GEOMETRY_VARIANT}
  vcx::interpolation
  vcx::poisson_${GEOMETRY_VARIANT}
  vcx::time_integration_${GEOMETRY_VARIANT}
  vcx::speciesinfo
  vcx::boltzmann_${GEOMETRY_VARIANT}
  vcx::advection
  vcx::paraconfpp
  vcx::rhs_${GEOMETRY_VARIANT}
  PDI::pdi
  paraconf::paraconf
  sll::splines
  )

install(TARGETS sheath_${GEOMETRY_VARIANT})

endforeach()
