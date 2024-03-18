# Copyright (C) The DDC development team, see COPYRIGHT.md file
#
# SPDX-License-Identifier: MIT

add_executable(heat_equation heat_equation.cpp)
target_link_libraries(heat_equation PUBLIC DDC::DDC)
if("${DDC_BUILD_PDI_WRAPPER}")
  target_link_libraries(heat_equation PUBLIC DDC::PDI_Wrapper)
endif()
if("${DDC_BUILD_KERNELS_FFT}")
  add_executable(heat_equation_spectral heat_equation_spectral.cpp)
  target_link_libraries(heat_equation_spectral PUBLIC DDC::DDC)
endif()
add_executable(game_of_life game_of_life.cpp)
target_link_libraries(game_of_life PUBLIC DDC::DDC)
if("${DDC_BUILD_KERNELS_SPLINES}")
  add_executable(characteristics_advection characteristics_advection.cpp)
  target_link_libraries(characteristics_advection PUBLIC DDC::DDC)
endif()
