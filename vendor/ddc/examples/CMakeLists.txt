# SPDX-License-Identifier: MIT

add_executable(heat_equation heat_equation.cpp)
target_link_libraries(heat_equation PUBLIC DDC::DDC)
add_executable(heat_equation_spectral heat_equation_spectral.cpp)
target_link_libraries(heat_equation_spectral PUBLIC DDC::DDC)
add_executable(game_of_life game_of_life.cpp)
target_link_libraries(game_of_life PUBLIC DDC::DDC)
