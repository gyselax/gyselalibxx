cmake_minimum_required(VERSION 3.18)

add_library(baselib
		src/cartesiandistribution2d.cpp
		src/distributed2dfield.cpp
		src/simulation.cpp
)
target_include_directories(baselib PUBLIC include/)
target_link_libraries(baselib PUBLIC MPI::MPI_CXX std::mdspan)
add_library(GLCS2020_project::baselib ALIAS baselib)
