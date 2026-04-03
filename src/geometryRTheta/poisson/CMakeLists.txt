# SPDX-License-Identifier: MIT


add_library(poisson_RTheta INTERFACE)

target_compile_features("poisson_RTheta" INTERFACE cxx_std_20)

target_sources("poisson_RTheta"
    INTERFACE
        FILE_SET HEADERS
            FILES
                iqnsolver.hpp
                poisson_like_rhs_function.hpp
)

target_link_libraries(poisson_RTheta INTERFACE
        DDC::core
        gslx::data_types
        gslx::geometry_RTheta
        gslx::speciesinfo
        gslx::utils
)

add_library("gslx::poisson_RTheta" ALIAS "poisson_RTheta")

install(TARGETS poisson_RTheta EXPORT GyselalibxxTargets FILE_SET HEADERS)
