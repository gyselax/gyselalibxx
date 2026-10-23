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
        gyselalibxx::data_types
        gyselalibxx::geometry_RTheta
        gyselalibxx::speciesinfo
        gyselalibxx::utils
)

add_library("gyselalibxx::poisson_RTheta" ALIAS "poisson_RTheta")

install(TARGETS poisson_RTheta EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryRTheta/poisson)
