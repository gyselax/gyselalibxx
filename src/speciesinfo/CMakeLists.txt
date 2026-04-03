# SPDX-License-Identifier: MIT

add_library(speciesinfo STATIC
    species_init.cpp
)

target_compile_features("speciesinfo" PUBLIC cxx_std_20)

target_sources("speciesinfo"
    INTERFACE
        FILE_SET HEADERS
            FILES
                species_info.hpp
                species_init.hpp
)

target_link_libraries(speciesinfo
    DDC::core
    gslx::io
    gslx::paraconfpp
    gslx::utils
)
add_library(gslx::speciesinfo ALIAS speciesinfo)

install(TARGETS speciesinfo EXPORT GyselalibxxTargets FILE_SET HEADERS)

add_library(moments INTERFACE)

target_sources("moments"
    INTERFACE
        FILE_SET HEADERS
            FILES
                moments.hpp
)

target_link_libraries(moments INTERFACE
    DDC::core)
add_library(gslx::moments ALIAS moments)
install(TARGETS moments EXPORT GyselalibxxTargets FILE_SET HEADERS)
