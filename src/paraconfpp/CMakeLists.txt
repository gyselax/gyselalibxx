# SPDX-License-Identifier: MIT

add_library(paraconfpp INTERFACE)

target_compile_features("paraconfpp" INTERFACE cxx_std_20)

target_sources("paraconfpp"
    INTERFACE
        FILE_SET HEADERS
            FILES
                paraconfpp.hpp
)

target_link_libraries(paraconfpp INTERFACE
        paraconf::paraconf
)
add_library(gslx::paraconfpp ALIAS paraconfpp)

install(TARGETS paraconfpp EXPORT GyselalibxxTargets FILE_SET HEADERS)
