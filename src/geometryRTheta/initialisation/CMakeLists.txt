
# DIOCOTRON SIMULATION ------------------------------------------
add_library(diocotron_init_eq STATIC
    diocotron_initialisation_equilibrium.cpp
)

target_compile_features("diocotron_init_eq" PUBLIC cxx_std_20)

target_sources("diocotron_init_eq"
    INTERFACE
        FILE_SET HEADERS
            FILES
                diocotron_initialisation_equilibrium.hpp
)

target_link_libraries(diocotron_init_eq
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::geometry_RTheta
        gslx::utils
)

add_library("gslx::diocotron_init_eq" ALIAS "diocotron_init_eq")

install(TARGETS diocotron_init_eq EXPORT GyselalibxxTargets FILE_SET HEADERS)



# VORTEX MERGER SIMULATION --------------------------------------
add_library(vortex_merger_init_eq INTERFACE)

target_compile_features("vortex_merger_init_eq" INTERFACE cxx_std_20)

target_sources("vortex_merger_init_eq"
    INTERFACE
        FILE_SET HEADERS
            FILES
                vortex_merger_equilibrium.hpp
                vortex_merger_initialisation.hpp
)
target_link_libraries(vortex_merger_init_eq
    INTERFACE
        DDC::core
        gslx::geometry_RTheta
        gslx::utils
        gslx::poisson_RTheta
)

add_library("gslx::vortex_merger_init_eq" ALIAS "vortex_merger_init_eq")

install(TARGETS vortex_merger_init_eq EXPORT GyselalibxxTargets FILE_SET HEADERS)
