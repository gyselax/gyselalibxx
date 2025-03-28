
# DIOCOTRON SIMULATION ------------------------------------------
add_library(diocotron_init_eq STATIC
    diocotron_initialisation_equilibrium.cpp
)

target_include_directories(diocotron_init_eq
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(diocotron_init_eq
    PUBLIC
        DDC::core
        DDC::pdi        
        Eigen3::Eigen
        
        gslx::geometry_RTheta
        gslx::utils
)

add_library("gslx::diocotron_init_eq" ALIAS "diocotron_init_eq")



# VORTEX MERGER SIMULATION --------------------------------------
add_library(vortex_merger_init_eq INTERFACE)
target_include_directories(vortex_merger_init_eq
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(vortex_merger_init_eq
    INTERFACE
        DDC::core
        #DDC::pdi        
        #Eigen3::Eigen
        
        gslx::geometry_RTheta
        gslx::utils
        gslx::poisson_RTheta
)

add_library("gslx::vortex_merger_init_eq" ALIAS "vortex_merger_init_eq")
