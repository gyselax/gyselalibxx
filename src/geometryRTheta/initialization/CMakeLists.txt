

add_library(diocotron_init_eq STATIC
    diocotron_initialization_equilibrium.cpp
)

target_compile_features(diocotron_init_eq
    PUBLIC
        cxx_std_17
)

target_include_directories(diocotron_init_eq
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(diocotron_init_eq
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper        
        Eigen3::Eigen
        
        sll::splines
        gslx::geometry_RTheta
)

add_library("gslx::diocotron_init_eq" ALIAS "diocotron_init_eq")