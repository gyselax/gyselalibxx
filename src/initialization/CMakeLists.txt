add_library(initialization STATIC
    singlemodeperturbinitialization.cpp
)

target_compile_features(initialization
    PUBLIC
        cxx_std_17
)

target_include_directories(initialization
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(initialization
    PUBLIC
        DDC::DDC
        vcx::speciesinfo
        vcx::geometry
)

add_library(vcx::initialization ALIAS initialization)
