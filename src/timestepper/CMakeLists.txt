
add_library(timestepper INTERFACE)

target_compile_features("timestepper" INTERFACE cxx_std_20)

target_sources("timestepper"
    INTERFACE
        FILE_SET HEADERS
            FILES
                crank_nicolson.hpp
                euler.hpp
                itimestepper.hpp
                rk2.hpp
                rk3.hpp
                rk4.hpp
)

target_link_libraries(timestepper INTERFACE
    DDC::core
    gslx::data_types
    gslx::math_tools
    gslx::multipatch_data_types
    gslx::multipatch_utils
    gslx::utils
)
add_library(gslx::timestepper ALIAS timestepper)

install(TARGETS timestepper EXPORT GyselalibxxTargets FILE_SET HEADERS)
