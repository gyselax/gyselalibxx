
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
    gyselalibxx::data_types
    gyselalibxx::math_tools
    gyselalibxx::multipatch_data_types
    gyselalibxx::multipatch_utils
    gyselalibxx::utils
)
add_library(gyselalibxx::timestepper ALIAS timestepper)

install(TARGETS timestepper EXPORT gyselalibxxTargets FILE_SET HEADERS)
