
add_library("multipatch_advection" INTERFACE)

target_compile_features("multipatch_advection" INTERFACE cxx_std_20)

target_sources("multipatch_advection"
    INTERFACE
        FILE_SET HEADERS
            FILES
                bsl_advection.hpp
                spline_foot_finder.hpp
)

target_link_libraries("multipatch_advection" INTERFACE
    DDC::core
    DDC::splines

    gyselalibxx::data_types
    gyselalibxx::utils
)
add_library("gyselalibxx::multipatch_advection" ALIAS "multipatch_advection")

install(TARGETS multipatch_advection EXPORT gyselalibxxTargets FILE_SET HEADERS)
