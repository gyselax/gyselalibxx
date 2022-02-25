# SPDX-License-Identifier: MIT

add_library(poisson STATIC
    chargedensitycalculator.cpp
    electricfield.cpp
    nullpoissonsolver.cpp
)

target_compile_features(poisson
    PUBLIC
        cxx_std_17
)

target_include_directories(poisson
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(poisson
    PUBLIC
        DDC::DDC
        sll::splines
        vcx::speciesinfo
        vcx::geometry
)

if(VOICEXX_ENABLE_PERIODIC_RDIMX)
    target_sources(poisson
        PRIVATE
            fftpoissonsolver.cpp
            femperiodicpoissonsolver.cpp
    )
    target_link_libraries(poisson
        PUBLIC
            vcx::fft)
else()
    target_sources(poisson
        PRIVATE
            femnonperiodicpoissonsolver.cpp
    )
endif()

add_library(vcx::poisson ALIAS poisson)
