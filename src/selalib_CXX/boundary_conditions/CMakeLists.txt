SET( BC_SRC
  sll_m_boundary_condition_descriptors.F90
  )

add_library( boundary_conditions OBJECT ${BC_SRC} )

set_property(TARGET boundary_conditions PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/boundary_conditions_FMOD")

target_include_directories( boundary_conditions
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}"
    INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/boundary_conditions_FMOD"
)
add_library( selalib::boundary_conditions ALIAS boundary_conditions )
