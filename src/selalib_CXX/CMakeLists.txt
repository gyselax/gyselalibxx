ADD_SUBDIRECTORY( matrix )
ADD_SUBDIRECTORY( splines )

# Parallelization tools
IF( MPI_ENABLED )
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   #   ADD_SUBDIRECTORY( parallelization )
ENDIF( MPI_ENABLED )
