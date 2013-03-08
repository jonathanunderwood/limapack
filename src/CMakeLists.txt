# Make sure hdf5 is properly used
include_directories(${HDF5_INCLUDE_DIRS})
add_definitions(${HDF5_DEFINITIONS})

# TODO - eventually these should becompiled into a shared library.
#add_library(limalib OBJECT
set (common_files
  memory.c
  config_gsl_complex.c  
  hdf5_gsl_complex.c
  laser.c  
  laser_polzn.c  
  laser_type1.c
  jmarray_int.c
  jmarray_double.c
  jkmarray_int.c
  jkmarray_double.c
  polarizability.c
  odesys.c
  molecule.c
  molecule_linear.c
  dmtxel.c
  dcmsq.c
  slurp.c
  asymrot_eigsys.c
  molecule_asymrot.c
)

# add_library(limalib OBJECT
#   ${common_files}
#   )

#target_link_libraries(limalib ${LIBCONFIG_LIBRARIES} gsl m)

add_executable(lima
  lima.c
  ${common_files}
#  $<TARGET_OBJECTS:limalib>
)
target_link_libraries(lima ${LIBCONFIG_LIBRARIES} ${HDF5_LIBRARIES} ${GSL_LIBRARIES} ${MATH_LIBRARIES})

if (${CMAKE_USE_PTHREADS_INIT})
  target_link_libraries (lima ${CMAKE_THREAD_LIBS_INIT})
#  set_target_properties (limalib lima PROPERTIES
  set_target_properties (lima PROPERTIES
    COMPILE_DEFINITIONS "BUILD_WITH_PTHREADS"
    )
endif ()


## Set up options for compiling against MPI if present
if (MPI_C_FOUND)
  # add_library(limalib_mpi OBJECT
  #   ${common_files}
  #   )

  add_executable(lima_mpi
    lima_mpi.c
    #    $<TARGET_OBJECTS:limalib>
    ${common_files}
    )

  target_link_libraries(lima_mpi ${LIBCONFIG_LIBRARIES} ${HDF5_LIBRARIES} ${GSL_LIBRARIES} ${MATH_LIBRARIES})
  # set_target_properties(limalib lima lima_mpi PROPERTIES
  #   COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}"
  #   COMPILE_DEFINITIONS "BUILD_WITH_MPI"
  #   LINK_FLAGS "${MPI_C_LINK_FLAGS}"
  #   )
#  set_property (TARGET limalib_mpi lima_mpi
  set_property (TARGET lima_mpi
  	        APPEND PROPERTY COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}"
		)
  set_property (TARGET lima_mpi
#  set_property (TARGET limalib_mpi lima_mpi
  	        APPEND PROPERTY COMPILE_DEFINITIONS "BUILD_WITH_MPI"
		)
  set_property (TARGET lima_mpi
#  set_property (TARGET limalib_mpi lima_mpi
  	        APPEND PROPERTY LINK_FLAGS "${MPI_C_LINK_FLAGS}"
		)
 
  include_directories(${MPI_C_INCLUDE_PATH})
#  target_link_libraries(lima ${MPI_C_LIBRARIES})
  target_link_libraries(lima_mpi ${MPI_C_LIBRARIES})
  # get_cmake_property(_variableNames VARIABLES)
  # foreach (_variableName ${_variableNames})
  #   message(STATUS "${_variableName}=${${_variableName}}")
  # endforeach()
endif()	


