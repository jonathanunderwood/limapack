add_library(limapack STATIC
  memory.c
  config_gsl_complex.c  
  laser.c  
  laser_polzn.c  
  laser_type1.c
  jkmarray.c
  jmarray.c
  polarizability.c
  odesys.c
  molecule.c
  molecule_linear.c
  dmtxel.c
)

target_link_libraries(limapack ${LIBCONFIG_LIBRARIES} gsl m)

add_executable(lima
	lima.c
)

target_link_libraries(lima limapack ${LIBCONFIG_LIBRARIES} gsl m)
