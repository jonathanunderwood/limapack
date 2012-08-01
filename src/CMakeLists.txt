add_library(limapack STATIC
  memory.c
  config_gsl_complex.c  
  laser.c  
  laser_cfg_parse.c  
  laser_polzn.c  
  laser_type1.c
  jkmarray.c
  jmarray.c
  polarizability.c
  odesys.c
)
target_link_libraries(limapack config gsl m)

