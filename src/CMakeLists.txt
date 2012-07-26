add_library(limapack STATIC
  memory.c
  config_gsl_complex.c  
  laser.c  
  laser_cfg_parse.c  
  laser_polzn.c  
  laser_type1.c
)
target_link_libraries(limapack libconfig gsl m)


add_executable(test_laser test_laser.c)
target_link_libraries(test_laser limapack)
