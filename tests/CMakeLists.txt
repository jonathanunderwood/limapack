include_directories ("${PROJECT_SOURCE_DIR}/src")

add_executable (test_laser test_laser.c)
target_link_libraries(test_laser limapack config gsl m)
add_test(NAME test_laser
  COMMAND test_laser)
