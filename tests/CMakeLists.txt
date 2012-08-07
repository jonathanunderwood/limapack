include_directories ("${PROJECT_SOURCE_DIR}/src")

add_executable (test_laser test_laser.c)
target_link_libraries(test_laser limapack)
add_test(NAME test_laser
  COMMAND test_laser)

add_executable (test_jmarray test_jmarray.c)
target_link_libraries(test_jmarray limapack)
add_test(NAME test_jmarray
  COMMAND test_jmarray)

add_executable (test_jkmarray test_jkmarray.c)
target_link_libraries(test_jkmarray limapack)
add_test(NAME test_jkmarray
  COMMAND test_jkmarray)

add_executable (test_polarizability test_polarizability.c)
target_link_libraries(test_polarizability limapack)
add_test(NAME test_polarizability
  COMMAND test_polarizability)

add_executable (test_ode test_ode.c)
target_link_libraries(test_ode limapack)
add_test(NAME test_ode
  COMMAND test_ode)
