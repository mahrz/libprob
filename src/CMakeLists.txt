install(FILES
  prob
  DESTINATION ${INCLUDE_INSTALL_DIR} COMPONENT Devel
  )

file(GLOB prob_directory_files_to_install "*.hpp")

install(FILES
  ${prob_directory_files_to_install}
  DESTINATION ${INCLUDE_INSTALL_DIR} COMPONENT Devel
  )

file(GLOB prob_directory_files_to_install "Util/*.hpp")

install(FILES
  ${prob_directory_files_to_install}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Util COMPONENT Devel
  )

file(GLOB prob_directory_files_to_install "InformationTheory/*.hpp")

install(FILES
  ${prob_directory_files_to_install}
  DESTINATION ${INCLUDE_INSTALL_DIR}/InformationTheory COMPONENT Devel
  )

add_executable(prob_test main.cpp)
target_link_libraries(prob_test nlopt_cxx)