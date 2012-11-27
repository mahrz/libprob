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