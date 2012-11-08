if(PROJECT_NAME)
  message(FATAL_ERROR "C++11 support has to be enabled before the first project definition")
else()
  set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_LIST_DIR}/CXX11Toolchain.cmake)

  include(${CMAKE_CURRENT_LIST_DIR}/CXX11Toolchain.cmake)

  if(CXX_LIB_ABI_ID)
    list(INSERT CMAKE_INCLUDE_PATH 0 /opt/${CXX_LIB_ABI_ID}/include)
    list(INSERT CMAKE_LIBRARY_PATH 0 /opt/${CXX_LIB_ABI_ID}/lib)

    set(CXX_LIB_ABI_ID ${CXX_LIB_ABI_ID} CACHE STRING
      "Identifies c++ standard libraries with incompatible ABIs")
    set(CXX_LIB_ABI_TAG "-${CXX_LIB_ABI_ID}")
    string(TOUPPER ${CXX_LIB_ABI_TAG} CXX_LIB_ABI_TAG_UPPER)
  endif()

  # TODO currently we do not actually check for working C++11 support
  set(CXX11_FOUND TRUE)
endif()
