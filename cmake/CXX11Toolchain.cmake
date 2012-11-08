# TODO find a nicer way to implement all of this

# Newer CMake versions need to write to CMAKE_PLATFORM_INFO_DIR to
# determine the C++ compiler
set(_platform_dir ${CMAKE_PLATFORM_INFO_DIR})
set(CMAKE_PLATFORM_INFO_DIR "${CMAKE_CURRENT_BINARY_DIR}")
include(CMakeDetermineCXXCompiler)
set(CMAKE_PLATFORM_INFO_DIR ${_platform_dir})

# Append -std=c++11 to COMPILER_ARG1. This ensures that the flag is
# passed even when general metadata about the compiler is queried,
# which is important for the eclipse and codeblocks generators and
# maybe other places
if("${CMAKE_CXX_COMPILER_ID}" MATCHES Clang OR
    "${CMAKE_CXX_COMPILER_ID}" MATCHES GNU OR
    "${CMAKE_CXX_COMPILER_ID}" MATCHES Intel) # TODO intel compiler untested
  # TODO also support gnu++11?
  if(NOT "${CMAKE_CXX_COMPILER_ARG1}" MATCHES "-std=c\\+\\+11")
    set(CMAKE_CXX_COMPILER_ARG1 "${CMAKE_CXX_COMPILER_ARG1} -std=c++11")
  endif()
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES Clang AND
    "${CMAKE_CXX_COMPILER_ARG1}" MATCHES "-stdlib=libc\\+\\+")
  # TODO this is only tested with arch linux's libc++-svn aur package
  # and assume libstdc++ provides the runtime
  if(NOT "${CMAKE_EXE_LINKER_FLAGS_INIT}" MATCHES "-Wl,-lstdc\\+\\+")
    set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -Wl,-lstdc++")
  endif()
  if(NOT CXX_LIB_ABI_ID)
    set(CXX_LIB_ABI_ID libc++-c++11)
  endif()
endif()

if(NOT CXX_LIB_ABI_ID)
  set(CXX_LIB_ABI_ID c++11)
endif()

