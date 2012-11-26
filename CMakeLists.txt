cmake_minimum_required(VERSION 2.6)

# Set to wherever your cmake modules are
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Make sure you have all required modules/libraries installed in the system
find_package(CXX11 REQUIRED)

project(probability_library)

if(INCLUDE_INSTALL_DIR)
else()
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

message(STATUS "Include Directory (INCLUDE_INSTALL_DIR):   ${INCLUDE_INSTALL_DIR}")

# Lets be pedantic
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")

# Testing support, recompile gtest due to C++11 
add_subdirectory(gtest-1.6.0)
enable_testing(true)

include_directories(src ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_subdirectory(src)

# Install directory
set(INCLUDE_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/include/prob"
    CACHE INTERNAL
    "Header file install directory"
  )

# Documentation generation
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

# Test binary
add_executable(test_random_variable test/Tests.cpp test/RandomVariableTest.cpp)
target_link_libraries(test_random_variable gtest gtest_main)
add_test(random_variable test_random_variable)

add_executable(test_splitter test/Tests.cpp test/SplitterTest.cpp)
target_link_libraries(test_splitter gtest gtest_main)
add_test(splitter test_splitter)

add_executable(test_distribution test/Tests.cpp test/DistributionTest.cpp)
target_link_libraries(test_distribution gtest gtest_main)
add_test(distribution test_distribution)

add_executable(test_algebra test/Tests.cpp test/AlgebraTest.cpp)
target_link_libraries(test_algebra gtest gtest_main)
add_test(algebra test_algebra)

add_executable(test_information test/Tests.cpp test/InformationTest.cpp)
target_link_libraries(test_information gtest gtest_main)
add_test(information test_information)