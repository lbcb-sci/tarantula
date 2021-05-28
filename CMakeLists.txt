cmake_minimum_required(VERSION 3.11)

project(tarantula VERSION 0.0.1
                  LANGUAGES CXX
                  DESCRIPTION "Tarantula is a Hi-C scaffolder.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(tarantula_main_project ON)
endif ()
option(tarantula_build_tests "Build unit tests" ${tarantula_main_project})

find_package(bioparser 3.0.13 QUIET)
find_package(cereal 1.3.0 QUIET)
find_package(ram 2.0.0 QUIET)

if (NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.13)

  FetchContent_GetProperties(bioparser)
  if (NOT bioparser_POPULATED)
    FetchContent_Populate(bioparser)
    add_subdirectory(
      ${bioparser_SOURCE_DIR}
      ${bioparser_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

if (NOT cereal_FOUND)
  FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal
    GIT_TAG v1.3.0)

  FetchContent_GetProperties(cereal)
  if (NOT cereal_POPULATED)
    FetchContent_Populate(cereal)
    add_subdirectory(
      ${cereal_SOURCE_DIR}
      ${cereal_BINARY_DIR}
      EXCLUDE_FROM_ALL)
    add_library(cereal::cereal ALIAS cereal)
  endif ()
endif ()

if (NOT ram_FOUND)
  FetchContent_Declare(
    ram
    GIT_REPOSITORY https://github.com/lbcb-sci/ram
    GIT_TAG 7597511)

  FetchContent_GetProperties(ram)
  if (NOT ram_POPULATED)
    FetchContent_Populate(ram)
    add_subdirectory(
      ${ram_SOURCE_DIR}
      ${ram_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

if (tarantula_build_tests)
  find_package(GTest 1.10.0 QUIET)
  if (NOT GTest_FOUND)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG release-1.10.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()
endif ()

add_library(tarantula
  src/graph.cpp)

target_link_libraries(tarantula
  cereal::cereal
  ram::ram)

target_include_directories(tarantula PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

add_executable(tarantula_exe
  src/main.cpp)

target_link_libraries(tarantula_exe
  tarantula
  bioparser::bioparser)

target_compile_definitions(tarantula_exe PRIVATE VERSION="${PROJECT_VERSION}")
set_property(TARGET tarantula_exe PROPERTY OUTPUT_NAME tarantula)

install(TARGETS tarantula_exe DESTINATION ${CMAKE_INSTALL_BINDIR})

if (tarantula_build_tests)
  add_executable(tarantula_test
    test/tarantula_test.cpp)

  target_link_libraries(tarantula_test
    tarantula
    bioparser::bioparser
    GTest::Main)

  target_compile_definitions(tarantula_test
    PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/")
endif ()
