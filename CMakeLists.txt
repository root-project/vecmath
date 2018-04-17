################################################################################
# VecCoreLib Utility Library                                                   #
################################################################################

cmake_minimum_required(VERSION 3.6.0)

include(cmake/RecordCmdLine.cmake)

project(VecCoreLib LANGUAGES C CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11 CACHE STRING "C++ ISO Standard")
set(CMAKE_CXX_STANDARD_REQUIRED True)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Select build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

option(VECRNG  "Enable support for VecRng" ON)

option(BUILD_BENCHMARKS "Build binaries for performance benchmarking")
option(BUILD_TESTING    "Build test binaries and create test target")

# Let's see if we can find VecCore
find_package(VecCore ${VecCore_VERSION} QUIET COMPONENTS ${VecCore_COMPONENTS})
if (NOT VecCore_FOUND)
  set(BUILTIN_VECCORE ON)
else()
  message("-- Found VecCore at: ${VecCore_ROOT_DIR}")
endif()

include_directories(include ${PROJECT_BINARY_DIR}/include)

configure_file(include/Config.in include/Config.h)

set(VecCoreLib_INSTALL_INCLUDEDIR "include")
set(VecCoreLib_INSTALL_CMAKEDIR "share/${PROJECT_NAME}/cmake")

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/VecCoreLibConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/VecCoreLibConfig.cmake
  INSTALL_DESTINATION "${VecCoreLib_INSTALL_CMAKEDIR}"
  PATH_VARS VecCoreLib_INSTALL_INCLUDEDIR
  )
configure_file(cmake/VecCoreLibConfigVersion.cmake.in cmake/VecCoreLibConfigVersion.cmake @ONLY)

install(FILES "${PROJECT_BINARY_DIR}/include/Config.h" DESTINATION ${VecCoreLib_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
#install(DIRECTORY include/${PROJECT_NAME} DESTINATION ${VecCoreLib_INSTALL_INCLUDEDIR})
install(DIRECTORY ${PROJECT_BINARY_DIR}/cmake/
  DESTINATION ${VecCoreLib_INSTALL_CMAKEDIR} FILES_MATCHING PATTERN *.cmake)

if (VECRNG)
  install(DIRECTORY include/rng DESTINATION ${VecCoreLib_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    PATTERN "test" EXCLUDE)
endif()

enable_testing()

if (BUILD_TESTING)
  add_subdirectory(test)
endif()

if (BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
