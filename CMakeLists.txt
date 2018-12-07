cmake_minimum_required(VERSION 3.6.0)

include(cmake/RecordCmdLine.cmake)

option(BUILD_TESTING    "Build test binaries and create test target"  OFF)
option(BUILD_BENCHMARKS "Build binaries for performance benchmarking" OFF)
set(VecCore_VERSION 0.4.2)

if(BUILD_TESTING)
  message("Vc backend selected for testing")
  set(VecCore_BACKEND Vc)
endif()

project(VecMath VERSION 0.1.0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

add_library(VecMath INTERFACE)
target_include_directories(VecMath INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

#Find VecCore with correct backend
find_package(VecCore ${VecCore_VERSION} REQUIRED COMPONENTS ${VecCore_BACKEND})
target_link_libraries(VecMath INTERFACE VecCore::VecCore) 

#Find vdt if available
find_package(Vdt)
if(VDT_FOUND)
  include_directories(${VDT_INCLUDE_DIR})
else()

  #GetVDT on cmake generation phase to prevent rebuilding when using this lib from inside your project
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/CMakeLists.txt
  "
cmake_minimum_required(VERSION 3.6.0)
include(ExternalProject)
ExternalProject_Add(vdt
        GIT_REPOSITORY https://github.com/dpiparo/vdt.git

        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/install

        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/install
        INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/include/VecMath/Private/vdt
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/install/include/vdt ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/include/VecMath/Private/vdt
	)
  "
  )
  execute_process(COMMAND ${CMAKE_COMMAND} . -G ${CMAKE_GENERATOR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/)
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/)
  execute_process(COMMAND ${CMAKE_COMMAND} --build . --target install
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/)

  #when building include vdt from downloaded repo
  target_include_directories(VecMath INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/include>)

  #installation configuration
  #Install vdt headers(without shared lib) to $INSTALL_DESTINATION/include/VecMath/Private/vdt/*.h when VecMath is installed
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/include/VecMath/Private DESTINATION include/VecMath FILES_MATCHING PATTERN "*.h")

endif()

install(DIRECTORY include/ DESTINATION include)

include(CMakePackageConfigHelpers)

#exported targets file
install(TARGETS VecMath EXPORT VecMathTargets
  INCLUDES DESTINATION include
)
export(EXPORT VecMathTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/VecMathTargets.cmake"
  NAMESPACE VecMath::
)

#package config file
set(VecMath_INSTALL_INCLUDEDIR include)
configure_package_config_file(cmake/VecMathConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/VecMathConfig.cmake
  INSTALL_DESTINATION lib/cmake/VecMath
  PATH_VARS VecMath_INSTALL_INCLUDEDIR
)

#version file
write_basic_package_version_file("VecMathConfigVersion.cmake"
	VERSION ${VecMath_VERSION}
	COMPATIBILITY SameMajorVersion
	)

#install package
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/VecMathConfigVersion.cmake
${CMAKE_CURRENT_BINARY_DIR}/VecMathConfig.cmake
${CMAKE_CURRENT_BINARY_DIR}/VecMathTargets.cmake
	DESTINATION lib/cmake/VecMath
	)


if (BUILD_TESTING)
  add_subdirectory(test)
endif()


if (BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
