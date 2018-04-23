cmake_minimum_required(VERSION 3.6.0)

include(cmake/RecordCmdLine.cmake)

option(BUILD_TESTING    "Build test binaries and create test target" OFF)
option(BUILD_BENCHMARKS "Build binaries for performance benchmarking")
set(VecCore_VERSION 0.4.2)

project(VecCoreLib VERSION 0.1.0)

add_library(VecCoreLib INTERFACE)
target_include_directories(VecCoreLib INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

#Find VecCore with correct backend
find_package(VecCore ${VecCore_VERSION} REQUIRED COMPONENTS ${VecCore_BACKEND})
target_link_libraries(VecCoreLib INTERFACE VecCore::VecCore) 


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
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/include/VecCoreLib/Private/vdt
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/install/include/vdt ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/include/VecCoreLib/Private/vdt
	)
")
execute_process(COMMAND ${CMAKE_COMMAND} . -G ${CMAKE_GENERATOR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/)
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/)
execute_process(COMMAND ${CMAKE_COMMAND} --build . --target install
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/)

#when building include vdt from downloaded repo
target_include_directories(VecCoreLib INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/include>)


#installation configuration
#Install vdt headers(without shared lib) to $INSTALL_DESTINATION/include/VecCoreLib/Private/vdt/*.h when VecCoreLib is installed
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builtin_vdt/include/VecCoreLib/Private DESTINATION include/VecCoreLib FILES_MATCHING PATTERN "*.h")
install(DIRECTORY include/ DESTINATION include)

include(CMakePackageConfigHelpers)

#exported targets file
install(TARGETS VecCoreLib EXPORT VecCoreLibTargets
  INCLUDES DESTINATION include
)
export(EXPORT VecCoreLibTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/VecCoreLibTargets.cmake"
  NAMESPACE VecCoreLib::
)

#package config file
set(VecCoreLib_INSTALL_INCLUDEDIR include)
configure_package_config_file(cmake/VecCoreLibConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/VecCoreLibConfig.cmake
  INSTALL_DESTINATION lib/cmake/VecCoreLib
  PATH_VARS VecCoreLib_INSTALL_INCLUDEDIR
)

#version file
write_basic_package_version_file("VecCoreLibConfigVersion.cmake"
	VERSION ${VecCoreLib_VERSION}
	COMPATIBILITY SameMajorVersion
	)

#install package
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/VecCoreLibConfigVersion.cmake
${CMAKE_CURRENT_BINARY_DIR}/VecCoreLibConfig.cmake
${CMAKE_CURRENT_BINARY_DIR}/VecCoreLibTargets.cmake
	DESTINATION lib/cmake/VecCoreLib
	)


if (BUILD_TESTING)
  add_subdirectory(test)
endif()


if (BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
