cmake_minimum_required(VERSION 3.6.0)

include(cmake/RecordCmdLine.cmake)

project(VecCoreLib VERSION 0.1.0)

add_library(VecCoreLib INTERFACE)
target_include_directories(VecCoreLib INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

option(BUILD_TESTING    "Build test binaries and create test target" OFF)
option(BUILD_BENCHMARKS "Build binaries for performance benchmarking")

if(NOT VecCore_BACKEND)
	set(VecCore_BACKEND Vc)
endif()


#Find VecCore with correct backend
find_package(VecCore 0.4.2 REQUIRED COMPONENTS ${VecCore_BACKEND})
target_include_directories(VecCoreLib INTERFACE ${VecCore_INCLUDE_DIRS})
target_link_libraries(VecCoreLib INTERFACE ${VecCore_LIBRARIES})
target_compile_definitions(VecCoreLib INTERFACE ${VecCore_DEFINITIONS})


#GetVDT
include(ExternalProject)
ExternalProject_Add(vdt
        GIT_REPOSITORY https://github.com/dpiparo/vdt.git

        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/builtin_vdt/install

        INSTALL_DIR ${CMAKE_BINARY_DIR}/builtin_vdt/install
        INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/builtin_vdt/include/VecCoreLib/Private/vdt
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/builtin_vdt/install/include/vdt ${CMAKE_BINARY_DIR}/builtin_vdt/include/VecCoreLib/Private/vdt
	)
add_dependencies(VecCoreLib PRIVATE vdt)
#when building include vdt from downloaded repo
target_include_directories(VecCoreLib INTERFACE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/builtin_vdt/include>)


#installation configuration
#Install vdt headers(without shared lib) to $INSTALL_DESTINATION/include/VecCoreLib/Private/vdt/*.h when VecCoreLib is installed
install(DIRECTORY ${CMAKE_BINARY_DIR}/builtin_vdt/include/VecCoreLib/Private DESTINATION include/VecCoreLib FILES_MATCHING PATTERN "*.h")
install(DIRECTORY include/ DESTINATION include)
install(TARGETS VecCoreLib EXPORT VecCoreLibConfig)
install(EXPORT VecCoreLibConfig
	DESTINATION lib/cmake/VecCoreLib
	)

#install version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file("VecCoreLibConfigVersion.cmake"
	VERSION ${VecCoreLib_VERSION}
	COMPATIBILITY SameMajorVersion
	)
install(FILES ${CMAKE_BINARY_DIR}/VecCoreLibConfigVersion.cmake
	DESTINATION lib/cmake/VecCoreLib
	)


if (BUILD_TESTING)
  add_subdirectory(test)
endif()


if (BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
