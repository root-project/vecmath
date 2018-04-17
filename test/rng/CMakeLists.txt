#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.1.0)

project(vecrng)
set(RNGTEST ON)

enable_language(CXX)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 11 CACHE STRING "C++ ISO Standard")
set(CMAKE_CXX_STANDARD_REQUIRED True)

###############################################################################
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)

# Setting the C++ compiler to use
set("${CMAKE_CXX_COMPILER_ID}" TRUE CACHE STRING "C++ Compiler")
if (NOT GNU AND NOT Clang AND NOT AppleClang AND NOT Intel)
  message(WARNING "Unsupported compiler.")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

include(MacroUtilities)
###############################################################################

# Configuration options

option(RNGTEST_CUDA "Enable compilation for CUDA." OFF)

set(RNGTEST_ARCH native sse2 ssse3 sse4.1 sse4.2 aes avx avx2 mic mic_avx512)

enum_option(RNGTEST_VECTOR
      DOC "Vector instruction set to be used"
      TYPE STRING
      VALUES ${RNGTEST_ARCH}
      CASE_INSENSITIVE
    )

if ("${BACKEND}" MATCHES "vc|VC")
  set(BACKEND "Vc")
endif()
if ("${BACKEND}" MATCHES "scalar|SCALAR")
  set(BACKEND "Scalar")
endif()

set(RNGTEST_BACKEND Scalar UMESIMD Vc)

enum_option(BACKEND
            DOC "Backend infrastructure to be used"
            TYPE STRING
            VALUES ${RNGTEST_BACKEND}
            CASE_SENSITIVE
           )

################################################################################

if (RNGTEST_CUDA)
  set(VecCore_COMPONENTS ${VecCore_COMPONENTS} "CUDA")
  set(CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "CUDA separable compilation on.")
endif()

# Determine backend

set(${BACKEND} True)
string(TOUPPER "${BACKEND}" _BACKEND_UP)
string(TOLOWER "${BACKEND}" _BACKEND_LO)

if("${BACKEND}" MATCHES "Scalar")
  set(RNGTEST_SCALAR True)
  message(STATUS "Configuring with Scalar backend")
elseif ("${BACKEND}" MATCHES "Vc|UMESIMD")
  set(RNGTEST_${_BACKEND_UP} True)
  set(VecCore_COMPONENTS ${VecCore_COMPONENTS} ${BACKEND})
  message(STATUS "Configuring with ${BACKEND} backend")
else()
  message(FATAL_ERROR "Unknown backend: '${BACKEND}'.")
endif()


# Find VecCore with selected components turned on (CUDA and backend)
find_package(VecCore REQUIRED COMPONENTS ${VecCore_COMPONENTS})

include_directories(SYSTEM ${VecCore_INCLUDE_DIRS})
set(RNGTEST_DEFINITIONS "${RNGTEST_DEFINITIONS} ${VecCore_DEFINITIONS}")
set(RNGTEST_DEFINITIONS "${RNGTEST_DEFINITIONS} -DRNGTEST_${_BACKEND_UP}")
set(RNGTEST_EXTERNAL_INCLUDES "${RNGTEST_EXTERNAL_INCLUDES};${VecCore_INCLUDE_DIRS}")
set(RNGTEST_LIBRARIES_EXTERNAL ${RNGTEST_LIBRARIES_EXTERNAL} ${VecCore_LIBRARIES})

################################################################################

# Set Default compiler flags for each build type

if (GNU)
  set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -Wall -fPIC")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -ggdb -O0")
  else()
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -finline-limit=10000000 -ffast-math -ftree-vectorize")
  endif()
  set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -fabi-version=0")
  set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -fabi-version=0")
elseif(Intel)
  # CMake fails to add -std option for icc, so add it manually
  if("${CMAKE_VERSION}" VERSION_LESS "3.6")
    if(CMAKE_CXX_STANDARD STREQUAL 11)
      list(APPEND CMAKE_CXX_FLAGS -std=c++11)
    elseif(CMAKE_CXX_STANDARD STREQUAL 14)
      list(APPEND CMAKE_CXX_FLAGS -std=c++14)
    else()
      message(FATAL_ERROR "Unsupported C++ standard requested")
    endif()
  endif()
  set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -Wall -fPIC -diag-disable 3438")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -g -O0")
  else()
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -fno-alias")
  endif()
elseif(Clang OR AppleClang)
  if (APPLE)
    # prefer LLVM's standard C++ library on Mac OS X
    set (RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -stdlib=libc++")
  endif()
  set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -Wall -fPIC")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS}  -ggdb -O0")
  else()
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -ffast-math -ftree-vectorize")
  endif()
else()
   message(WARNING "No compiler type matched")
endif()

if (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  set(RNGTEST_NVCC_COMPILATION_FLAGS -g -G ${RNGTEST_NVCC_COMPILATION_FLAGS})
endif()
if (CMAKE_BUILD_TYPE MATCHES MinSizeRel)
  set(RNGTEST_NVCC_COMPILATION_FLAGS -use_fast_math ${RNGTEST_NVCC_COMPILATION_FLAGS})
endif()

string(TOLOWER ${RNGTEST_VECTOR} _arch_lo)
string(TOUPPER ${RNGTEST_VECTOR} _arch_up)

if (${_arch_lo} MATCHES mic_avx512)
  if (Intel)
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -xMIC_AVX512")
    set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -xMIC_AVX512")
  else()
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -march=knl")
    set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -march=knl")
  endif()
elseif (${_arch_lo} MATCHES mic)
  if (Intel)
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -m${_arch_lo}")
    set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -m${_arch_lo}")
    set(MIC TRUE)
  else()
   message(FATAL_ERROR "MIC native or offload compilation requires Intel Compiler.")
  endif()
elseif (${_arch_lo} MATCHES native)
  if (Intel)
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -xHost")
    set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -xHost")
  else()
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -march=${_arch_lo}")
    set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -march=${_arch_lo}")
  endif()
else()
  if (Intel)
    if (${_arch_lo} MATCHES avx2)
      set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -xCORE-${_arch_up}")
      set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -xCORE-${_arch_up}")
    else()
      set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -x${_arch_up}")
      set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -x${_arch_up}")
    endif()
  else()
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -m${_arch_lo}")
    set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -m${_arch_lo}")
  endif()
endif()
message(STATUS "Compiling for ${_arch_up} SIMD architecture")

################################################################################

# Set C++ compiler flags. Will also be forwarded by CUDA when compiling C++.
if (RNGTEST_CUDA)
  set(RNGTEST_COMPILATION_FLAGS "${RNGTEST_COMPILATION_FLAGS} -DRNGTEST_CUDA")

  set(CUDA_ARCH 30 CACHE STRING "CUDA device architecture.")
  set(CUDA_ARCH "-arch=sm_${CUDA_ARCH}")
  set(RNGTEST_DEFINITIONS "${RNGTEST_DEFINITIONS} -DRNGTEST_ENABLE_CUDA")
  set(RNGTEST_NVCC_COMPILATION_FLAGS ${RNGTEST_NVCC_COMPILATION_FLAGS}
      -std=c++11
      -Xcompiler -Wno-unused-function
      -Xcudafe "--diag_suppress=code_is_unreachable"
      -Xcudafe "--diag_suppress=initialization_not_reachable")
  set(RNGTEST_NVCC_CONFIGS
        DEBUG -g -G
        RELEASE  -O3 -use_fast_math
        RELWITHDEBINFO -g -G -O3 -use_fast_math
        MINSIZEREL -use_fast_math )
  if (NOT (CMAKE_VERSION LESS 3.3))
     cmake_policy(SET CMP0058 OLD)
  endif()
endif()

if (Clang OR AppleClang)
  set(RNGTEST_LIBRARIES_EXTERNAL ${RNGTEST_LIBRARIES_EXTERNAL} dl)
endif()

# Enable pretty diagnostics coloring if gcc 4.9+ is being used
if (GNU)
  # gcc needs -ldl to be added explicitly
  set(RNGTEST_LIBRARIES_EXTERNAL ${RNGTEST_LIBRARIES_EXTERNAL} dl)
  set(RNGTEST_ERROR_LIMIT 20 CACHE STRING
      "Limit number of errors output by diagnostics.")
  set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -fmax-errors=${RNGTEST_ERROR_LIMIT}")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)
  if (GCC_VERSION VERSION_GREATER 4.9 OR GCC_VERSION VERSION_EQUAL 4.9)
    set(RNGTEST_CXX_FLAGS "${RNGTEST_CXX_FLAGS} -fdiagnostics-color=auto")
  endif()
endif()

################################################################################

# see if we have CCACHE ( to avoid recompilation on branch switches )
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  message(STATUS "found ccache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

if (RNGTEST_CUDA)
  find_package(CUDA REQUIRED)
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${RNGTEST_NVCC_COMPILATION_FLAGS})
endif()

# Pass flags to compilers
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}} ${RNGTEST_CXX_FLAGS} ${RNGTEST_DEFINITIONS} ${RNGTEST_COMPILATION_FLAGS}")
string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

message(STATUS "Compiling with C++ flags: ${CMAKE_CXX_FLAGS}")
if (RNGTEST_CUDA)
  message(STATUS "Compiling with NVCC flags: ${CUDA_NVCC_FLAGS}")
endif()

################################################################################

include_directories(${PROJECT_SOURCE_DIR} inc)

set(SRC_CPP_RELATIVE
  src/RngBenchmarker.cc
  src/RngBenchmarker_cpu.cc

  src/PdfBenchmarker.cc
  src/PdfBenchmarker_cpu.cc
)

if (RNGTEST_CUDA)
#-------------------------------------------------------------------------------
# Copy source files to .cu files in order to compile them with NVCC
#
  foreach(SRC_FILE ${SRC_CPP})
    get_filename_component(SRC_FILENAME ${SRC_FILE} NAME_WE)

    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      COMMAND ${CMAKE_COMMAND} -E copy ${SRC_FILE}
      ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      DEPENDS ${SRC_FILE}
      )
    
    set(
      SRC_CUDA ${SRC_CUDA}
      ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      )
  endforeach()

  set(SRC_CUDA ${SRC_CUDA}
    src/RngBenchmarker.cu
    src/RngBenchmarker_gpu.cu

    src/PdfBenchmarker.cu
    src/PdfBenchmarker_gpu.cu
  )
endif()

foreach(SRC ${SRC_CPP_RELATIVE})
  set(SRC_CPP ${SRC_CPP} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC})
endforeach()

set(RNGTEST_SRCS
  ${RNGTEST_SRCS}
  ${RNGTEST_CUDA_SRCS}
)

foreach(SRC ${RNGTEST_SRCS})
  set(SRC_CPP ${SRC_CPP} ${CMAKE_SOURCE_DIR}/${SRC})
endforeach()

################################################################################
# Build libraries

add_library(rngtest ${SRC_CPP})
target_link_libraries(rngtest ${RNGTEST_LIBRARIES}
                      ${RNGTEST_LIBRARIES_EXTERNAL})

set(RNGTEST_LIBRARIES ${RNGTEST_LIBRARIES} rngtest)

if (APPLE)
   set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
else()
  set(RNGTEST_LIBRARIES_EXTERNAL ${RNGTEST_LIBRARIES_EXTERNAL} rt)
endif()

if (BUILD_SHARED_LIBS)
    set(RNGTEST_LIBNAME "librngtest${CMAKE_SHARED_LIBRARY_SUFFIX}")
else()
    set(RNGTEST_LIBNAME "librngtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif()

# build the CUDA version of the library
if (RNGTEST_CUDA)
  cuda_add_library(
    rngtestcuda
    ${SRC_CUDA}
    SHARED
    OPTIONS ${CUDA_ARCH}
    ${RNGTEST_NVCC_CONFIGS}
  )
  target_link_libraries(rngtestcuda ${RNGTEST_LIBRARIES}
                        ${RNGTEST_LIBRARIES_EXTERNAL})
  set(RNGTEST_LIBRARIES ${RNGTEST_LIBRARIES} rngtestcuda)
endif()

################################################################################
# Define executables
set(TEST_EXECUTABLES
  ${PROJECT_SOURCE_DIR}/src/RngBenchmark.cc
  ${PROJECT_SOURCE_DIR}/src/PdfBenchmark.cc
)

# Build executables
set(RNGTEST_EXECUTABLES)

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    target_link_libraries(${TARGET_NAME} ${RNGTEST_LIBRARIES}
                          ${RNGTEST_LIBRARIES_EXTERNAL})
    set(RNGTEST_EXECUTABLES ${RNGTEST_EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

build_executables("${TEST_EXECUTABLES}")
if (RNGTEST_CUDA)
  build_executables("${TEST_EXECUTABLES_CUDA}")
endif()

################################################################################
