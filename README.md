# Intro
  The vecmath project provides a collection of vectorized based math utilities based on VecCore library.
  
# Requirements

+ VecCore(starting from cee2341)
+ Optional: Vc, UMESIMD as VecCore backends

# How to use.

## 1) Use VecMath installed package

### Install this library.

```bash
mkdir build/ && cd build/
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
	      -DVecCore_DIR="path to dir with VecCoreConfig.cmake"  #if VecCore is not installed inside prefix
make install
```

### Use it from your project
In your CMakeLists.txt:
```cmake
find_package(VecMath 0.1.0 REQUIRED COMPONENTS Vc) #or any other VecCore backend

#...

target_link_libraries(YOUR_TARGET VecMath::VecMath)
```

or if you can't use new cmake `target_link_libraries` use old interface
```cmake
include_directories(${VecMath_INCLUDE_DIRS})
add_definitions(${VecMath_DEFINITIONS})

target_link_libraries(YOUR_TARGET ${VecMath_LIBRARIES})
```

Pass this to cmake generation of your project to find VecMath(if it is not in CMAKE_INSTALL_PREFIX):
`-DVecMath_DIR="Path to dir where VecMathConfig.cmake is installed"`


## 2) Build VecMath in source as a git submodule

Create git submodule in YOU_PROJECT/VecMath pointing to this repository
```
cd YOUR_PROJECT
git submodule add https://gitlab.cern.ch/GeantV/VecMath.git
```

In your CMakeLists.txt:
```cmake
set(VecCore_BACKEND Vc) #or any other VecCore backend
add_subdirectory(VecMath)

#...

target_link_libraries(YOUR_TARGET VecMath)
```


## API
Inside you .cxx file 
```cpp
#include <VecMath/Math.h>
#include <VecMath/Rng.h>
```
