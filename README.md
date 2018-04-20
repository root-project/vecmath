# Intro
  The VecCoreLib project gathers VecCore-based utilities developed by the GeantV R&D, but having no GeantV dependencies. This is a temporary solution for avoiding developing only in GeantV general code usable in other areas than simulation, such as RND or fast math.
  
# Requirements

VecCore(starting from bb702eb8b90098237a43b415ef69475c697bee26)
Vc, UMESIMD as a VecCore backends

# How to use.

## 1) Build VecCoreLib in source

Copy VecCoreLib to your project

In your CMakeLists.txt:
```cmake
set(VecCore_BACKEND Vc) #optional Vc is default backend
add_subdirectory(VecCoreLib)

...

target_link_libraries(YOUR_TARGET VecCoreLib)
```

## 2) Use VecCoreLib installed package

### Install this library.

```bash
mkdir build/ && cd build/
cmake ../ -DVecCore_BACKEND=Vc -DCMAKE_INSTALL_PREFIX=$PREFIX 
          -DVc_DIR="path to dir with VcConfig.cmake"  #if Vc is not installed inside prefix
          -DVecCore_DIR="path to dir with VecCoreConfig.cmake"  #if VecCore is not installed inside prefix
make install
```

Pass this to cmake generation of your project to find VecCoreLib(if it is not in CMAKE_INSTALL_PREFIX):
`-DVecCoreLib_DIR="Path to dir where VecCoreLibConfig.cmake is installed"`


### Use it as a package
In your CMakeLists.txt:
```cmake
find_package(VecCoreLib 0.1.0 REQUIRED COMPONENTS Vc) #or any other VecCore backends

#...

target_link_libraries(YOUR_TARGET VecCoreLib::VecCoreLib)
```

Pass this to cmake generation of your project to find VecCoreLib(if it is not in CMAKE_INSTALL_PREFIX):
`-DVecCoreLib_DIR="Path to dir where VecCoreLibConfig.cmake is installed"`

## API
Inside you .cxx file 
```cpp
#include <VecCoreLib/Math.h>
#include <VecCoreLib/Rng.h>
```
