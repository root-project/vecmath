# Intro
  The VecCoreLib project gathers VecCore-based utilities developed by the GeantV R&D, but having no GeantV dependencies. This is a temporary solution for avoiding developing only in GeantV general code usable in other areas than simulation, such as RND or fast math.
  
# Requirements

+ VecCore(starting from bb702eb8b90098237a43b415ef69475c697bee26)
+ Optional: Vc, UMESIMD as a VecCore backends

# How to use.

## 1) Use VecCoreLib installed package

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
find_package(VecCoreLib 0.1.0 REQUIRED COMPONENTS Vc) #or any other VecCore backend

#...

target_link_libraries(YOUR_TARGET VecCoreLib::VecCoreLib)
```

Pass this to cmake generation of your project to find VecCoreLib(if it is not in CMAKE_INSTALL_PREFIX):
`-DVecCoreLib_DIR="Path to dir where VecCoreLibConfig.cmake is installed"`


## 2) Build VecCoreLib in source as a git submodule

Create git submodule in YOU_PROJECT/VecCoreLib pointing to this repository
```
cd YOUR_PROJECT
git submodule add https://gitlab.cern.ch/GeantV/VecCoreLib.git
```

In your CMakeLists.txt:
```cmake
set(VecCore_BACKEND Vc) #or any other VecCore backend
add_subdirectory(VecCoreLib)

#...

target_link_libraries(YOUR_TARGET VecCoreLib)
```


## API
Inside you .cxx file 
```cpp
#include <VecCoreLib/Math.h>
#include <VecCoreLib/Rng.h>
```
