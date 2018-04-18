# Intro
  The VecCoreLib project gathers VecCore-based utilities developed by the GeantV R&D, but having no GeantV dependencies. This is a temporary solution for avoiding developing only in GeantV general code usable in other areas than simulation, such as RND or fast math.
  

# How to use.

## Install this library.

```bash
mkdir build/ && cd build/
cmake ../ -DVecCore_BACKEND=Vc -DCMAKE_INSTALL_PREFIX=$PREFIX 
          -DVc_DIR="path to dir with VcConfig.cmake"  #if Vc is not installed inside prefix
          -DVecCore_DIR="path to dir with VecCoreConfig.cmake"  #if VecCore is not installed inside prefix
make install
```

## Use from your project.

In your CMakeLists.txt:
```cmake
find_package(VecCoreLib REQUIRED)

#...

# after this you will be able to use VecCoreLib includes and VecCore with correct 
# backend(backend is select when you install VecCoreLib)

target_link_libraries(YOUR_TARGET VecCoreLib)
```

Inside you .cxx file 
```cpp
#include <VecCoreLib/Math.h>
#include <VecCoreLib/Rng.h>
```

To find project pass to your cmake generation only:
`-DVecCoreLib_DIR="Path to dir where VecCoreLibConfig.cmake is installed"`
