#ifndef VECRNG_R123_H
#define VECRNG_R123_H 1

/**
  R123: NameSpace for SIMD/SIMT Random124 classes based on Random123
 */

// namespace for Random-123 parameters

#include "RngDefs.h"
#include <array>
#include <stdint.h>

namespace vecRng {
inline namespace VECRNG_IMPL_NAMESPACE {

namespace R123 {

// typedef

template <typename BackendT, unsigned int N>
using array_t = typename BackendT::UInt32_v[N];

// Threefry parameters

VECRNG_GLOBAL unsigned int R_64x4[8][2] = {{14, 16}, {52, 57}, {23, 40},
                                           {5, 37},  {25, 33}, {46, 12},
                                           {58, 22}, {32, 32}};

VECRNG_GLOBAL unsigned int R_32x4[8][2] = {{10, 26}, {11, 21}, {13, 27},
                                           {23, 5},  {6, 20},  {17, 11},
                                           {25, 10}, {18, 20}};

VECRNG_GLOBAL long long THREEFRY_SKIP_AHEAD = 2 ^ 64;

// Philox parameters
VECRNG_GLOBAL uint32_t PHILOX_M4_W32_0 = (uint32_t)0xD2511F53;
VECRNG_GLOBAL uint32_t PHILOX_M4_W32_1 = (uint32_t)0xCD9E8D57;

VECRNG_GLOBAL uint32_t PHILOX_N4_W32_0 = (uint32_t)0x9E3779B9;
VECRNG_GLOBAL uint32_t PHILOX_N4_W32_1 = (uint32_t)0xBB67AE85;

VECRNG_GLOBAL long long PHILOX_SKIP_AHEAD = 2 ^ 64;

} // namespace R123

} // namespace VECRNG_IMPL_NAMESPACE
} // end namespace vecRng

#endif
