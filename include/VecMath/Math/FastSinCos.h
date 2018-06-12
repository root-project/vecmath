#ifndef VECMATH_MATH_FASTSINCOS_H
#define VECMATH_MATH_FASTSINCOS_H

#include "VecMath/Private/vdt/sincos.h"
#include <VecCore/VecCore>

namespace vecMath {
template <typename R> inline void FastSinCos(R, R &, R &) {
  static_assert(sizeof(R) == 0, "Backend is not supported");
} // Does not compile

template <> inline void FastSinCos(double x, double &s, double &c) {
  return vdt::fast_sincos(x, s, c);
}

#ifdef VECCORE_ENABLE_VC
template <>
inline void FastSinCos(Vc::double_v x, Vc::double_v &s, Vc::double_v &c) {
  return Vc::sincos(x, &s, &c);
}
#endif

template <typename R> inline R FastSin(R x) {
  R s, c;
  FastSinCos(x, s, c);
  return s;
}

template <typename R> inline R FastCos(R x) {
  R s, c;
  FastSinCos(x, s, c);
  return c;
}
}

#endif // VECMATH_MATH_FASTSINCOS_H
