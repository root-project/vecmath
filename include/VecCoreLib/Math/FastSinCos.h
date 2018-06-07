#ifndef VECCORELIB_MATH_FASTSINCOS_H
#define VECCORELIB_MATH_FASTSINCOS_H

#include "VecCoreLib/Private/vdt/sincos.h"
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

#ifdef VECCORE_ENABLE_AGNER
template <> inline void FastSinCos(vcl::Vec4d x, vcl::Vec4d &s, vcl::Vec4d &c) {
  s = vcl::sincos(&c, x);
}
template <> inline void FastSinCos(vcl::Vec8f x, vcl::Vec8f &s, vcl::Vec8f &c) {
  s = vcl::sincos(&c, x);
}
template <> inline void FastSinCos(vcl::Vec8d x, vcl::Vec8d &s, vcl::Vec8d &c) {
  s = vcl::sincos(&c, x);
}
template <>
inline void FastSinCos(vcl::Vec16f x, vcl::Vec16f &s, vcl::Vec16f &c) {
  s = vcl::sincos(&c, x);
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
} // namespace vecMath

#endif // VECCORELIB_MATH_FASTSINCOS_H
