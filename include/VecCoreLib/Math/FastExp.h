#ifndef VECCORELIB_MATH_FASTEXP_H
#define VECCORELIB_MATH_FASTEXP_H

#include "VecCoreLib/Private/vdt/exp.h"
#include <VecCore/VecCore>

namespace vecMath {
template <typename R> inline R FastExp(R x) {
  static_assert(sizeof(R) == 0, "Backend is not supported");
  return x;
} // Does not compile

template <> inline double FastExp(double x) { return vdt::fast_exp(x); }

#ifdef VECCORE_ENABLE_VC
template <> inline Vc::double_v FastExp(Vc::double_v x) { return Vc::exp(x); }
#endif

#ifdef VECCORE_ENABLE_AGNER
template <> inline vcl::Vec4d FastExp(vcl::Vec4d x) { return vcl::exp(x); }
template <> inline vcl::Vec8f FastExp(vcl::Vec8f x) { return vcl::exp(x); }
template <> inline vcl::Vec8d FastExp(vcl::Vec8d x) { return vcl::exp(x); }
template <> inline vcl::Vec16f FastExp(vcl::Vec16f x) { return vcl::exp(x); }
#endif
} // namespace vecMath

#endif // VECCORELIB_MATH_FASTEXP_H
