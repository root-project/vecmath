#ifndef VECMATH_MATH_FASTEXP_H
#define VECMATH_MATH_FASTEXP_H

#include "vdt/exp.h"
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
} // namespace vecMath

#endif // VECMATH_MATH_FASTEXP_H
