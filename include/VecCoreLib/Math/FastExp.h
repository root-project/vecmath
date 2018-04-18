#ifndef VECCORELIB_MATH_FASTEXP_H
#define VECCORELIB_MATH_FASTEXP_H

#include "VecCoreLib/Private/vdt/exp.h"
#include <VecCore/VecCore>

namespace vecMath {
template <typename R> R FastExp(R x) {
  static_assert(sizeof(R) == 0);
  return x;
} // Does not compile

template <> double FastExp(double x) { return vdt::fast_exp(x); }

#ifdef VECCORE_ENABLE_VC
template <> Vc::double_v FastExp(Vc::double_v x) { return Vc::exp(x); }
#endif
}

#endif // VECCORELIB_MATH_FASTEXP_H
