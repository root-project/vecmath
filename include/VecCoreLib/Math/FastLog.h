#ifndef VECCORELIB_MATH_FASTLOG_H
#define VECCORELIB_MATH_FASTLOG_H

#include "VecCoreLib/Private/vdt/log.h"
#include <VecCore/VecCore>

namespace vecMath {
template <typename R> R FastLog(R x) {
  static_assert(sizeof(R) == 0);
  return x;
} // Does not compile

template <> double FastLog(double x) { return vdt::fast_log(x); }

#ifdef VECCORE_ENABLE_VC
template <> Vc::double_v FastLog(Vc::double_v x) { return Vc::log(x); }
#endif
}

#endif // VECCORELIB_MATH_FASTLOG_H
