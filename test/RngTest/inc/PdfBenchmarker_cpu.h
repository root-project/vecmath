#ifndef PDF_BENCHMARKER_CPU_H
#define PDF_BENCHMARKER_CPU_H 1

#include "RngTest.h"

namespace vecrng {

typedef Real_t (*ScalarKernelFunc_t)(int nsample, double& result);

// Scalar

Real_t ScalarMRG32k3aExp(int nsample, double& result);
Real_t ScalarMRG32k3aNormal(int nsample, double& result);

 ScalarKernelFunc_t ScalarKernelFunc[] = {ScalarMRG32k3aExp,
                                         ScalarMRG32k3aNormal};

// Vector

typedef Real_t (*VectorKernelFunc_t)(int nsample, double& result);

Real_t VectorMRG32k3aExp(int nsample, double& result);
Real_t VectorMRG32k3aNormal(int nsample, double& result);

VectorKernelFunc_t VectorKernelFunc[] = {VectorMRG32k3aExp,
                                         VectorMRG32k3aNormal};

} // end namespace vecrng

#endif
