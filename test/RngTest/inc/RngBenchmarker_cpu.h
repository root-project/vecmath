#ifndef RNG_BENCHMARKER_CPU_H
#define RNG_BENCHMARKER_CPU_H 1

#include "RngTest.h"
#ifdef RNGTEST_MKL
#include "mkl_vsl.h"
#endif

namespace vecrng {

typedef Real_t (*ScalarKernelFunc_t)(int nsample, double& result);

// Scalar

Real_t ScalarMRG32k3a(int nsample, double& result);
Real_t ScalarThreefry(int nsample, double& result);
Real_t ScalarPhilox(int nsample, double& result);

ScalarKernelFunc_t ScalarKernelFunc[] = {ScalarMRG32k3a, ScalarThreefry, ScalarPhilox};

// Vector

typedef Real_t (*VectorKernelFunc_t)(int nsample, double& result);

Real_t VectorMRG32k3a(int nsample, double& result);
Real_t VectorThreefry(int nsample, double& result);
Real_t VectorPhilox(int nsample, double& result);

Real_t VectorJoiningMRG32k3a(int nsample, double& result);
// Real_t VectorJoiningRng_MRG32k3a(int nsample, double& result);
// Real_t VectorJoiningRng_ThreeFry(int nsample, double& result);

VectorKernelFunc_t VectorKernelFunc[] = {VectorMRG32k3a,
                                         VectorJoiningMRG32k3a,
                                         VectorThreefry,
                                         VectorPhilox};

using FuncAndName = std::pair< VectorKernelFunc_t, std::string>;

// std::vector< FuncAndName >
FuncAndName VectorKnlFuncAndName[] = {
   { VectorMRG32k3a,            "VectorMRG32k3a" },
   { VectorJoiningMRG32k3a,     "VectorJoiningMRG32k3a" },
// #ifdef NEW_TEST_JOINING_PROXY   
//    { VectorJoiningRng_MRG32k3a, "VectorJoiningRng_MRG32k3a" },
// //  { VectorJoiningRng_ThreeFry, "VectorJoiningRng_ThreeFry" },
// #endif   
   { VectorThreefry,            "VectorThreefry" },
   { VectorPhilox,              "VectorPhilox" }
 };
constexpr int numVecKnlFuncs= 4; // Number of elements of VectorKnlFuncAndName[];

// for (auto funcAndName : VectorKnlFuncAndName )

// Scalar-States (comparision to GPU code)

typedef Real_t (*StateKernelFunc_t)(int nsample, double& result);

Real_t StateMRG32k3a(int nsample, double& result);
Real_t StateThreefry(int nsample, double& result);
Real_t StatePhilox(int nsample, double& result);

StateKernelFunc_t StateKernelFunc[] = {StateMRG32k3a, StateThreefry, StatePhilox};

#ifdef RNGTEST_MKL
// MKL/VSL (comparision to Intel Vector Libraries)
typedef Real_t (*VSLKernelFunc_t)(int nsample, double& result);

Real_t VSLRngTest(int nsample, double& result, int rngType);

Real_t VSLMRG32k3a(int nsample, double& result) { return VSLRngTest(nsample, result, VSL_BRNG_MRG32K3A); }
Real_t VSLMCG59(int nsample, double& result)    { return VSLRngTest(nsample, result, VSL_BRNG_MCG59); }
Real_t VSLMT19937(int nsample, double& result)  { return VSLRngTest(nsample, result, VSL_BRNG_MT19937); }
Real_t VSLSMT19937(int nsample, double& result) { return VSLRngTest(nsample, result, VSL_BRNG_SFMT19937); }
Real_t VSLGFSR250(int nsample, double& result)  { return VSLRngTest(nsample, result, VSL_BRNG_R250); }
Real_t VSLSOBOL32(int nsample, double& result)  { return VSLRngTest(nsample, result, VSL_BRNG_SOBOL); }

//OMP
//Real_t OPENMPMRG(int nsample, double& result);

VSLKernelFunc_t VSLKernelFunc[] = {VSLMRG32k3a, VSLMCG59, VSLMT19937, VSLSMT19937, VSLGFSR250, VSLSOBOL32};
#endif

} // end namespace vecrng

#endif
