#ifndef PDF_BENCHMARKER_GPU_H
#define PDF_BENCHMARKER_GPU_H 1

#include "VecMath/Rng/MRG32k3a.h"

#include <curand_kernel.h>

namespace vecrng {

  // Cuda Wrapper
  void CudaMRG32k3aGauss(vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t *devStates, double *result, 
		    int nsample, int blocksPerGrid, int threadsPerBlock);

  // Curand Test 
  void CurandMRG32k3aGauss(curandStateMRG32k3a *devStates, double *result,
		      int nsample, int blocksPerGrid, int threadsPerBlock);

  void curand_setup_gpu_gauss(curandStateMRG32k3a* devStates, unsigned long seed,  
                        int blocksPerGrid, int threadsPerBlock);

} // end namespace vecrng

#endif
